#include "subsumed_search.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace ckm {

SubsumedSearch::SubsumedSearch(const IncrementalGraph& graph)
    : graph_(graph) {}

void SubsumedSearch::validateBatchEdges() const {
    for (EdgeId edgeId : batchEdges_) {
        if (edgeId < 0 || edgeId >= graph_.edgeCount()) {
            throw std::invalid_argument("batch contains an invalid edge ID");
        }
        // Inserted edges are in the graph by now; deleted ones are gone.
        if (graph_.isActive(edgeId) == recovering_) {
            throw std::invalid_argument(
                recovering_
                    ? "recovery requires the deleted edges to be out of the graph"
                    : "subsumed search requires active batch edges");
        }
    }
}

void SubsumedSearch::prepareBatch(
    std::span<const EdgeId> batch,
    SubsumedSearchStats& stats) {
    ++stats.batches;
    batchEdges_.assign(batch.begin(), batch.end());
    validateBatchEdges();
    batchPrepared_ = true;
    newCliqueHashes_.clear();
    probeBuffer_.clear();
    probeArena_.clear();
}

void SubsumedSearch::consumeNewClique(
    const Clique& canonicalNewClique,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    if (!batchPrepared_) {
        throw std::logic_error("the subsumed batch must be prepared first");
    }
    if (!std::ranges::is_sorted(canonicalNewClique) ||
        std::ranges::adjacent_find(canonicalNewClique) !=
            canonicalNewClique.end()) {
        throw std::invalid_argument("streamed new clique must be canonical");
    }
    consumeNewCliqueKernelEnum(canonicalNewClique, oldCliques, stats, output);
}

void SubsumedSearch::finishBatch(
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    flushDeferredProbes(oldCliques, stats, output);
    if (bitmapExtractionPrepared_) {
        releaseBitmapExtraction();
    }
    batchPrepared_ = false;
}

void SubsumedSearch::flushDeferredProbes(
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    if (probeBuffer_.empty()) {
        return;
    }
    stats.maximumDeferredProbes = std::max(
        stats.maximumDeferredProbes, probeBuffer_.size());
    // Nothing in this loop depends on an earlier probe, so the whole buffer
    // can be kept in flight.  Inline probing could hide one cache miss at a
    // time; this fills the machine's memory queue.
    //
    // One probe is not one access, though.  A probe that misses is settled by
    // the table, but a probe that hits then reads its record and then that
    // record's vertices, and neither address is known until the read before it
    // has landed -- so prefetching a probe's table slot hides one miss of
    // three.  Splitting the buffer into stages breaks that chain: every probe
    // is resolved against the table first, which is where the misses can all
    // be issued at once, and only the few that survive go on to the record and
    // the arena, each stage again prefetched ahead of itself.
    constexpr std::size_t lookahead = 16U;
    const std::size_t count = probeBuffer_.size();
    resolvedProbes_.clear();
    for (std::size_t index = 0; index < std::min(lookahead, count); ++index) {
        oldCliques.prefetch(probeBuffer_[index].setHash);
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (index + lookahead < count) {
            oldCliques.prefetch(probeBuffer_[index + lookahead].setHash);
        }
        ++stats.indexLookups;
        const std::size_t slot =
            oldCliques.firstTagMatch(probeBuffer_[index].setHash);
        if (slot == MaximalCliqueIndex::noSlot) {
            continue;
        }
        const std::uint32_t record = oldCliques.recordAtSlot(slot);
        oldCliques.prefetchRecord(record);
        // The record is captured here rather than re-read from the slot
        // later: by then an earlier duplicate may have tombstoned the slot,
        // and a tombstone is not a record id.  Records themselves are stable
        // for the whole flush, since only an insertion can renumber them.
        resolvedProbes_.push_back(
            ResolvedProbe{slot, static_cast<std::uint32_t>(index), record});
    }

    // The record's own line is already on its way for the first few; from
    // there on each iteration starts the arena fetch of one further ahead.
    constexpr std::size_t verifyLookahead = 8U;
    const std::size_t resolved = resolvedProbes_.size();
    for (std::size_t index = 0; index < std::min(verifyLookahead, resolved);
         ++index) {
        oldCliques.prefetchRecordVertices(resolvedProbes_[index].record);
    }
    for (std::size_t index = 0; index < resolved; ++index) {
        if (index + verifyLookahead < resolved) {
            oldCliques.prefetchRecordVertices(
                resolvedProbes_[index + verifyLookahead].record);
        }
        const DeferredProbe& probe =
            probeBuffer_[resolvedProbes_[index].probeIndex];
        std::optional<Clique> erased = oldCliques.extractMaskedAtSlot(
            resolvedProbes_[index].slot,
            std::span<const Vertex>(
                probeArena_.data() + probe.offset, probe.size),
            std::span<const std::uint64_t>(&probe.mask, 1U),
            probe.memberCount,
            probe.setHash);
        if (erased.has_value()) {
            ++stats.emittedCliques;
            output.emitCanonical(std::move(*erased));
        }
    }
    resolvedProbes_.clear();
    probeBuffer_.clear();
    probeArena_.clear();
}

// One candidate, either buffered for the batch's flush or probed on the spot.
// Buffering copies the producing clique into the arena at most once, however
// many of its candidates end up deferred.
void SubsumedSearch::submitProbe(
    std::span<const Vertex> clique,
    std::span<const std::uint64_t> maskWords,
    std::size_t memberCount,
    std::uint64_t setHash,
    std::uint32_t& arenaOffset,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    if (recovering_) {
        acceptRecovered(clique, maskWords, memberCount, setHash, stats, output);
        return;
    }
    // A buffered probe carries one mask word, which is what keeps it at 24
    // bytes; the handful of cliques wider than that probe inline.
    if (maskWords.size() == 1U) {
        if (arenaOffset == noArenaOffset) {
            arenaOffset = static_cast<std::uint32_t>(probeArena_.size());
            probeArena_.insert(
                probeArena_.end(), clique.begin(), clique.end());
        }
        probeBuffer_.push_back(DeferredProbe{
            setHash,
            maskWords[0],
            arenaOffset,
            static_cast<std::uint16_t>(clique.size()),
            static_cast<std::uint16_t>(memberCount)});
        ++stats.deferredProbes;
        return;
    }
    ++stats.indexLookups;
    ++stats.inlineProbes;
    std::optional<Clique> erased = oldCliques.extractMaskedHashed(
        clique, maskWords, memberCount, setHash);
    if (erased.has_value()) {
        ++stats.emittedCliques;
        output.emitCanonical(std::move(*erased));
    }
}

// Both candidates of a one-edge kernel are C minus one endpoint: the two
// masks, the member count and the two hashes are immediate, and the two
// probes are independent random accesses that can be started together.  The
// two hashes follow from the clique's own by one XOR each.
void SubsumedSearch::probeOneEdgeKernel(
    const Clique& clique,
    std::size_t uIndex,
    std::size_t vIndex,
    std::uint64_t setAccumulator,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    ++stats.oneEdgeKernels;
    const std::size_t memberCount = clique.size() - 1U;
    if (!candidateSizePossible(oldCliques, memberCount)) {
        stats.sizeGateSkippedProbes += 2;
        return;
    }
    const std::size_t dropIndex[2] = {uIndex, vIndex};
    const std::uint64_t setHashes[2] = {
        finishSetHash(setAccumulator ^ vertexSetWord(clique[uIndex]), memberCount),
        finishSetHash(setAccumulator ^ vertexSetWord(clique[vIndex]), memberCount)};
    // Started now, so that the table lines are on their way by the time the
    // probes -- or the batch's flush -- reach them.  A recovery probes no
    // table there.
    if (!recovering_) {
        oldCliques.prefetch(setHashes[0]);
        oldCliques.prefetch(setHashes[1]);
    }
    const std::uint64_t full = clique.size() == 64U
        ? std::numeric_limits<std::uint64_t>::max()
        : (std::uint64_t{1} << clique.size()) - 1U;
    std::uint32_t arenaOffset = noArenaOffset;
    for (std::size_t choice = 0; choice < 2U; ++choice) {
        const std::uint64_t bits =
            full & ~(std::uint64_t{1} << dropIndex[choice]);
        submitProbe(
            clique,
            std::span<const std::uint64_t>(&bits, 1U),
            memberCount,
            setHashes[choice],
            arenaOffset,
            oldCliques,
            stats,
            output);
    }
}

std::vector<Clique> SubsumedSearch::enumerateBatch(
    std::span<const EdgeId> batch,
    std::span<const Clique> newCliques,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats) {
    prepareBatch(batch, stats);

    CliqueBatchOutput output(true);
    // Every clique in the batch is its own heap allocation, so the very first
    // thing done with it -- the pass that looks up its vertices -- misses.
    // The vector of cliques is walked in order, so those misses are known far
    // enough ahead to be started now.
    constexpr std::size_t cliqueLookahead = 16U;
    for (std::size_t index = 0;
         index < std::min(cliqueLookahead, newCliques.size()); ++index) {
        __builtin_prefetch(newCliques[index].data(), 0, 1);
    }
    for (std::size_t index = 0; index < newCliques.size(); ++index) {
        if (index + cliqueLookahead < newCliques.size()) {
            __builtin_prefetch(newCliques[index + cliqueLookahead].data(), 0, 1);
        }
        const Clique& rawClique = newCliques[index];
        if (std::ranges::is_sorted(rawClique) &&
            std::ranges::adjacent_find(rawClique) == rawClique.end()) {
            consumeNewCliqueKernelEnum(rawClique, oldCliques, stats, output);
            continue;
        }
        const Clique clique = canonicalClique(rawClique);
        consumeNewCliqueKernelEnum(clique, oldCliques, stats, output);
    }

    finishBatch(oldCliques, stats, output);
    // Handed back in the order the probes found them, and not checked for
    // repeats.  No caller of this timed path needs the order (the binary sorts
    // by its own output key, and only when it writes anything), and no clique
    // can be emitted twice -- a hit erases it, so a later probe for the same
    // set only meets the tombstone.  The unit tests sort and check the output
    // themselves before comparing it.
    return std::move(output).takeCliques();
}

// ------------------------------------------------------------- deletion

std::vector<Clique> SubsumedSearch::recoverBatch(
    std::span<const EdgeId> deletedBatch,
    std::span<const Clique> destroyed,
    MaximalCliqueIndex& cliques,
    SubsumedSearchStats& stats) {
    recovering_ = true;
    filterGathered_ = false;
    acceptedRecords_.clear();
    acceptedVertices_.clear();
    if (++acceptedEpoch_ == 0U) {
        std::ranges::fill(acceptedSlots_, 0U);
        acceptedEpoch_ = 1U;
    }
    if (acceptedSlots_.empty()) {
        acceptedSlots_.assign(1024U, 0U);
    }
    struct Reset {
        bool& recovering;
        bool& knownMaximal;
        std::vector<GatheredCandidate>& gathered;
        std::vector<Vertex>& gatheredVertices;
        ~Reset() {
            recovering = false;
            knownMaximal = false;
            gathered.clear();
            gatheredVertices.clear();
        }
    } reset{recovering_, knownMaximal_, gathered_, gatheredVertices_};
    prepareBatch(deletedBatch, stats);
    CliqueBatchOutput output(true);
    // A destroyed clique's buffer, its members' list headers, then the lists'
    // first entries, fetched ahead of it, in three steps that each wait for
    // the one before; its candidates are checked on those lists, or
    // gathered by their lengths.  (With the same for the deferred candidates
    // checked one by one: 0.93 of wikipedia's KernelEnum recovery, 0.99
    // of youtube's and sx-stackoverflow's.)
    constexpr std::size_t cliqueAhead = 12U;
    constexpr std::size_t headerAhead = 8U;
    constexpr std::size_t listAhead = 4U;
    constexpr std::size_t prefetchedSize = 64U;
    const std::size_t destroyedCount = destroyed.size();
    for (std::size_t index = 0U; index < destroyedCount; ++index) {
        if (index + cliqueAhead < destroyedCount) {
            __builtin_prefetch(destroyed[index + cliqueAhead].data(), 0, 1);
        }
        if (index + headerAhead < destroyedCount &&
            destroyed[index + headerAhead].size() <= prefetchedSize) {
            for (const Vertex vertex : destroyed[index + headerAhead]) {
                graph_.prefetchListHeader(vertex);
            }
        }
        if (index + listAhead < destroyedCount &&
            destroyed[index + listAhead].size() <= prefetchedSize) {
            for (const Vertex vertex : destroyed[index + listAhead]) {
                graph_.prefetchListStart(vertex);
            }
        }
        const Clique& rawClique = destroyed[index];
        if (std::ranges::is_sorted(rawClique) &&
            std::ranges::adjacent_find(rawClique) == rawClique.end()) {
            consumedClique_ = rawClique;
            consumeNewCliqueKernelEnum(rawClique, cliques, stats, output);
            continue;
        }
        const Clique clique = canonicalClique(rawClique);
        consumedClique_ = clique;
        consumeNewCliqueKernelEnum(clique, cliques, stats, output);
    }
    consumedClique_ = {};
    checkGatheredCandidates(stats, output);
    insertAcceptedInBatch(cliques);
    // Nothing was held back, so this only releases the batch structures.
    finishBatch(cliques, stats, output);
    return std::move(output).takeCliques();
}

namespace {

[[nodiscard]] bool adjacentInList(
    std::span<const AdjacencyEntry> list,
    Vertex vertex) noexcept {
    const auto found = std::lower_bound(
        list.begin(), list.end(), vertex,
        [](const AdjacencyEntry& entry, Vertex value) {
            return entry.neighbor < value;
        });
    return found != list.end() && found->neighbor == vertex;
}

// A member's row over a centre's neighbourhood is read by looking the
// member's list up in the centre's bits, an entry at a time, unless the
// member's list is so much the longer that searching it for each of the
// centre's vertices reads less (`groupRow`).
[[nodiscard]] bool rowBySearch(std::size_t centreDegree, std::size_t degree) noexcept {
    return 2U * centreDegree * (static_cast<std::size_t>(std::bit_width(degree)) + 1U) <
        degree;
}

// What a row costs to read when deciding which groups go on rows, in the
// entries a check one by one reads: the least of a merge of the two lists,
// a search of the member's for each of the centre's vertices, and a lookup
// of the member's in the centre's numbering, counted three.  (Rows were
// read the cheapest of these three ways when the rule was timed group by
// group; reading every row by lookup since, at an entry or so each, and
// counting it so, put three times the groups on rows for no gain:
// sx-stackoverflow 0.975, youtube 1.008, wikipedia 1.142.)
[[nodiscard]] long long rowReadCost(std::size_t centreDegree, std::size_t degree) noexcept {
    const auto merge = static_cast<long long>(centreDegree + degree);
    const auto search = static_cast<long long>(
        centreDegree * (static_cast<std::size_t>(std::bit_width(degree)) + 1U));
    const auto lookup = static_cast<long long>(3U * degree);
    return std::min({merge, search, lookup});
}

// Moves `value` to a fresh epoch, clearing `epochs` when the counter wraps.
void advanceEpoch(VertexScratch& epochs, std::uint32_t& value) noexcept {
    if (++value == 0U) {
        epochs.clear();
        value = 1U;
    }
}

// A candidate waits for the end of the batch when every member has at
// least this degree: a check from a member of lower degree reads few
// entries, and gathering it would cost about as much as the check (64 since
// the checks got faster; it was 32).
constexpr int deferDegree = 64;

// A group of fewer deferred candidates is checked one by one without
// weighing rows: so few rarely share their rows enough to pay for reading
// them, and weighing them costs about what it could save.
constexpr std::uint32_t minRowGroup = 8U;

// The most vertices the deferred candidates may hold before they are
// checked: 256 MB.  Without it coauthor's window of batch 200000 gathered
// 5.2 billion in one batch, more than 32-bit offsets reach.
constexpr std::size_t gatherLimit = std::size_t{1} << 26U;

// KernelEnum's search over a destroyed clique goes under extension masks once
// the lists its candidates' checks would scan add up to this many times the
// shortest list of the clique's universal part, which the masks scan whole.
// (Masking every such clique took 2.43 times as long over the paper's
// deletion batches.)
constexpr std::size_t maskFactor = 8U;

} // namespace

void VertexScratch::ensure(std::size_t size) {
    if (size <= size_) {
        return;
    }
    auto* grown = static_cast<std::uint32_t*>(std::calloc(size, sizeof(std::uint32_t)));
    if (grown == nullptr) {
        throw std::bad_alloc();
    }
    if (size_ != 0U) {
        std::memcpy(grown, data_.get(), size_ * sizeof(std::uint32_t));
    }
    data_.reset(grown);
    size_ = size;
}

void VertexScratch::clear() noexcept {
    if (size_ != 0U) {
        std::memset(data_.get(), 0, size_ * sizeof(std::uint32_t));
    }
}

// `gatheredKeys_` in order by their high halves, the keys of one high half
// in the order they came: two 16-bit radix passes over that half.  A pass
// costs its 65,536 counters whatever the count, which a few keys do not
// repay; a comparison sort of the whole keys then gives the same order, the
// low halves being positions that come in increasing order.
void SubsumedSearch::sortKeysByHighHalf() {
    if (gatheredKeys_.size() < 4096U) {
        std::ranges::sort(gatheredKeys_);
        return;
    }
    radixScratch_.resize(gatheredKeys_.size());
    for (const unsigned shift : {32U, 48U}) {
        radixCounts_.assign(std::size_t{1} << 16U, 0U);
        for (const std::uint64_t key : gatheredKeys_) {
            ++radixCounts_[static_cast<std::size_t>((key >> shift) & 0xFFFFU)];
        }
        std::uint32_t sum = 0U;
        for (std::uint32_t& bucket : radixCounts_) {
            const std::uint32_t here = bucket;
            bucket = sum;
            sum += here;
        }
        for (const std::uint64_t key : gatheredKeys_) {
            radixScratch_[radixCounts_[static_cast<std::size_t>((key >> shift) & 0xFFFFU)]++] = key;
        }
        gatheredKeys_.swap(radixScratch_);
    }
}

void SubsumedSearch::checkGatheredCandidates(
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    if (gathered_.empty()) {
        return;
    }

    // Equal candidates -- the same set from different destroyed cliques --
    // share a hash.  Only equal hashes need to be adjacent, so the keys are
    // the hash's low half above the position, put in order by two 16-bit
    // radix passes over that half; within a run of one low half (almost
    // always a single hash) the whole hash, size and vertices settle which
    // repeat an earlier one, and the first of each set is checked.  (A
    // comparison sort of the positions by their candidates' hashes missed
    // the cache on most comparisons, 17% of flickr's recovery, and of
    // (hash, position) pairs still took 8%.)
    const std::size_t count = gathered_.size();
    gatheredKeys_.resize(count);
    for (std::size_t index = 0; index < count; ++index) {
        gatheredKeys_[index] =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(gathered_[index].setHash)) << 32U) |
            static_cast<std::uint64_t>(index);
    }
    sortKeysByHighHalf();
    gatheredOrder_.clear();
    for (std::size_t begin = 0U; begin < count;) {
        std::size_t end = begin + 1U;
        while (end < count && (gatheredKeys_[end] >> 32U) == (gatheredKeys_[begin] >> 32U)) {
            ++end;
        }
        const std::size_t firstOfRun = gatheredOrder_.size();
        for (std::size_t key = begin; key < end; ++key) {
            const auto position = static_cast<std::uint32_t>(gatheredKeys_[key]);
            const GatheredCandidate& next = gathered_[position];
            bool repeated = false;
            for (std::size_t held = firstOfRun; held < gatheredOrder_.size(); ++held) {
                const GatheredCandidate& earlier = gathered_[gatheredOrder_[held]];
                if (earlier.setHash == next.setHash && earlier.size == next.size &&
                    std::ranges::equal(gatheredVertices(earlier), gatheredVertices(next))) {
                    repeated = true;
                    break;
                }
            }
            if (repeated) {
                ++stats.recovery.duplicates;
                continue;
            }
            gatheredOrder_.push_back(position);
        }
        begin = end;
    }
    std::size_t unique = gatheredOrder_.size();
    // Ones the batch already accepted -- in an earlier round of these, or
    // from a search under extension masks -- are in its table.  Only those
    // two can equal a deferred candidate (an undeferred one has a member of
    // lower degree), so without either the table is not read (about 3% of
    // flickr's recovery).
    if (filterGathered_) {
        std::size_t kept = 0U;
        for (const std::uint32_t index : gatheredOrder_) {
            const GatheredCandidate& candidate = gathered_[index];
            if (acceptedInBatch(gatheredVertices(candidate), candidate.setHash)) {
                ++stats.recovery.duplicates;
                continue;
            }
            gatheredOrder_[kept++] = index;
        }
        gatheredOrder_.resize(kept);
        unique = kept;
    }
    // A later round of this batch must look for this one's.
    filterGathered_ = true;
    // Checked in the order they were made: a destroyed clique's candidates
    // follow one another, and the lists its members' checks read are still
    // in cache from the one before.  Marked in a bitmap over the positions
    // and read back in order, instead of sorted.
    keptPositions_.assign((count + 63U) / 64U, 0U);
    for (const std::uint32_t position : gatheredOrder_) {
        keptPositions_[position >> 6U] |= std::uint64_t{1} << (position & 63U);
    }
    gatheredOrder_.clear();
    for (std::size_t word = 0U; word < keptPositions_.size(); ++word) {
        for (std::uint64_t bits = keptPositions_[word]; bits != 0U; bits &= bits - 1U) {
            gatheredOrder_.push_back(static_cast<std::uint32_t>(
                word * 64U + static_cast<std::size_t>(std::countr_zero(bits))));
        }
    }

    const auto settle = [&](std::uint32_t index, bool maximal) {
        const GatheredCandidate& candidate = gathered_[index];
        ++stats.recovery.checks;
        if (!maximal) {
            ++stats.recovery.rejected;
            return;
        }
        const std::span<const Vertex> vertices = gatheredVertices(candidate);
        recordAccepted(vertices, candidate.setHash);
        ++stats.recovery.accepted;
        output.emitCanonical(Clique(vertices.begin(), vertices.end()));
    };

    // One by one, each candidate's member list headers and then the lists'
    // first entries fetched ahead of it, as for the destroyed cliques.
    const auto checkOneByOne = [&](std::span<const std::uint32_t> order) {
        constexpr std::size_t headerAhead = 8U;
        constexpr std::size_t listAhead = 4U;
        const std::size_t total = order.size();
        for (std::size_t at = 0U; at < total; ++at) {
            if (at + headerAhead < total) {
                for (const Vertex vertex : gatheredVertices(gathered_[order[at + headerAhead]])) {
                    graph_.prefetchListHeader(vertex);
                }
            }
            if (at + listAhead < total) {
                for (const Vertex vertex : gatheredVertices(gathered_[order[at + listAhead]])) {
                    graph_.prefetchListStart(vertex);
                }
            }
            settle(order[at], maximalInGraph(gatheredVertices(gathered_[order[at]]), stats));
        }
    };
    const auto allOneByOne = [&] {
        checkOneByOne(gatheredOrder_);
        gathered_.clear();
        gatheredVertices_.clear();
        gatheredOrder_.clear();
    };
    const auto vertexCount = static_cast<std::size_t>(graph_.vertexCount());
    holdCount_.ensure(vertexCount);
    holdEpoch_.ensure(vertexCount);
    rowOfVertex_.ensure(vertexCount);
    rowEpoch_.ensure(vertexCount);

    // Each candidate's centre: its member that the most candidates hold,
    // the one of smaller degree (a narrower neighbourhood) on a tie.  A
    // group cannot outgrow the candidates holding its centre, so where no
    // vertex is held by `minRowGroup` of them nothing is grouped, and a
    // candidate none of whose members is takes its centre without the tie.
    advanceEpoch(holdEpoch_, holdEpochValue_);
    std::uint32_t mostHeld = 0U;
    for (const std::uint32_t index : gatheredOrder_) {
        for (const Vertex vertex : gatheredVertices(gathered_[index])) {
            const auto v = static_cast<std::size_t>(vertex);
            if (holdEpoch_[v] != holdEpochValue_) {
                holdEpoch_[v] = holdEpochValue_;
                holdCount_[v] = 0U;
            }
            mostHeld = std::max(mostHeld, ++holdCount_[v]);
        }
    }
    if (mostHeld < minRowGroup) {
        allOneByOne();
        return;
    }
    // By centre, in the order they were made within one: the centre above
    // the position, in order as the hashes above.  (A comparison sort of
    // (centre, position) pairs took most of the 10% of sx-stackoverflow's
    // recovery spent on the centres.)
    gatheredKeys_.clear();
    for (const std::uint32_t index : gatheredOrder_) {
        GatheredCandidate& candidate = gathered_[index];
        const std::span<const Vertex> vertices = gatheredVertices(candidate);
        Vertex centre = vertices[0];
        std::uint32_t best = holdCount_[static_cast<std::size_t>(centre)];
        for (const Vertex vertex : vertices.subspan(1)) {
            const std::uint32_t held = holdCount_[static_cast<std::size_t>(vertex)];
            if (held > best ||
                (held == best && held >= minRowGroup &&
                 graph_.degree(vertex) < graph_.degree(centre))) {
                centre = vertex;
                best = held;
            }
        }
        candidate.centre = centre;
        gatheredKeys_.push_back(
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(centre)) << 32U) | index);
    }
    sortKeysByHighHalf();
    centreOrder_.resize(gatheredKeys_.size());
    for (std::size_t key = 0U; key < gatheredKeys_.size(); ++key) {
        centreOrder_[key] = static_cast<std::uint32_t>(gatheredKeys_[key]);
    }

    // Which groups go on rows.  Reads either way: one by one, each
    // candidate scans its member of least degree and tests what it finds;
    // on rows, the centre's list once, each distinct member's row once, and
    // the words the ANDs touch.
    rowGroups_.clear();
    for (std::size_t begin = 0U; begin < unique;) {
        const Vertex centre = gathered_[centreOrder_[begin]].centre;
        std::size_t end = begin + 1U;
        while (end < unique && gathered_[centreOrder_[end]].centre == centre) {
            ++end;
        }
        const std::span<const std::uint32_t> group(
            centreOrder_.data() + begin, end - begin);
        const std::size_t groupBegin = begin;
        begin = end;
        if (group.size() < minRowGroup) {
            continue;
        }
        const auto centreDegree = static_cast<std::size_t>(graph_.degree(centre));
        const std::size_t words = (centreDegree + 63U) / 64U;
        long long leastSum = 0;
        auto rowReads = static_cast<long long>(centreDegree);
        long long andWords = 0;
        advanceEpoch(rowEpoch_, rowEpochValue_);
        for (const std::uint32_t index : group) {
            const std::span<const Vertex> vertices = gatheredVertices(gathered_[index]);
            long long least = std::numeric_limits<long long>::max();
            for (const Vertex vertex : vertices) {
                const auto degree = static_cast<std::size_t>(graph_.degree(vertex));
                least = std::min(least, static_cast<long long>(degree));
                const auto v = static_cast<std::size_t>(vertex);
                if (vertex != centre && rowEpoch_[v] != rowEpochValue_) {
                    rowEpoch_[v] = rowEpochValue_;
                    rowReads += rowReadCost(centreDegree, degree);
                }
            }
            leastSum += least;
            andWords += static_cast<long long>((vertices.size() - 1U) * words);
        }
        // Timed group by group both ways on nine cells: a check
        // one by one costs about two entries of its least list, a row about
        // its read, and an AND a thirty-second of an entry per word -- it
        // was counted a whole one, which kept flickr's groups off rows (3
        // times the best there, 1.03-1.12 times now).
        if (32 * rowReads + andWords >= 64 * leastSum) {
            continue;
        }
        rowGroups_.emplace_back(groupBegin, end);
        for (const std::uint32_t index : group) {
            gathered_[index].onRows = true;
        }
    }

    // The rest one by one, in the order they were made.  (By centre instead,
    // to keep the centre's list in cache from one check to the next, took
    // 1.003-1.019 of the time on sx-stackoverflow, youtube and wikipedia.)
    oneByOneOrder_.clear();
    for (const std::uint32_t index : gatheredOrder_) {
        if (!gathered_[index].onRows) {
            oneByOneOrder_.push_back(index);
        }
    }
    checkOneByOne(oneByOneOrder_);

    // On rows over the centre's neighbourhood, numbered by its list: for a
    // candidate the AND of its other members' rows -- the vertices of that
    // neighbourhood adjacent to all of them, members excluded since none is
    // its own neighbour.  Empty means nothing extends it.
    for (const auto& [groupBegin, groupEnd] : rowGroups_) {
        const std::span<const std::uint32_t> group(
            centreOrder_.data() + groupBegin, groupEnd - groupBegin);
        const Vertex centre = gathered_[group[0]].centre;
        const std::size_t words =
            (static_cast<std::size_t>(graph_.degree(centre)) + 63U) / 64U;
        ++stats.recovery.groups;
        advanceEpoch(rowEpoch_, rowEpochValue_);
        rowArena_.clear();
        rowScratch_.assign(words, 0U);
        liveWords_.resize(words);
        // Vertices that extended the group's latest candidates, as positions
        // in the centre's list, the newest first.  A group's candidates are
        // alike, and most are not maximal: one of these extends a candidate
        // when every other member's row holds it, read a bit per row and
        // given up at the first row without it, before any AND.
        constexpr std::size_t witnessSlots = 4U;
        std::array<std::uint32_t, witnessSlots> witnesses{};
        std::size_t witnessCount = 0U;
        for (const std::uint32_t index : group) {
            const std::span<const Vertex> members = gatheredVertices(gathered_[index]);
            bool witnessed = false;
            for (std::size_t slot = 0U; slot < witnessCount && !witnessed; ++slot) {
                const std::uint32_t at = witnesses[slot];
                witnessed = true;
                for (const Vertex vertex : members) {
                    if (vertex == centre) {
                        continue;
                    }
                    const std::uint64_t* row = groupRow(centre, vertex, words, stats);
                    if (((row[at >> 6U] >> (at & 63U)) & 1U) == 0U) {
                        witnessed = false;
                        break;
                    }
                }
            }
            if (witnessed) {
                ++stats.recovery.groupedChecks;
                settle(index, false);
                continue;
            }
            // Only the words still set are carried: the first row's nonzero
            // words (positions in `liveWords_`, values in `rowScratch_`),
            // thinned by every further row.  After two or three rows few are
            // left, where the whole width was ANDed for each row (47% of
            // flickr's recovery).  (In blocks of 32 words through
            // every row, stopping at the first block that keeps a word,
            // which cut the words a candidate that is not maximal ANDs by 90%
            // on wikipedia, took 1.003-1.058 of the time.)
            std::size_t live = 0U;
            bool first = true;
            bool empty = false;
            for (const Vertex vertex : members) {
                if (vertex == centre) {
                    continue;
                }
                const std::uint64_t* row = groupRow(centre, vertex, words, stats);
                if (first) {
                    first = false;
                    for (std::size_t word = 0U; word < words; ++word) {
                        liveWords_[live] = static_cast<std::uint32_t>(word);
                        rowScratch_[live] = row[word];
                        live += row[word] != 0U ? 1U : 0U;
                    }
                } else {
                    std::size_t remaining = 0U;
                    for (std::size_t at = 0U; at < live; ++at) {
                        const std::uint32_t word = liveWords_[at];
                        const std::uint64_t value = rowScratch_[at] & row[word];
                        liveWords_[remaining] = word;
                        rowScratch_[remaining] = value;
                        remaining += value != 0U ? 1U : 0U;
                    }
                    live = remaining;
                }
                if (live == 0U) {
                    empty = true;
                    break;
                }
            }
            if (!empty) {
                // The vertex of the lowest bit still set extends it.
                const auto at = static_cast<std::uint32_t>(
                    liveWords_[0] * 64U +
                    static_cast<std::uint32_t>(std::countr_zero(rowScratch_[0])));
                const std::size_t shift = std::min(witnessCount, witnessSlots - 1U);
                for (std::size_t slot = shift; slot > 0U; --slot) {
                    witnesses[slot] = witnesses[slot - 1U];
                }
                witnesses[0] = at;
                witnessCount = std::min(witnessCount + 1U, witnessSlots);
            }
            ++stats.recovery.groupedChecks;
            settle(index, empty);
        }
        if (centreNumbered_) {
            for (const AdjacencyEntry& entry : graph_.neighbors(centre)) {
                centreBits_[static_cast<std::size_t>(entry.neighbor) >> 6U] = 0U;
            }
            centreNumbered_ = false;
        }
    }
    gathered_.clear();
    gatheredVertices_.clear();
    gatheredOrder_.clear();
}

const std::uint64_t* SubsumedSearch::groupRow(
    Vertex centre,
    Vertex vertex,
    std::size_t words,
    SubsumedSearchStats& stats) {
    const auto v = static_cast<std::size_t>(vertex);
    if (rowEpoch_[v] == rowEpochValue_) {
        return rowArena_.data() + static_cast<std::size_t>(rowOfVertex_[v]) * words;
    }
    rowEpoch_[v] = rowEpochValue_;
    rowOfVertex_[v] = static_cast<std::uint32_t>(rowArena_.size() / words);
    rowArena_.resize(rowArena_.size() + words, 0U);
    std::uint64_t* row = rowArena_.data() + static_cast<std::size_t>(rowOfVertex_[v]) * words;
    const std::span<const AdjacencyEntry> list = graph_.neighbors(vertex);
    const std::span<const AdjacencyEntry> centreList = graph_.neighbors(centre);
    const auto mark = [row](std::size_t at) {
        row[at >> 6U] |= std::uint64_t{1} << (at & 63U);
    };
    // Each vertex of the centre's list searched for in the member's, when
    // that list is the longer by far; otherwise the member's list looked up
    // in the centre's neighbourhood.  (A merge of the two lists read both
    // and branched on every comparison: on youtube's reverse deletion the
    // rows took a fifth of the recovery, most of it in merges; reading them
    // all by lookup took 0.92-0.97 of the recovery on sx-stackoverflow,
    // youtube and wikipedia.)
    if (rowBySearch(centreList.size(), list.size())) {
        for (std::size_t at = 0U; at < centreList.size(); ++at) {
            ++stats.recovery.adjacencyTests;
            if (adjacentInList(list, centreList[at].neighbor)) {
                mark(at);
            }
        }
    } else {
        // On a bit per vertex, a thirty-second of a position per vertex,
        // which stays in cache where those positions missed it (youtube's
        // 3.2 million vertices: 0.955 of the recovery).
        if (!centreNumbered_) {
            centreNumbered_ = true;
            const std::size_t bitWords =
                (static_cast<std::size_t>(graph_.vertexCount()) + 63U) / 64U;
            if (centreBits_.size() != bitWords) {
                centreBits_.assign(bitWords, 0U);
                centreRank_.assign(bitWords, 0U);
            }
            std::uint32_t position = 0U;
            for (const AdjacencyEntry& entry : centreList) {
                const auto w = static_cast<std::size_t>(entry.neighbor);
                std::uint64_t& bits = centreBits_[w >> 6U];
                if (bits == 0U) {
                    centreRank_[w >> 6U] = position;
                }
                bits |= std::uint64_t{1} << (w & 63U);
                ++position;
            }
            stats.recovery.rowReads += static_cast<long long>(centreList.size());
        }
        // Without a branch on the bit: each entry's position is written to
        // `rowHits_` and kept by advancing over it only when the bit is set.
        rowHits_.resize(list.size() + 1U);
        std::uint32_t* hits = rowHits_.data();
        const std::uint64_t* centreBits = centreBits_.data();
        const std::uint32_t* centreRank = centreRank_.data();
        std::size_t hitCount = 0U;
        for (const AdjacencyEntry& entry : list) {
            const auto w = static_cast<std::uint32_t>(entry.neighbor);
            const std::uint64_t bits = centreBits[w >> 6U];
            const unsigned shift = w & 63U;
            hits[hitCount] = centreRank[w >> 6U] + static_cast<std::uint32_t>(
                std::popcount(bits & ((std::uint64_t{1} << shift) - 1U)));
            hitCount += static_cast<std::size_t>((bits >> shift) & 1U);
        }
        for (std::size_t hit = 0U; hit < hitCount; ++hit) {
            mark(hits[hit]);
        }
        stats.recovery.rowReads += static_cast<long long>(list.size());
    }
    return row;
}

void SubsumedSearch::acceptRecovered(
    std::span<const Vertex> clique,
    std::span<const std::uint64_t> maskWords,
    std::size_t memberCount,
    std::uint64_t setHash,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    recoveredScratch_.clear();
    for (std::size_t word = 0; word < maskWords.size(); ++word) {
        std::uint64_t bits = maskWords[word];
        while (bits != 0U) {
            const std::size_t index =
                word * 64U + static_cast<std::size_t>(std::countr_zero(bits));
            bits &= bits - 1U;
            if (index < clique.size()) {
                recoveredScratch_.push_back(clique[index]);
            }
        }
    }
    if (recoveredScratch_.size() != memberCount) {
        throw std::logic_error("a candidate mask disagrees with its member count");
    }
    acceptRecoveredClique(recoveredScratch_, setHash, stats, output);
}

// A candidate is a clique of the graph after the deletion: it holds no
// deleted edge, every one inside its destroyed clique having been split.  So
// it is a new maximal clique exactly when nothing extends it.
void SubsumedSearch::acceptRecoveredClique(
    std::span<const Vertex> candidate,
    std::uint64_t setHash,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    ++stats.recovery.candidates;
    if (candidate.size() < 2U) {
        // A vertex that keeps an edge is in a clique of two, and one that
        // keeps none is the maximal clique {v}, which the maintenance engine
        // lists from the vertex's degree (see MaintenanceEngine); the index
        // holds cliques of two or more vertices only.
        ++stats.recovery.belowTwo;
        return;
    }
    // A search under extension masks emits only maximal candidates; one can
    // still repeat a candidate another destroyed clique gave.
    if (knownMaximal_) {
        if (acceptedInBatch(candidate, setHash)) {
            ++stats.recovery.duplicates;
            return;
        }
        recordAccepted(candidate, setHash);
        filterGathered_ = true;
        ++stats.recovery.accepted;
        output.emitCanonical(Clique(candidate.begin(), candidate.end()));
        return;
    }
    // A candidate whose every member has degree `deferDegree` or more is
    // checked with the rest of such candidates at the end of the batch, or
    // when `gatherLimit` vertices are waiting.  Nothing in the index can
    // equal one -- it is a proper subset of a destroyed clique and holds no
    // inserted edge -- so the duplicates are the batch's own, found in its
    // table of accepted candidates then.
    if (readMemberLists(candidate) >= static_cast<std::size_t>(deferDegree)) {
        gathered_.push_back(GatheredCandidate{
            setHash,
            gatheredVertices_.size(),
            static_cast<std::uint32_t>(candidate.size()),
            0});
        gatheredVertices_.insert(
            gatheredVertices_.end(), candidate.begin(), candidate.end());
        if (gatheredVertices_.size() >= gatherLimit) {
            checkGatheredCandidates(stats, output);
        }
        return;
    }
    // The duplicates are found in the batch's table of accepted candidates,
    // the clique index being far larger and colder.
    if (acceptedInBatch(candidate, setHash)) {
        ++stats.recovery.duplicates;
        return;
    }
    ++stats.recovery.checks;
    // The candidates are the maximal cliques of the destroyed clique D in
    // the graph after the deletion, so a vertex of D left out of one is joined
    // to some member by a deleted edge and cannot extend it: the scan skips
    // all of D.
    const bool maximal = maximalOnMemberLists(
        candidate, !consumedClique_.empty() ? consumedClique_ : candidate, stats);
    if (!maximal) {
        ++stats.recovery.rejected;
        return;
    }
    // It enters the index at the end of the batch (`insertAcceptedInBatch`):
    // nothing reads the index before then.
    recordAccepted(candidate, setHash);
    ++stats.recovery.accepted;
    output.emitCanonical(Clique(candidate.begin(), candidate.end()));
}

// Any vertex extending the clique is a neighbour of every member.  So each
// neighbour of the member of least degree, the clique itself aside, is
// tested against the other members, those of smaller degree first (their
// lists are the shorter searches and the likelier to miss), and the first
// one adjacent to all of them settles it: the check of a candidate that is
// not maximal stops there instead of finding every vertex that extends it.
bool SubsumedSearch::maximalInGraph(
    std::span<const Vertex> candidate,
    SubsumedSearchStats& stats) {
    readMemberLists(candidate);
    return maximalOnMemberLists(candidate, candidate, stats);
}

std::size_t SubsumedSearch::readMemberLists(
    std::span<const Vertex> candidate) {
    memberLists_.resize(candidate.size());
    std::size_t least = std::numeric_limits<std::size_t>::max();
    for (std::size_t index = 0U; index < candidate.size(); ++index) {
        memberLists_[index] = graph_.neighborsUnchecked(candidate[index]);
        least = std::min(least, memberLists_[index].size());
    }
    return least;
}

namespace {

// Whether `list` holds `vertex`, from `cursor` on, leaving `cursor` at the
// first entry not below it.  The vertices one cursor is asked about only
// grow, so it gallops forward from where it stopped instead of searching
// the whole list again.
[[nodiscard]] bool advanceTo(
    std::span<const AdjacencyEntry> list,
    std::size_t& cursor,
    Vertex vertex) noexcept {
    const std::size_t size = list.size();
    std::size_t low = cursor;
    if (low >= size) {
        return false;
    }
    if (list[low].neighbor >= vertex) {
        return list[low].neighbor == vertex;
    }
    std::size_t step = 1U;
    std::size_t high = low + 1U;
    while (high < size && list[high].neighbor < vertex) {
        low = high;
        step <<= 1U;
        high = low + step;
    }
    high = std::min(high, size);
    // The first entry not below `vertex` is in [low + 1, high]; halved
    // without a branch on the comparison, which a random list mispredicts
    // half the time.
    const AdjacencyEntry* base = list.data() + low + 1U;
    std::size_t length = high - low - 1U;
    while (length > 1U) {
        const std::size_t half = length / 2U;
        base = base[half].neighbor < vertex ? base + half : base;
        length -= half;
    }
    std::size_t position = static_cast<std::size_t>(base - list.data());
    if (length == 1U && base->neighbor < vertex) {
        ++position;
    }
    cursor = position;
    return position < size && list[position].neighbor == vertex;
}

} // namespace

// Any vertex extending the clique is a neighbour of every member.  So each
// neighbour of the member of least degree, `skipped` aside, is tested
// against the other members, those of smaller degree
// first (their lists are the shorter searches and the likelier to miss),
// and the first one adjacent to all of them settles it: the check of a
// candidate that is not maximal stops there instead of finding every
// vertex that extends it.  The lists are those `readMemberLists` read,
// each searched by a galloping cursor.  (Trying the latest checks'
// extenders first, by binary searches of the members' lists, settled 22-30%
// of the checks one by one on youtube and wikipedia but took 1.000-1.029 of
// the time: they settle the ones whose scans stop early anyway.)
bool SubsumedSearch::maximalOnMemberLists(
    std::span<const Vertex> candidate,
    std::span<const Vertex> skipped,
    SubsumedSearchStats& stats) {
    const std::size_t count = candidate.size();
    memberOrder_.clear();
    for (std::size_t index = 0U; index < count; ++index) {
        memberOrder_.emplace_back(
            memberLists_[index].size(), static_cast<std::uint32_t>(index));
    }
    // A candidate holds a handful of members, too few for std::sort.
    if (count <= 32U) {
        for (std::size_t index = 1U; index < count; ++index) {
            const auto entry = memberOrder_[index];
            std::size_t place = index;
            while (place > 0U && entry < memberOrder_[place - 1U]) {
                memberOrder_[place] = memberOrder_[place - 1U];
                --place;
            }
            memberOrder_[place] = entry;
        }
    } else {
        std::ranges::sort(memberOrder_);
    }
    memberCursors_.assign(count, 0U);
    long long reads = 0;
    long long tests = 0;
    bool maximal = true;
    auto member = skipped.begin();
    for (const AdjacencyEntry& entry : memberLists_[memberOrder_[0].second]) {
        ++reads;
        while (member != skipped.end() && *member < entry.neighbor) {
            ++member;
        }
        if (member != skipped.end() && *member == entry.neighbor) {
            continue;
        }
        bool extends = true;
        for (std::size_t rank = 1U; rank < count; ++rank) {
            ++tests;
            if (!advanceTo(memberLists_[memberOrder_[rank].second],
                           memberCursors_[rank], entry.neighbor)) {
                extends = false;
                break;
            }
        }
        if (extends) {
            maximal = false;
            break;
        }
    }
    stats.recovery.neighborReads += reads;
    stats.recovery.adjacencyTests += tests;
    return maximal;
}

// KernelEnum's recovery.  Every candidate of a destroyed clique D holds its
// universal part U_D, and no vertex of D outside a candidate extends it (the
// candidates are the maximal cliques of D in the graph after the deletion).
// So a vertex extending a
// candidate is outside D and adjacent to all of U_D, and it extends exactly
// the candidates that keep none of the kernel vertices it misses: a
// candidate is maximal exactly when it keeps, for each such vertex, one of
// those.  The shortest list of U_D is scanned and its other lists tested by
// galloping cursors, as in a check; each vertex found is then tested against
// the lists of K_D.
template <typename OnExtender>
void SubsumedSearch::scanUniversalExtenders(
    const Clique& clique,
    std::span<const std::uint32_t> universal,
    std::span<const std::uint32_t> kernel,
    SubsumedSearchStats& stats,
    OnExtender&& onExtender) {
    universalOrder_.clear();
    for (const std::uint32_t index : universal) {
        universalOrder_.emplace_back(
            graph_.neighborsUnchecked(clique[index]).size(), index);
    }
    const std::size_t count = universalOrder_.size();
    if (count <= 32U) {
        for (std::size_t index = 1U; index < count; ++index) {
            const auto entry = universalOrder_[index];
            std::size_t place = index;
            while (place > 0U && entry < universalOrder_[place - 1U]) {
                universalOrder_[place] = universalOrder_[place - 1U];
                --place;
            }
            universalOrder_[place] = entry;
        }
    } else {
        std::ranges::sort(universalOrder_);
    }
    universalLists_.resize(count);
    for (std::size_t rank = 0U; rank < count; ++rank) {
        universalLists_[rank] = graph_.neighborsUnchecked(clique[universalOrder_[rank].second]);
    }
    universalCursors_.assign(count, 0U);
    const std::size_t kernelSize = kernel.size();
    kernelLists_.resize(kernelSize);
    for (std::size_t local = 0U; local < kernelSize; ++local) {
        kernelLists_[local] = graph_.neighborsUnchecked(clique[kernel[local]]);
    }
    kernelCursors_.assign(kernelSize, 0U);

    long long reads = 0;
    long long tests = 0;
    long long extenders = 0;
    auto member = clique.begin();
    for (const AdjacencyEntry& entry : universalLists_[0]) {
        ++reads;
        const Vertex vertex = entry.neighbor;
        while (member != clique.end() && *member < vertex) {
            ++member;
        }
        if (member != clique.end() && *member == vertex) {
            continue;
        }
        bool extendsUniversal = true;
        for (std::size_t rank = 1U; rank < count; ++rank) {
            ++tests;
            if (!advanceTo(universalLists_[rank], universalCursors_[rank], vertex)) {
                extendsUniversal = false;
                break;
            }
        }
        if (!extendsUniversal) {
            continue;
        }
        ++extenders;
        std::uint64_t missed = 0U;
        for (std::size_t local = 0U; local < kernelSize; ++local) {
            ++tests;
            if (!advanceTo(kernelLists_[local], kernelCursors_[local], vertex)) {
                missed |= std::uint64_t{1} << local;
            }
        }
        if (onExtender(missed)) {
            break;
        }
    }
    stats.recovery.neighborReads += reads;
    stats.recovery.adjacencyTests += tests;
    stats.recovery.universalExtenders += extenders;
}

bool SubsumedSearch::collectExtensionMasks(
    const Clique& clique,
    std::span<const std::uint32_t> universal,
    std::span<const std::uint32_t> kernel,
    SubsumedSearchStats& stats) {
    extensionMasks_.clear();
    bool extendsWhole = false;
    scanUniversalExtenders(clique, universal, kernel, stats, [&](std::uint64_t missed) {
        if (missed == 0U) {
            extendsWhole = true;
            return true;
        }
        extensionMasks_.push_back(missed);
        return false;
    });
    if (extendsWhole) {
        return false;
    }
    // The minimal masks, the smallest first: a candidate keeping a vertex of
    // a mask keeps one of every mask holding it.
    std::ranges::sort(extensionMasks_, [](std::uint64_t a, std::uint64_t b) {
        const int left = std::popcount(a);
        const int right = std::popcount(b);
        return left != right ? left < right : a < b;
    });
    std::size_t kept = 0U;
    for (std::size_t index = 0U; index < extensionMasks_.size(); ++index) {
        const std::uint64_t mask = extensionMasks_[index];
        bool implied = false;
        for (std::size_t smaller = 0U; smaller < kept; ++smaller) {
            if ((extensionMasks_[smaller] & ~mask) == 0U) {
                implied = true;
                break;
            }
        }
        if (!implied) {
            extensionMasks_[kept++] = mask;
        }
    }
    extensionMasks_.resize(kept);
    stats.recovery.extensionMasks += static_cast<long long>(kept);
    return true;
}

bool SubsumedSearch::acceptedInBatch(
    std::span<const Vertex> candidate,
    std::uint64_t setHash) const noexcept {
    const std::size_t mask = acceptedSlots_.size() - 1U;
    for (std::size_t slot = static_cast<std::size_t>(setHash) & mask;;
         slot = (slot + 1U) & mask) {
        const std::uint64_t value = acceptedSlots_[slot];
        if ((value >> 32U) != acceptedEpoch_) {
            return false;
        }
        const AcceptedRecord& record =
            acceptedRecords_[static_cast<std::size_t>(value & 0xffffffffU)];
        if (record.setHash == setHash && record.size == candidate.size() &&
            std::ranges::equal(
                std::span<const Vertex>(
                    acceptedVertices_.data() + record.offset, record.size),
                candidate)) {
            return true;
        }
    }
}

// The candidates accepted on the spot go into the index together, each
// slot prefetched 32 insertions ahead: one at a time, every insertion
// waited on a cache miss into a table far larger than the cache.
void SubsumedSearch::insertAcceptedInBatch(MaximalCliqueIndex& cliques) {
    constexpr std::size_t ahead = 32U;
    const std::size_t count = acceptedRecords_.size();
    for (std::size_t record = 0U; record < std::min(ahead, count); ++record) {
        cliques.prefetchForInsert(acceptedRecords_[record].setHash);
    }
    for (std::size_t record = 0U; record < count; ++record) {
        if (record + ahead < count) {
            cliques.prefetchForInsert(acceptedRecords_[record + ahead].setHash);
        }
        const AcceptedRecord& accepted = acceptedRecords_[record];
        if (!cliques.insertHashed(
                std::span<const Vertex>(
                    acceptedVertices_.data() + accepted.offset, accepted.size),
                accepted.setHash)) {
            throw std::logic_error("a recovered clique is already in the index");
        }
    }
}

void SubsumedSearch::recordAccepted(
    std::span<const Vertex> candidate,
    std::uint64_t setHash) {
    const auto place = [this](std::size_t record) {
        const std::size_t mask = acceptedSlots_.size() - 1U;
        std::size_t slot =
            static_cast<std::size_t>(acceptedRecords_[record].setHash) & mask;
        while ((acceptedSlots_[slot] >> 32U) == acceptedEpoch_) {
            slot = (slot + 1U) & mask;
        }
        acceptedSlots_[slot] =
            (static_cast<std::uint64_t>(acceptedEpoch_) << 32U) | record;
    };
    acceptedRecords_.push_back(AcceptedRecord{
        setHash,
        acceptedVertices_.size(),
        static_cast<std::uint32_t>(candidate.size())});
    acceptedVertices_.insert(acceptedVertices_.end(), candidate.begin(), candidate.end());
    if (2U * acceptedRecords_.size() > acceptedSlots_.size()) {
        // Doubled, which leaves every slot of an older batch: the epoch
        // tags tell them apart, and this batch's records go in again.
        acceptedSlots_.assign(2U * acceptedSlots_.size(), 0U);
        for (std::size_t record = 0U; record < acceptedRecords_.size(); ++record) {
            place(record);
        }
        return;
    }
    place(acceptedRecords_.size() - 1U);
}

// ------------------------------------------------- the bitmap extraction

void SubsumedSearch::prepareBitmapExtraction() {
    const std::size_t words =
        (static_cast<std::size_t>(graph_.vertexCount()) + 63U) / 64U;
    if (endpointBits_.size() != words) {
        endpointBits_.assign(words, 0U);
    }
    endpointList_.clear();
    std::size_t capacity = 64U;
    while (capacity < 4U * batchEdges_.size()) {
        capacity <<= 1U;
    }
    // The table is cleared slot by slot at release, so a batch pays for its
    // own edges only, not for the table's capacity.
    if (batchEdgeSet_.size() < capacity) {
        batchEdgeSet_.assign(capacity, 0U);
    } else {
        capacity = batchEdgeSet_.size();
    }
    batchEdgeSlots_.clear();
    batchEdgeSetShift_ =
        64U - static_cast<unsigned>(std::countr_zero(capacity));
    const std::size_t mask = capacity - 1U;
    for (std::size_t position = 0; position < batchEdges_.size(); ++position) {
        const EdgeId edgeId = batchEdges_[position];
        const Edge& edge = graph_.edge(edgeId);
        for (const Vertex endpoint : {edge.u, edge.v}) {
            const auto v = static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(endpoint));
            std::uint64_t& word = endpointBits_[v >> 6U];
            const std::uint64_t bit = std::uint64_t{1} << (v & 63U);
            if ((word & bit) == 0U) {
                word |= bit;
                endpointList_.push_back(endpoint);
            }
        }
        const Vertex low = std::min(edge.u, edge.v);
        const Vertex high = std::max(edge.u, edge.v);
        const std::uint64_t key =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(low)) << 32U) |
            static_cast<std::uint32_t>(high);
        std::size_t slot = static_cast<std::size_t>(
            (key * 0x9E3779B97F4A7C15ULL) >> batchEdgeSetShift_);
        while (batchEdgeSet_[slot] != 0U && batchEdgeSet_[slot] != key) {
            slot = (slot + 1U) & mask;
        }
        if (batchEdgeSet_[slot] == 0U) {
            batchEdgeSet_[slot] = key;
            batchEdgeSlots_.push_back(static_cast<std::uint32_t>(slot));
        }
    }
    bitmapExtractionPrepared_ = true;
}

void SubsumedSearch::releaseBitmapExtraction() {
    for (const Vertex endpoint : endpointList_) {
        const auto v = static_cast<std::uint64_t>(
            static_cast<std::uint32_t>(endpoint));
        endpointBits_[v >> 6U] &= ~(std::uint64_t{1} << (v & 63U));
    }
    endpointList_.clear();
    for (const std::uint32_t slot : batchEdgeSlots_) {
        batchEdgeSet_[slot] = 0U;
    }
    batchEdgeSlots_.clear();
    bitmapExtractionPrepared_ = false;
}

std::size_t SubsumedSearch::findBatchEdge(
    Vertex u, Vertex v) const noexcept {
    const Vertex low = std::min(u, v);
    const Vertex high = std::max(u, v);
    const std::uint64_t key =
        (static_cast<std::uint64_t>(static_cast<std::uint32_t>(low)) << 32U) |
        static_cast<std::uint32_t>(high);
    const std::size_t mask = batchEdgeSet_.size() - 1U;
    std::size_t slot = static_cast<std::size_t>(
        (key * 0x9E3779B97F4A7C15ULL) >> batchEdgeSetShift_);
    while (true) {
        const std::uint64_t entry = batchEdgeSet_[slot];
        if (entry == key) {
            return slot;
        }
        if (entry == 0U) {
            return noBatchEdgeSlot;
        }
        slot = (slot + 1U) & mask;
    }
}

bool SubsumedSearch::batchEdgeSetContains(
    Vertex u, Vertex v) const noexcept {
    return findBatchEdge(u, v) != noBatchEdgeSlot;
}

// ----------------------------------------------------------- KernelEnum

// One pass over C: its set hash (which the insertion reuses) and the batch
// endpoints S in C, by one bit test per vertex.
void SubsumedSearch::consumeNewCliqueKernelEnum(
    const Clique& clique,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    ++stats.newCliquesProcessed;
    const std::size_t n = clique.size();
    std::uint64_t setAccumulator = 0U;
    if (n - 1U < minimumCandidateSize(oldCliques)) {
        if (!recovering_) {
            for (Vertex vertex : clique) {
                setAccumulator ^= vertexSetWord(vertex);
            }
            const std::uint64_t setHash = finishSetHash(setAccumulator, n);
            newCliqueHashes_.push_back(setHash);
            oldCliques.prefetchForInsert(setHash);
        }
        ++stats.sizeGateSkippedCliques;
        return;
    }
    ensureBitmapExtraction();
    batchEndpointPositions_.clear();
    const std::uint64_t* bits = endpointBits_.data();
    for (std::size_t index = 0; index < n; ++index) {
        const Vertex vertex = clique[index];
        setAccumulator ^= vertexSetWord(vertex);
        const auto v = static_cast<std::uint64_t>(static_cast<std::uint32_t>(vertex));
        if (((bits[v >> 6U] >> (v & 63U)) & 1U) != 0U) {
            batchEndpointPositions_.push_back(static_cast<std::uint32_t>(index));
        }
    }
    // The clique's own hash serves its insertion after the batch, which a
    // recovery never makes.
    if (!recovering_) {
        const std::uint64_t setHash = finishSetHash(setAccumulator, n);
        newCliqueHashes_.push_back(setHash);
        oldCliques.prefetchForInsert(setHash);
    }
    findNewEdgesInside(
        clique, setAccumulator, oldCliques, stats, output);
}

// The new edges F_C inside C, the batch edges between the vertices of S:
// into `newEdges_`, and their endpoints' clique positions into
// `newEdgePositions_`.
void SubsumedSearch::findNewEdgesInside(
    const Clique& clique,
    std::uint64_t setAccumulator,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    const std::size_t n = clique.size();
    const std::size_t batchEndpoints = batchEndpointPositions_.size();
    if (batchEndpoints < 2U) {
        throw std::logic_error("a new maximal clique must contain a batch edge");
    }
    if (batchEndpoints == 2U) {
        // A new clique holds at least one new edge, and both of its
        // endpoints are batch endpoints in C.  With only two of them in C,
        // that edge is all of F_C -- no pair is looked up.
        ++stats.internalBatchEdges;
        const std::size_t a = batchEndpointPositions_[0];
        const std::size_t b = batchEndpointPositions_[1];
        if (n <= 64U) {
            probeOneEdgeKernel(
                clique, a, b, setAccumulator, oldCliques, stats, output);
            return;
        }
        newEdges_.assign(1U, InternalBatchEdge{clique[a], clique[b]});
        newEdgePositions_.assign(1U, {static_cast<std::uint32_t>(a),
                                         static_cast<std::uint32_t>(b)});
    } else {
        newEdges_.clear();
        newEdgePositions_.clear();
        // Every pair of S is looked up in the batch's edge table.  (Walking
        // the batch edges anchored at S instead took 1.38 times as long on
        // insertion batches, and gained nothing for the S of more than 64
        // vertices alone; in a recovery it cut youtube's pair tests by 39%
        // but took 1.030 of the time: a test of the edge table, a megabyte,
        // is cheap.)
        for (std::size_t i = 0; i + 1U < batchEndpoints; ++i) {
            const Vertex u = clique[batchEndpointPositions_[i]];
            for (std::size_t j = i + 1U; j < batchEndpoints; ++j) {
                ++stats.batchEdgeVisits;
                const Vertex v = clique[batchEndpointPositions_[j]];
                if (batchEdgeSetContains(u, v)) {
                    newEdges_.push_back(InternalBatchEdge{u, v});
                    newEdgePositions_.emplace_back(
                        batchEndpointPositions_[i], batchEndpointPositions_[j]);
                }
            }
        }
        stats.internalBatchEdges +=
            static_cast<long long>(newEdges_.size());
        if (newEdges_.empty()) {
            throw std::logic_error(
                "a new maximal clique must contain a batch edge");
        }
    }
    enumerateCandidates(clique, setAccumulator, oldCliques, stats, output);
}

// The candidates of C once F_C is known, by the size of C: a single new
// edge's two, `enumerateCandidatesNarrow` for a clique of at most 64
// vertices (a narrow one: its sets fit in one word), or
// `enumerateCandidatesWide`.
void SubsumedSearch::enumerateCandidates(
    const Clique& clique,
    std::uint64_t setAccumulator,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    const std::size_t n = clique.size();
    if (newEdges_.size() == 1U && n <= 64U) {
        probeOneEdgeKernel(
            clique, newEdgePositions_[0].first,
            newEdgePositions_[0].second, setAccumulator, oldCliques, stats,
            output);
        return;
    }
    if (n <= 64U) {
        std::array<std::uint64_t, 64U> newEdgeRows;  // rows of K_C only
        std::uint64_t kernel = 0U;
        for (const auto& [a, b] : newEdgePositions_) {
            kernel |= (std::uint64_t{1} << a) | (std::uint64_t{1} << b);
        }
        for (auto rest = kernel; rest != 0U; rest &= rest - 1U) {
            newEdgeRows[static_cast<std::size_t>(std::countr_zero(rest))] = 0U;
        }
        for (const auto& [a, b] : newEdgePositions_) {
            newEdgeRows[a] |= std::uint64_t{1} << b;
            newEdgeRows[b] |= std::uint64_t{1} << a;
        }
        enumerateCandidatesNarrow(clique, newEdgeRows, kernel, newEdges_.size(),
                            setAccumulator, oldCliques, stats, output);
        return;
    }
    enumerateCandidatesWide(clique, setAccumulator, oldCliques, stats, output);
}

void SubsumedSearch::enumerateCandidatesWide(
    const Clique& clique,
    std::uint64_t setAccumulator,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    const std::size_t n = clique.size();
    std::uint32_t arenaOffset = noArenaOffset;
    const auto emit = [&](std::span<const std::uint64_t> mask,
                          std::size_t count, std::uint64_t omittedHash) {
        ++stats.kernelCandidates;
        if (!candidateSizePossible(oldCliques, count)) {
            ++stats.sizeGateSkippedProbes;
            return;
        }
        submitProbe(clique, mask, count,
            finishSetHash(setAccumulator ^ omittedHash, count),
            arenaOffset, oldCliques, stats, output);
    };
    // New edges joining every pair of K_C have the closed form of
    // `enumerateCandidatesNarrow` here too: U_C with one kernel vertex each.
    // (These went through BK_pivot: on coauthor's sliding window the wide
    // cliques took 65% of the recovery, four in five of them with such
    // kernels.)
    newEdgeDegree_.assign(n, 0U);
    for (const auto& [a, b] : newEdgePositions_) {
        ++newEdgeDegree_[a];
        ++newEdgeDegree_[b];
    }
    kernelPositions_.clear();
    for (std::size_t position = 0U; position < n; ++position) {
        if (newEdgeDegree_[position] != 0U) {
            kernelPositions_.push_back(static_cast<std::uint32_t>(position));
        }
    }
    const std::size_t t = kernelPositions_.size();
    if (newEdgePositions_.size() == t * (t - 1U) / 2U) {
        ++stats.edgelessKernels;
        candidateMask_.assign((n + 63U) / 64U, ~std::uint64_t{0});
        if (n % 64U != 0U) {
            candidateMask_.back() = (std::uint64_t{1} << (n % 64U)) - 1U;
        }
        std::uint64_t kernelHash = 0U;
        for (const std::uint32_t position : kernelPositions_) {
            candidateMask_[position / 64U] &= ~(std::uint64_t{1} << (position % 64U));
            kernelHash ^= vertexSetWord(clique[position]);
        }
        // The kernel vertices a maximal candidate keeps: all of them, or, in
        // a recovery whose candidates' checks would read more than
        // `maskFactor` times the shortest list of U_C, those every vertex
        // extending U_C misses -- each candidate holds all of U_C, so one
        // scan settles them all, and what it emits needs no check.
        std::uint64_t alive = t >= 64U ? ~std::uint64_t{0} : (std::uint64_t{1} << t) - 1U;
        bool masked = false;
        if (recovering_ && t <= 64U && t < n) {
            universalPositions_.clear();
            std::size_t universalLeast = std::numeric_limits<std::size_t>::max();
            for (std::size_t position = 0U; position < n; ++position) {
                if (newEdgeDegree_[position] == 0U) {
                    universalPositions_.push_back(static_cast<std::uint32_t>(position));
                    universalLeast = std::min(
                        universalLeast, graph_.neighborsUnchecked(clique[position]).size());
                }
            }
            std::size_t spent = 0U;
            for (const std::uint32_t position : kernelPositions_) {
                spent += std::min(
                    universalLeast, graph_.neighborsUnchecked(clique[position]).size());
            }
            if (spent > maskFactor * universalLeast) {
                masked = true;
                ++stats.recovery.maskedSearches;
                scanUniversalExtenders(clique, universalPositions_, kernelPositions_, stats,
                    [&](std::uint64_t missed) {
                        alive &= missed;
                        return alive == 0U;
                    });
            }
        }
        knownMaximal_ = masked;
        for (std::size_t local = 0U; local < t; ++local) {
            if (t <= 64U && ((alive >> local) & 1U) == 0U) {
                continue;
            }
            const std::uint32_t position = kernelPositions_[local];
            candidateMask_[position / 64U] |= std::uint64_t{1} << (position % 64U);
            emit(std::span<const std::uint64_t>(candidateMask_), n - t + 1U,
                 kernelHash ^ vertexSetWord(clique[position]));
            candidateMask_[position / 64U] &= ~(std::uint64_t{1} << (position % 64U));
        }
        knownMaximal_ = false;
        return;
    }
    KernelEnumBkPivotWide<decltype(emit)> search(
        clique, std::span<const InternalBatchEdge>(newEdges_),
        minimumCandidateSize(oldCliques), stats.bkPivotWide, emit);
    search.run();
}

void SubsumedSearch::enumerateCandidatesNarrow(
    const Clique& clique,
    const std::array<std::uint64_t, 64U>& newEdgeRows,
    std::uint64_t kernel,
    std::size_t edgeCount,
    std::uint64_t setAccumulator,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    const auto t = static_cast<std::size_t>(std::popcount(kernel));
    const std::uint64_t full = clique.size() == 64U
        ? std::numeric_limits<std::uint64_t>::max()
        : (std::uint64_t{1} << clique.size()) - 1U;
    const std::uint64_t universal = full & ~kernel;
    const std::size_t universalSize = clique.size() - t;
    std::uint64_t kernelHash = 0U;
    for (auto rest = kernel; rest != 0U; rest &= rest - 1U) {
        kernelHash ^= vertexSetWord(clique[static_cast<std::size_t>(std::countr_zero(rest))]);
    }
    std::uint32_t arenaOffset = noArenaOffset;
    // `kept` is the part of K_C a candidate keeps; its hash is the clique's
    // with the dropped kernel vertices XORed back out.
    const auto emitCandidate = [&](std::uint64_t kept, std::uint64_t keptHash) {
        ++stats.kernelCandidates;
        const std::size_t count =
            universalSize + static_cast<std::size_t>(std::popcount(kept));
        if (!candidateSizePossible(oldCliques, count)) {
            ++stats.sizeGateSkippedProbes;
            return;
        }
        const std::uint64_t mask = universal | kept;
        submitProbe(clique, std::span<const std::uint64_t>(&mask, 1U), count,
            finishSetHash(setAccumulator ^ kernelHash ^ keptHash, count),
            arenaOffset, oldCliques, stats, output);
    };
    // The closed forms.  New edges joining every pair of K_C leave G[K_C]
    // without an edge: the candidates are U_C with one kernel vertex each.
    // New edges sharing one endpoint c leave G[K_C] the vertex c and the
    // clique K_C \ {c}: the candidates are U_C u {c} and U_C u (K_C \ {c}).
    if (edgeCount == t * (t - 1U) / 2U) {
        ++stats.edgelessKernels;
        // The kernel vertices a maximal candidate keeps: all, or in a
        // recovery whose candidates' checks would read more than `maskFactor`
        // times the shortest list of U_C, those every vertex extending U_C
        // misses, as for the wide cliques.
        std::uint64_t alive = kernel;
        bool masked = false;
        if (recovering_ && universalSize != 0U) {
            std::size_t universalLeast = std::numeric_limits<std::size_t>::max();
            universalPositions_.clear();
            for (auto rest = universal; rest != 0U; rest &= rest - 1U) {
                const auto i = static_cast<std::uint32_t>(std::countr_zero(rest));
                universalPositions_.push_back(i);
                universalLeast = std::min(universalLeast, graph_.neighborsUnchecked(clique[i]).size());
            }
            std::size_t spent = 0U;
            kernelPositions_.clear();
            for (auto rest = kernel; rest != 0U; rest &= rest - 1U) {
                const auto i = static_cast<std::uint32_t>(std::countr_zero(rest));
                kernelPositions_.push_back(i);
                spent += std::min(universalLeast, graph_.neighborsUnchecked(clique[i]).size());
            }
            if (spent > maskFactor * universalLeast) {
                masked = true;
                ++stats.recovery.maskedSearches;
                std::uint64_t kept = t >= 64U ? ~std::uint64_t{0} : (std::uint64_t{1} << t) - 1U;
                scanUniversalExtenders(clique, universalPositions_, kernelPositions_, stats,
                    [&](std::uint64_t missed) {
                        kept &= missed;
                        return kept == 0U;
                    });
                alive = 0U;
                for (auto rest = kept; rest != 0U; rest &= rest - 1U) {
                    alive |= std::uint64_t{1}
                        << kernelPositions_[static_cast<std::size_t>(std::countr_zero(rest))];
                }
            }
        }
        knownMaximal_ = masked;
        for (auto rest = alive; rest != 0U; rest &= rest - 1U) {
            const auto i = static_cast<std::size_t>(std::countr_zero(rest));
            emitCandidate(std::uint64_t{1} << i, vertexSetWord(clique[i]));
        }
        knownMaximal_ = false;
        return;
    }
    if (edgeCount == t - 1U) {
        for (auto rest = kernel; rest != 0U; rest &= rest - 1U) {
            const auto i = static_cast<std::size_t>(std::countr_zero(rest));
            if (static_cast<std::size_t>(std::popcount(newEdgeRows[i])) == t - 1U) {
                ++stats.sharedEndpointKernels;
                const std::uint64_t center = std::uint64_t{1} << i;
                const std::uint64_t centerHash = vertexSetWord(clique[i]);
                emitCandidate(center, centerHash);
                emitCandidate(kernel & ~center, kernelHash ^ centerHash);
                return;
            }
        }
    }
    kernelEnumBkPivotNarrow(clique, newEdgeRows, kernel, universal, kernelHash,
                        setAccumulator, oldCliques, stats, output);
}

// KernelEnum's BK_pivot on G[K_C] for a clique of at most 64 vertices.  The
// kernel is numbered 0..t-1 and every set is one word, all on the stack;
// a candidate goes out as one word over the positions in C.
void SubsumedSearch::kernelEnumBkPivotNarrow(
    const Clique& clique,
    const std::array<std::uint64_t, 64U>& newEdgeRows,
    std::uint64_t kernel,
    std::uint64_t universal,
    std::uint64_t kernelHash,
    std::uint64_t setAccumulator,
    MaximalCliqueIndex& oldCliques,
    SubsumedSearchStats& stats,
    CliqueBatchOutput& output) {
    const auto t = static_cast<std::size_t>(std::popcount(kernel));
    const std::size_t universalSize = clique.size() - t;
    std::uint32_t arenaOffset = noArenaOffset;
    // `kept` is the part of K_C a candidate keeps; its hash is the clique's
    // with the dropped kernel vertices XORed back out.
    const auto emitCandidate = [&](std::uint64_t kept, std::uint64_t keptHash) {
        ++stats.kernelCandidates;
        const std::size_t count =
            universalSize + static_cast<std::size_t>(std::popcount(kept));
        if (!candidateSizePossible(oldCliques, count)) {
            ++stats.sizeGateSkippedProbes;
            return;
        }
        const std::uint64_t mask = universal | kept;
        submitProbe(clique, std::span<const std::uint64_t>(&mask, 1U), count,
            finishSetHash(setAccumulator ^ kernelHash ^ keptHash, count),
            arenaOffset, oldCliques, stats, output);
    };
    // Local indexing of K_C (t vertices).  Only the first t entries of
    // `position`/`word` and the entries of `local` at kernel positions are
    // ever read.
    std::array<std::uint8_t, 64U> position;
    std::array<std::uint64_t, 64U> word;
    std::array<std::uint8_t, 64U> local;
    std::size_t next = 0U;
    for (auto rest = kernel; rest != 0U; rest &= rest - 1U) {
        const auto i = static_cast<std::size_t>(std::countr_zero(rest));
        position[next] = static_cast<std::uint8_t>(i);
        word[next] = vertexSetWord(clique[i]);
        local[i] = static_cast<std::uint8_t>(next);
        ++next;
    }
    // BK_pivot on G[K_C]: two kernel vertices are adjacent in G unless a new
    // edge joins them, so a vertex's row is K_C less itself and its new-edge
    // neighbours.  The same search as `KernelEnumBkPivotWide`, here in one
    // word.
    ++stats.kernelSearches;
    const std::uint64_t localFull = t == 64U
        ? std::numeric_limits<std::uint64_t>::max()
        : (std::uint64_t{1} << t) - 1U;
    std::array<std::uint64_t, 64U> kernelRows{};
    for (std::size_t i = 0U; i < t; ++i) {
        std::uint64_t row = localFull & ~(std::uint64_t{1} << i);
        for (auto rest = newEdgeRows[position[i]]; rest != 0U; rest &= rest - 1U) {
            row &= ~(std::uint64_t{1} << local[static_cast<std::size_t>(std::countr_zero(rest))]);
        }
        kernelRows[i] = row;
    }
    // KernelEnum's recovery with a non-empty U_C.  Checked one by one, a
    // candidate costs about the list of its member of least degree.  Under
    // extension masks -- a candidate is maximal exactly when it keeps a
    // vertex of each -- the destroyed clique costs one scan of the shortest
    // list of U_C, whole, and the search then cuts each branch that can no
    // longer keep a vertex of some mask and emits only maximal candidates.
    // So the search runs plain first, holding its candidates, and starts
    // again under the masks only once their lists add up to `maskFactor`
    // times that of U_C (the masks on every such clique cost 1.04-15 times
    // the checks on flickr, youtube, sx-stackoverflow and wikipedia, where
    // U_C has a vertex of high degree, and saved 2.8 times on coauthor's
    // exploding batches).
    bool masked = false;
    bool holding = false;
    bool overBudget = false;
    std::size_t budget = 0U;
    std::size_t spent = 0U;
    std::size_t universalLeast = std::numeric_limits<std::size_t>::max();
    std::array<std::size_t, 64U> kernelDegree;
    if (recovering_ && universalSize != 0U) {
        for (auto rest = universal; rest != 0U; rest &= rest - 1U) {
            universalLeast = std::min(universalLeast, graph_.neighborsUnchecked(
                clique[static_cast<std::size_t>(std::countr_zero(rest))]).size());
        }
        for (std::size_t i = 0U; i < t; ++i) {
            kernelDegree[i] = graph_.neighborsUnchecked(clique[position[i]]).size();
        }
        budget = maskFactor * universalLeast;
        holding = true;
        heldLeaves_.clear();
    }
    const std::size_t minimumSize = minimumCandidateSize(oldCliques);
    const auto visit = [&](auto&& self, std::uint64_t selected, std::uint64_t p,
                           std::uint64_t x) -> void {
        if (overBudget) {
            return;
        }
        ++stats.kernelSearchNodes;
        // |P| serves the pivot's early exits below.
        const int pCount = std::popcount(p);
        if (universalSize + static_cast<std::size_t>(std::popcount(selected) + pCount) <
            minimumSize) {
            return;
        }
        if (masked) {
            const std::uint64_t reachable = selected | p;
            for (const std::uint64_t missed : extensionMasks_) {
                if ((missed & reachable) == 0U) {
                    ++stats.recovery.maskCuts;
                    return;
                }
            }
        }
        if (p == 0U) {
            // Nothing left to add: `selected` is maximal exactly when X is
            // empty.  A non-empty X would only have made the pivot scan find
            // no extension.
            if (x != 0U) {
                return;
            }
            std::uint64_t kept = 0U;
            std::uint64_t keptHash = 0U;
            std::size_t least = universalLeast;
            for (auto rest = selected; rest != 0U; rest &= rest - 1U) {
                const auto i = static_cast<std::size_t>(std::countr_zero(rest));
                kept |= std::uint64_t{1} << position[i];
                keptHash ^= word[i];
                if (holding) {
                    least = std::min(least, kernelDegree[i]);
                }
            }
            if (holding) {
                spent += least;
                if (spent > budget) {
                    overBudget = true;
                    return;
                }
                heldLeaves_.emplace_back(kept, keptHash);
                return;
            }
            emitCandidate(kept, keptHash);
            return;
        }
        // Pivot: the first vertex of P u X, in bit order, with the most
        // neighbours in P.  Two scores end the scan early without changing
        // the pivot or anything emitted.  |P| is only reachable from X (no
        // row holds its own vertex) and means that vertex extends every set
        // below this call, so the full scan would find no extension: the
        // call ends.  |P|-1 is the most a vertex of P can score, so with X
        // empty no later vertex can beat the first one that reaches it.
        std::size_t pivot = 0U;
        int best = -1;
        for (auto rest = p | x; rest != 0U; rest &= rest - 1U) {
            const auto i = static_cast<std::size_t>(std::countr_zero(rest));
            const int count = std::popcount(p & kernelRows[i]);
            if (count == pCount) {
                return;
            }
            if (count > best) {
                best = count;
                pivot = i;
                if (x == 0U && count == pCount - 1) {
                    break;
                }
            }
        }
        auto extensions = p & ~kernelRows[pivot];
        while (extensions != 0U) {
            const std::uint64_t bit = extensions & (~extensions + 1U);
            extensions ^= bit;
            const auto i = static_cast<std::size_t>(std::countr_zero(bit));
            self(self, selected | bit, p & kernelRows[i], x & kernelRows[i]);
            p &= ~bit;
            x |= bit;
        }
    };
    visit(visit, 0U, localFull, 0U);
    if (!holding) {
        return;
    }
    holding = false;
    if (!overBudget) {
        for (const auto& [kept, keptHash] : heldLeaves_) {
            emitCandidate(kept, keptHash);
        }
        return;
    }
    universalPositions_.clear();
    for (auto rest = universal; rest != 0U; rest &= rest - 1U) {
        universalPositions_.push_back(static_cast<std::uint32_t>(std::countr_zero(rest)));
    }
    kernelPositions_.assign(
        position.begin(), position.begin() + static_cast<std::ptrdiff_t>(t));
    if (!collectExtensionMasks(clique, universalPositions_, kernelPositions_, stats)) {
        ++stats.recovery.extendedWhole;
        return;
    }
    ++stats.recovery.maskedSearches;
    masked = true;
    overBudget = false;
    knownMaximal_ = true;
    visit(visit, 0U, localFull, 0U);
    knownMaximal_ = false;
}

} // namespace ckm
