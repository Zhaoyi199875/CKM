#include "cover_seed.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ckm {
namespace {

constexpr int earlierSeedsOnly = std::numeric_limits<int>::min();

void collectTreeVertices(
    const Clique& current,
    const LocalVertexSet& candidates,
    const LocalVertexSet& excluded,
    std::vector<Vertex>& vertices) {
    vertices.clear();
    std::ranges::merge(
        candidates.vertices, excluded.vertices, std::back_inserter(vertices));
    for (Vertex member : current) {
        insertSortedUnique(vertices, member);
    }
}

}

CoverSeed::CoverSeed(
    const IncrementalGraph& graph,
    CoverSeedConfiguration configuration)
    : graph_(graph),
      configuration_(configuration),
      batchIndex_(graph),
      localGraph_(graph.vertexCount()),
      treeGraph_(graph.vertexCount()) {
    const auto vertexCount = static_cast<std::size_t>(
        graph.vertexCount() < 0 ? 0 : graph.vertexCount());
    seedIndexByVertex_.assign(vertexCount, -1);
    endpointRankByVertex_.assign(vertexCount, -1);
    epochByVertex_.assign(vertexCount, 0U);
    unionEpochByVertex_.assign(vertexCount, 0U);
    localGraph_.setTrackDegrees(true);
}

std::vector<Clique> CoverSeed::enumerateBatch(
    std::span<const EdgeId> batch,
    CoverSeedStats& stats) {
    CliqueBatchOutput output;
    enumerateBatchInto(batch, stats, output);
    return std::move(output).takeCliques();
}

std::uint32_t CoverSeed::beginEpoch() {
    if (epoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::ranges::fill(epochByVertex_, 0U);
        epoch_ = 0U;
    }
    return ++epoch_;
}

std::uint32_t CoverSeed::beginUnionEpoch() {
    if (unionEpoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::ranges::fill(unionEpochByVertex_, 0U);
        unionEpoch_ = 0U;
    }
    return ++unionEpoch_;
}

const LocalVertexSet& CoverSeed::partnerMaskFor(
    const LocalSearchGraph& graph) const noexcept {
    if (&graph == &treeGraph_) {
        return treePartnerMask_;
    }
    for (std::size_t level = 0; level < deepLevel_; ++level) {
        if (&graph == deepGraphs_[level].get()) {
            return deepMasks_[level];
        }
    }
    return partnerMask_;
}

bool CoverSeed::reindexPays(
    const LocalSearchGraph& graph,
    const Clique& current,
    const LocalVertexSet& candidates,
    const LocalVertexSet& excluded) const noexcept {
    const std::size_t treeWords =
        (candidates.size() + excluded.size() + current.size() + 63U) / 64U;
    return configuration_.reindex && graph.wordsPerRow() > treeWords;
}

template <typename Recurse>
void CoverSeed::descendInto(
    const LocalSearchGraph& graph,
    const Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    CoverSeedStats& stats,
    const Recurse& recurse) {
    if (!reindexPays(graph, current, candidates, excluded)) {
        recurse(graph, candidates, excluded);
        return;
    }
    if (deepLevel_ == deepGraphs_.size()) {
        deepGraphs_.push_back(
            std::make_unique<LocalSearchGraph>(graph_.vertexCount()));
        deepMasks_.emplace_back();
        deepVertices_.emplace_back();
    }
    const std::size_t level = deepLevel_;
    LocalSearchGraph& compact = *deepGraphs_[level];
    std::vector<Vertex>& vertices = deepVertices_[level];
    collectTreeVertices(current, candidates, excluded, vertices);
    compact.rebuildFrom(graph, vertices, stats.localGraph);
    if (maskActive_) {
        const LocalVertexSet& parentMask = partnerMaskFor(graph);
        if (parentMask.bits.size() == graph.wordsPerRow()) {
            compact.projectSetFrom(parentMask, deepMasks_[level]);
        } else {
            deepMasks_[level].vertices.clear();
            deepMasks_[level].bits.assign(compact.wordsPerRow(), 0U);
        }
    }
    ++stats.deepReindexedSubtrees;
    LocalVertexSet compactCandidates = compact.makeSet(candidates.vertices);
    LocalVertexSet compactExcluded = compact.makeSet(excluded.vertices);
    ++deepLevel_;
    struct Release {
        std::size_t& level;
        ~Release() { --level; }
    } release{deepLevel_};
    recurse(compact, compactCandidates, compactExcluded);
}

void CoverSeed::enumerateBatchInto(
    std::span<const EdgeId> batch,
    CoverSeedStats& stats,
    CliqueBatchOutput& output) {
    if (batch.empty()) {
        return;
    }
    auto batchScope = batchIndex_.begin(batch);
    batch_ = batch;
    ++stats.batches;
    stats.batchEdges += static_cast<long long>(batch.size());

    const auto assignStart = std::chrono::steady_clock::now();
    assignEdges(batch, stats);
    stats.assignMilliseconds += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - assignStart).count();

    for (std::size_t seedIndex = 0; seedIndex < seeds_.size(); ++seedIndex) {
        const Vertex seed = seeds_[seedIndex];
        const int seedRank = static_cast<int>(seedIndex);
        if (seedOffsets_[seedIndex + 1U] - seedOffsets_[seedIndex] == 1U) {
            ++stats.singleEdgeSeeds;
            ++stats.singleEdgeSeedEdges;
            searchSingleEdgeSeed(seed, seedIndex, seedRank, stats, output);
            continue;
        }
        buildDomain(seed, seedIndex, stats);
        localGraph_.rebuild(graph_, domain_, {}, {}, stats.localGraph);
        searchSeed(seed, seedIndex, seedRank, stats, output);
    }

    for (Vertex seed : seeds_) {
        seedIndexByVertex_[static_cast<std::size_t>(seed)] = -1;
    }
    for (Vertex endpoint : endpoints_) {
        endpointRankByVertex_[static_cast<std::size_t>(endpoint)] = -1;
    }
    batch_ = {};
}

void CoverSeed::assignEdges(
    std::span<const EdgeId> batch,
    CoverSeedStats& stats) {
    seeds_.clear();
    ownedCounts_.clear();
    endpoints_.clear();
    batchDegrees_.clear();
    ownerRankByPosition_.assign(batch.size(), -1);

    for (std::size_t position = 0; position < batch.size(); ++position) {
        const Edge& edge = graph_.edge(batch[position]);
        for (const Vertex endpoint : {edge.u, edge.v}) {
            int& rank = endpointRankByVertex_[static_cast<std::size_t>(endpoint)];
            if (rank < 0) {
                rank = static_cast<int>(endpoints_.size());
                endpoints_.push_back(endpoint);
                batchDegrees_.push_back(0U);
            }
            ++batchDegrees_[static_cast<std::size_t>(rank)];
        }
    }

    const auto rankOf = [this](Vertex vertex) {
        return static_cast<std::size_t>(
            endpointRankByVertex_[static_cast<std::size_t>(vertex)]);
    };
    inCover_.assign(endpoints_.size(), 0U);
    for (std::size_t position = 0; position < batch.size(); ++position) {
        const Edge& edge = graph_.edge(batch[position]);
        const std::size_t uRank = rankOf(edge.u);
        const std::size_t vRank = rankOf(edge.v);
        if (inCover_[uRank] == 0U && inCover_[vRank] == 0U) {
            inCover_[uRank] = 1U;
            inCover_[vRank] = 1U;
        }
    }
    for (const unsigned char member : inCover_) {
        stats.coverVertices += member != 0U ? 1 : 0;
    }

    for (std::size_t position = 0; position < batch.size(); ++position) {
        const Edge& edge = graph_.edge(batch[position]);
        const std::size_t uRank = rankOf(edge.u);
        const std::size_t vRank = rankOf(edge.v);
        const bool uCovered = inCover_[uRank] != 0U;
        const bool vCovered = inCover_[vRank] != 0U;
        if (!uCovered && !vCovered) {
            throw std::logic_error(
                "the seed cover left a batch edge uncovered");
        }
        Vertex seed = uCovered ? edge.u : edge.v;
        if (uCovered && vCovered) {
            ++stats.coverBothEndpoints;
            seed = batchDegrees_[vRank] > batchDegrees_[uRank] ? edge.v : edge.u;
        }
        int& rank = seedIndexByVertex_[static_cast<std::size_t>(seed)];
        if (rank < 0) {
            rank = static_cast<int>(seeds_.size());
            seeds_.push_back(seed);
            ownedCounts_.push_back(0U);
        }
        ownerRankByPosition_[position] = rank;
        ++ownedCounts_[static_cast<std::size_t>(rank)];
    }

    seedOffsets_.assign(seeds_.size() + 1U, 0U);
    for (std::size_t index = 0; index < ownedCounts_.size(); ++index) {
        seedOffsets_[index + 1U] = seedOffsets_[index] + ownedCounts_[index];
    }
    seedFill_.assign(seedOffsets_.begin(), seedOffsets_.end() - 1);
    seedEdgePositions_.assign(batch.size(), -1);
    orderByPosition_.assign(batch.size(), -1);
    for (std::size_t position = 0; position < batch.size(); ++position) {
        const auto owner =
            static_cast<std::size_t>(ownerRankByPosition_[position]);
        const std::size_t slot = seedFill_[owner]++;
        seedEdgePositions_[slot] = static_cast<int>(position);
        orderByPosition_[position] =
            static_cast<int>(slot - seedOffsets_[owner]);
    }

    buildIncidenceIndex(batch);

    stats.seeds += static_cast<long long>(seeds_.size());
    for (std::size_t index = 0; index < seeds_.size(); ++index) {
        const std::size_t owned = seedOffsets_[index + 1U] - seedOffsets_[index];
        stats.assignedEdges += static_cast<long long>(owned);
        stats.maximumSeedEdges = std::max(
            stats.maximumSeedEdges, static_cast<long long>(owned));
    }
}

void CoverSeed::buildIncidenceIndex(std::span<const EdgeId> batch) {
    incidenceOffsets_.assign(endpoints_.size() + 1U, 0U);
    for (std::size_t index = 0; index < batchDegrees_.size(); ++index) {
        incidenceOffsets_[index + 1U] =
            incidenceOffsets_[index] + batchDegrees_[index];
    }
    incidenceEntries_.assign(incidenceOffsets_.back(), SeedIncidence{});
    incidenceFill_.assign(
        incidenceOffsets_.begin(), incidenceOffsets_.end() - 1);
    for (std::size_t rank = 0; rank < seeds_.size(); ++rank) {
        const std::size_t begin = seedOffsets_[rank];
        const std::size_t end = seedOffsets_[rank + 1U];
        for (std::size_t slot = begin; slot < end; ++slot) {
            const auto position =
                static_cast<std::size_t>(seedEdgePositions_[slot]);
            const Edge& edge = graph_.edge(batch[position]);
            const SeedIncidence base{
                0, static_cast<int>(rank), static_cast<int>(slot - begin)};
            SeedIncidence toV = base;
            toV.neighbor = edge.v;
            SeedIncidence toU = base;
            toU.neighbor = edge.u;
            incidenceEntries_[incidenceFill_[
                static_cast<std::size_t>(endpointRankByVertex_[
                    static_cast<std::size_t>(edge.u)])]++] = toV;
            incidenceEntries_[incidenceFill_[
                static_cast<std::size_t>(endpointRankByVertex_[
                    static_cast<std::size_t>(edge.v)])]++] = toU;
        }
    }
}

bool CoverSeed::claimedEarlier(
    EdgeId edgeId,
    int thresholdRank,
    int thresholdOrder) const noexcept {
    const int position = batchIndex_.edgeRank(edgeId);
    if (position < 0) {
        return false;
    }
    const int owner = ownerRankByPosition_[static_cast<std::size_t>(position)];
    if (owner != thresholdRank) {
        return owner < thresholdRank;
    }
    return orderByPosition_[static_cast<std::size_t>(position)] <
        thresholdOrder;
}

Vertex CoverSeed::partnerAt(
    std::size_t position,
    Vertex seed) const {
    const Edge& edge = graph_.edge(batch_[position]);
    return edge.u == seed ? edge.v : edge.u;
}

void CoverSeed::buildDomain(
    Vertex seed,
    std::size_t seedIndex,
    CoverSeedStats& stats) {
    const std::uint32_t epoch = beginEpoch();
    domain_.clear();
    partners_.clear();
    const auto mark = [&](Vertex vertex) {
        std::uint32_t& marked = epochByVertex_[static_cast<std::size_t>(vertex)];
        if (marked != epoch) {
            marked = epoch;
            domain_.push_back(vertex);
        }
    };
    mark(seed);

    const std::size_t begin = seedOffsets_[seedIndex];
    const std::size_t end = seedOffsets_[seedIndex + 1U];
    const std::uint32_t unionEpoch = beginUnionEpoch();
    for (std::size_t slot = begin; slot < end; ++slot) {
        const auto position =
            static_cast<std::size_t>(seedEdgePositions_[slot]);
        const Vertex partner = partnerAt(position, seed);
        partners_.push_back(partner);
        for (const AdjacencyEntry& entry : graph_.neighbors(partner)) {
            unionEpochByVertex_[
                static_cast<std::size_t>(entry.neighbor)] = unionEpoch;
        }
    }
    for (Vertex partner : partners_) {
        mark(partner);
    }
    for (const AdjacencyEntry& entry : graph_.neighbors(seed)) {
        if (unionEpochByVertex_[
                static_cast<std::size_t>(entry.neighbor)] == unionEpoch) {
            mark(entry.neighbor);
        }
    }

    stats.domainVertexOccurrences += static_cast<long long>(domain_.size());
    stats.maximumDomainVertices =
        std::max(stats.maximumDomainVertices, domain_.size());
}

void CoverSeed::excludeClaimedNeighbors(
    const LocalSearchGraph& localGraph,
    Vertex vertex,
    int thresholdRank,
    int thresholdOrder,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    CoverSeedStats& stats) const {
    const int vertexRank =
        endpointRankByVertex_[static_cast<std::size_t>(vertex)];
    if (vertexRank < 0) {
        return;
    }
    const auto rank = static_cast<std::size_t>(vertexRank);
    const std::size_t end = incidenceOffsets_[rank + 1U];
    for (std::size_t index = incidenceOffsets_[rank]; index < end; ++index) {
        const SeedIncidence& entry = incidenceEntries_[index];
        if (entry.ownerRank > thresholdRank ||
            (entry.ownerRank == thresholdRank &&
             entry.order >= thresholdOrder)) {
            break;
        }
        ++stats.ownershipProbes;
        if (!localGraph.contains(candidates, entry.neighbor)) {
            continue;
        }
        localGraph.erase(candidates, entry.neighbor);
        localGraph.insert(excluded, entry.neighbor);
        ++stats.ownershipMoves;
    }
}

void CoverSeed::searchSingleEdgeSeed(
    Vertex seed,
    std::size_t seedIndex,
    int seedRank,
    CoverSeedStats& stats,
    CliqueBatchOutput& output) {
    const auto position =
        static_cast<std::size_t>(seedEdgePositions_[seedOffsets_[seedIndex]]);
    const Vertex partner = partnerAt(position, seed);

    candidateVertices_.clear();
    excludedVertices_.clear();
    graph_.forEachCommonNeighbor(
        seed,
        partner,
        [&](Vertex common, EdgeId fromSeed, EdgeId fromPartner) {
            if (claimedEarlier(fromSeed, seedRank, earlierSeedsOnly) ||
                claimedEarlier(fromPartner, seedRank, earlierSeedsOnly)) {
                excludedVertices_.push_back(common);
            } else {
                candidateVertices_.push_back(common);
            }
        });

    Clique current{seed, partner};
    localGraph_.rebuild(
        graph_,
        current,
        candidateVertices_,
        excludedVertices_,
        stats.localGraph);
    stats.domainVertexOccurrences += static_cast<long long>(
        localGraph_.localVertices().size());
    stats.maximumDomainVertices = std::max(
        stats.maximumDomainVertices, localGraph_.localVertices().size());

    LocalVertexSet candidates = localGraph_.makeSet(candidateVertices_);
    LocalVertexSet excluded = localGraph_.makeSet(excludedVertices_);
    branchAnchored(
        localGraph_, current, candidates, excluded, false, seedRank, stats,
        output);
}

void CoverSeed::prepareRoot(
    Vertex seed,
    int seedRank,
    CoverSeedStats& stats) {
    allVertices_.vertices.clear();
    allVertices_.bits.assign(localGraph_.wordsPerRow(), ~std::uint64_t{0});
    localGraph_.intersectInto(allVertices_, seed, rootCandidates_);
    localGraph_.setBits({}, rootExcluded_);
    excludeClaimedNeighbors(
        localGraph_,
        seed,
        seedRank,
        earlierSeedsOnly,
        rootCandidates_,
        rootExcluded_,
        stats);
}

void CoverSeed::searchPartnerRoots(
    Vertex seed,
    std::size_t seedIndex,
    int seedRank,
    CoverSeedStats& stats,
    CliqueBatchOutput& output) {
    localGraph_.setBits({}, finishedPartners_);
    const std::span<const Vertex> localVertices = localGraph_.localVertices();
    const std::size_t begin = seedOffsets_[seedIndex];
    const std::size_t end = seedOffsets_[seedIndex + 1U];
    for (std::size_t slot = begin; slot < end; ++slot) {
        stats.anchorRootEdges += 1;
        const Vertex partner = partnerAt(
            static_cast<std::size_t>(seedEdgePositions_[slot]), seed);
        localGraph_.intersectInto(rootCandidates_, partner, partnerCandidates_);
        localGraph_.intersectInto(rootExcluded_, partner, partnerExcluded_);
        for (std::size_t word = 0; word < partnerCandidates_.bits.size();
             ++word) {
            std::uint64_t claimed =
                partnerCandidates_.bits[word] & finishedPartners_.bits[word];
            while (claimed != 0U) {
                const std::size_t local = word * 64U +
                    static_cast<std::size_t>(std::countr_zero(claimed));
                claimed &= claimed - 1U;
                ++stats.ownershipProbes;
                localGraph_.erase(partnerCandidates_, localVertices[local]);
                localGraph_.insert(partnerExcluded_, localVertices[local]);
                ++stats.ownershipMoves;
            }
        }
        excludeClaimedNeighbors(
            localGraph_,
            partner,
            seedRank,
            earlierSeedsOnly,
            partnerCandidates_,
            partnerExcluded_,
            stats);

        partnerClique_.assign({seed, partner});
        branchAnchored(
            localGraph_, partnerClique_, partnerCandidates_, partnerExcluded_,
            true, seedRank, stats, output);

        const auto local =
            static_cast<std::size_t>(localGraph_.indexOf(partner));
        finishedPartners_.bits[local / 64U] |= std::uint64_t{1} << (local % 64U);
    }
}

void CoverSeed::searchSeed(
    Vertex seed,
    std::size_t seedIndex,
    int seedRank,
    CoverSeedStats& stats,
    CliqueBatchOutput& output) {
    prepareRoot(seed, seedRank, stats);

    struct ScanCounters {
        long long recursiveCalls = 0;
        long long pivotCalls = 0;
        long long pivotCandidatesEvaluated = 0;
        long long nodeWords = 0;
        long long pivotWords = 0;
    } scan;
    const std::size_t anchorCount = partners_.size();

    const long long mark = static_cast<long long>(rootCandidates_.size()) -
        static_cast<long long>(anchorCount);
    bool dense = false;
    for (const LocalVertexSet* set : {&rootExcluded_, &rootCandidates_}) {
        for (std::size_t word = 0; word < set->bits.size() && !dense; ++word) {
            std::uint64_t bits = set->bits[word];
            while (bits != 0U) {
                const std::size_t local = word * 64U +
                    static_cast<std::size_t>(std::countr_zero(bits));
                bits &= bits - 1U;
                if (static_cast<long long>(localGraph_.localDegree(local)) - 1 >
                    mark) {
                    dense = true;
                    break;
                }
            }
        }
    }
    bool anchorBranching = true;
    Vertex pivot = -1;
    std::size_t branches = 0U;
    if (dense) {
        pivot = chooseLocalBkPivot(
            localGraph_, rootCandidates_, rootExcluded_, scan);
        const std::span<const std::uint64_t> pivotRow = localGraph_.rowAt(
            static_cast<std::size_t>(localGraph_.indexOf(pivot)));
        for (std::size_t word = 0; word < rootCandidates_.bits.size(); ++word) {
            branches += static_cast<std::size_t>(
                std::popcount(rootCandidates_.bits[word] & ~pivotRow[word]));
        }
        anchorBranching = anchorCount <= branches;
        ++stats.rootPivotScans;
        stats.rootPivotBranches += static_cast<long long>(branches);
    } else {
        ++stats.rootDegreeDecisions;
    }
    stats.rootAnchorBranches += static_cast<long long>(anchorCount);
    stats.rootChosenBranches +=
        static_cast<long long>(anchorBranching ? anchorCount : branches);

    if (anchorBranching) {
        ++stats.anchorRootSeeds;
        searchPartnerRoots(seed, seedIndex, seedRank, stats, output);
        return;
    }
    ++stats.pivotRootSeeds;
    stats.pivotRootEdges += static_cast<long long>(partners_.size());
    stats.pivotCalls += scan.pivotCalls;
    stats.pivotCandidatesEvaluated += scan.pivotCandidatesEvaluated;
    stats.pivotWords += scan.pivotWords;
    localGraph_.setBits(partners_, partnerMask_);
    struct MaskScope {
        bool& active;
        explicit MaskScope(bool& flag) : active(flag) { active = true; }
        ~MaskScope() { active = false; }
        MaskScope(const MaskScope&) = delete;
        MaskScope& operator=(const MaskScope&) = delete;
    };
    const MaskScope maskScope(maskActive_);
    Clique current{seed};
    branchUnanchored(
        localGraph_, current, rootCandidates_, rootExcluded_, true,
        BranchPartition::Pivot, pivot, seedRank, stats, output);
}

void CoverSeed::branchUnanchored(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    bool atRoot,
    BranchPartition partition,
    Vertex knownPivot,
    int seedRank,
    CoverSeedStats& stats,
    CliqueBatchOutput& output) {
    ++stats.recursiveCalls;
    ++stats.unanchoredNodes;
    stats.nodeWords += static_cast<long long>(graph.wordsPerRow());
    const LocalVertexSet& mask = partnerMaskFor(graph);
    std::size_t anchorCount = 0U;
    for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
        anchorCount += static_cast<std::size_t>(
            std::popcount(candidates.bits[word] & mask.bits[word]));
    }
    if (anchorCount == 0U) {
        ++stats.liveTestPrunes;
        return;
    }

    if (unanchoredDepth_ == unanchoredFrames_.size()) {
        unanchoredFrames_.emplace_back();
    }
    UnanchoredFrame& frame = unanchoredFrames_[unanchoredDepth_];
    ++unanchoredDepth_;
    struct FrameRelease {
        std::size_t& depth;
        ~FrameRelease() { --depth; }
    } frameRelease{unanchoredDepth_};
    std::vector<Vertex>& extensions = frame.extensions;

    Vertex pivot = knownPivot;
    bool anchorBranching = false;
    if (partition == BranchPartition::AnchorOrPivot) {
        if (anchorCount == 1U) {
            anchorBranching = true;
            ++stats.singlePartnerAnchors;
        } else {
            if (pivot < 0) {
                pivot = chooseLocalBkPivot(graph, candidates, excluded, stats);
            }
            std::size_t branches = candidates.size();
            if (pivot >= 0) {
                const std::span<const std::uint64_t> pivotRow = graph.rowAt(
                    static_cast<std::size_t>(graph.indexOf(pivot)));
                branches = 0U;
                for (std::size_t word = 0; word < candidates.bits.size();
                     ++word) {
                    branches += static_cast<std::size_t>(std::popcount(
                        candidates.bits[word] & ~pivotRow[word]));
                }
            }
            anchorBranching = anchorCount <= branches;
        }
        if (anchorBranching) {
            ++stats.unanchoredAnchorChoices;
        }
    }
    constexpr BranchPartition below = BranchPartition::AnchorOrPivot;

    LocalVertexSet& childCandidates = frame.childCandidates;
    LocalVertexSet& childExcluded = frame.childExcluded;
    const auto openChild = [&](Vertex vertex) {
        graph.intersectInto(candidates, vertex, childCandidates);
        graph.intersectInto(excluded, vertex, childExcluded);
        excludeClaimedNeighbors(
            graph,
            vertex,
            seedRank,
            earlierSeedsOnly,
            childCandidates,
            childExcluded,
            stats);
        current.push_back(vertex);
    };
    const auto closeChild = [&](Vertex vertex) {
        current.pop_back();
        graph.moveToExcludedBits(candidates, excluded, vertex);
    };

    const auto descendAnchored = [&] {
        descendInto(
            graph, current, childCandidates, childExcluded, stats,
            [&](const LocalSearchGraph& childGraph,
                LocalVertexSet& descendantCandidates,
                LocalVertexSet& descendantExcluded) {
                branchAnchored(
                    childGraph, current, descendantCandidates,
                    descendantExcluded, false, seedRank, stats, output);
            });
    };
    const auto descendUnanchored = [&] {
        descendInto(
            graph, current, childCandidates, childExcluded, stats,
            [&](const LocalSearchGraph& childGraph,
                LocalVertexSet& descendantCandidates,
                LocalVertexSet& descendantExcluded) {
                branchUnanchored(
                    childGraph, current, descendantCandidates,
                    descendantExcluded, false, below, -1, seedRank, stats,
                    output);
            });
    };

    if (anchorBranching) {
        extensions.clear();
        const std::span<const Vertex> locals = graph.localVertices();
        for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
            std::uint64_t bits = candidates.bits[word] & mask.bits[word];
            while (bits != 0U) {
                extensions.push_back(locals[word * 64U +
                    static_cast<std::size_t>(std::countr_zero(bits))]);
                bits &= bits - 1U;
            }
        }
        for (Vertex partner : extensions) {
            openChild(partner);
            descendAnchored();
            closeChild(partner);
        }
        return;
    }

    if (pivot < 0) {
        pivot = chooseLocalBkPivot(graph, candidates, excluded, stats);
    }
    if (pivot < 0) {
        extensions.assign(candidates.vertices.begin(), candidates.vertices.end());
    } else {
        graph.differenceVerticesInto(candidates, pivot, extensions);
    }
    for (Vertex vertex : extensions) {
        if (!graph.contains(candidates, vertex)) {
            continue;
        }
        openChild(vertex);
        if (graph.contains(mask, vertex)) {
            if (atRoot) {
                branchAnchored(
                    graph, current, childCandidates, childExcluded, true,
                    seedRank, stats, output);
            } else {
                descendAnchored();
            }
        } else if (atRoot && &graph == &localGraph_ &&
                   reindexPays(graph, current, childCandidates, childExcluded)) {
            collectTreeVertices(
                current, childCandidates, childExcluded, treeVertices_);
            treeGraph_.rebuildFrom(localGraph_, treeVertices_, stats.localGraph);
            treeGraph_.projectSetFrom(partnerMask_, treePartnerMask_);
            ++stats.rootChildReindexedSubtrees;
            LocalVertexSet compactCandidates =
                treeGraph_.makeSet(childCandidates.vertices);
            LocalVertexSet compactExcluded =
                treeGraph_.makeSet(childExcluded.vertices);
            branchUnanchored(
                treeGraph_, current, compactCandidates, compactExcluded, false,
                below, -1, seedRank, stats, output);
        } else {
            descendUnanchored();
        }
        closeChild(vertex);
    }
}

void CoverSeed::branchAnchored(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    bool mayReindex,
    int seedRank,
    CoverSeedStats& stats,
    CliqueBatchOutput& output) {
    const auto refineChild = [&, this](
        const LocalSearchGraph& childGraph,
        Vertex vertex,
        LocalVertexSet& childCandidates,
        LocalVertexSet& childExcluded) {
        excludeClaimedNeighbors(
            childGraph,
            vertex,
            seedRank,
            earlierSeedsOnly,
            childCandidates,
            childExcluded,
            stats);
    };
    const auto onLeaf = [&](std::span<const Vertex> clique) {
        ++stats.leaves;
        output.emit(clique);
        ++stats.emittedCliques;
    };
    const auto descent = [this, &stats](
        const LocalSearchGraph& childGraph,
        const Clique& clique,
        LocalVertexSet& childCandidates,
        LocalVertexSet& childExcluded,
        const auto& recurse) {
        descendInto(
            childGraph, clique, childCandidates, childExcluded, stats, recurse);
    };
    const auto search = [&](
        const LocalSearchGraph& on,
        LocalVertexSet& rootCandidates,
        LocalVertexSet& rootExcluded) {
        bronKerboschLocalInPlace(
            on,
            current,
            rootCandidates,
            rootExcluded,
            refineChild,
            onLeaf,
            descent,
            stats);
    };
    if (!mayReindex || &graph != &localGraph_ ||
        !reindexPays(graph, current, candidates, excluded)) {
        search(graph, candidates, excluded);
        return;
    }
    collectTreeVertices(current, candidates, excluded, treeVertices_);
    treeGraph_.rebuildFrom(localGraph_, treeVertices_, stats.localGraph);
    if (maskActive_) {
        treeGraph_.projectSetFrom(partnerMask_, treePartnerMask_);
    }
    ++stats.reindexedTrees;
    compactCandidates_.vertices.swap(candidates.vertices);
    compactExcluded_.vertices.swap(excluded.vertices);
    treeGraph_.synchronize(compactCandidates_);
    treeGraph_.synchronize(compactExcluded_);
    search(treeGraph_, compactCandidates_, compactExcluded_);
}

}
