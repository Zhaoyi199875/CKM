#include "star_clique_search.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>

namespace immce {
namespace {

// The ownership threshold below the roots: a walk under it stops at the
// star's own edges and so moves only what earlier stars claim. That is all a
// walk can move there anyway: an own-star edge at a vertex w of a tree is
// (u, w) with w a later partner -- the earlier ones never leave the excluded
// set -- and u is already in the clique.
constexpr int earlierStarsOnly = std::numeric_limits<int>::min();

// The universe a subtree needs on a graph of its own: its candidates, its
// excluded set and the clique it extends, sorted.
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

} // namespace

const char* starBranchScopeName(StarBranchScope scope) noexcept {
    return scope == StarBranchScope::Root ? "root" : "unbound";
}

const char* starBranchingName(StarBranching branching) noexcept {
    switch (branching) {
    case StarBranching::Anchor:
        return "anchor";
    case StarBranching::Pivot:
        return "pivot";
    case StarBranching::Rule:
        break;
    }
    return "rule";
}

StarCliqueEnumerator::StarCliqueEnumerator(
    const IncrementalGraph& graph,
    StarSearchConfiguration configuration)
    : graph_(graph),
      configuration_(configuration),
      batchIndex_(graph),
      localGraph_(graph.vertexCount()),
      treeGraph_(graph.vertexCount()) {
    const auto vertexCount = static_cast<std::size_t>(
        graph.vertexCount() < 0 ? 0 : graph.vertexCount());
    starIndexByVertex_.assign(vertexCount, -1);
    endpointRankByVertex_.assign(vertexCount, -1);
    epochByVertex_.assign(vertexCount, 0U);
    unionEpochByVertex_.assign(vertexCount, 0U);
    // The branch-count rule bounds a vertex's cover by its degree in the
    // star's graph before it scans any row.
    localGraph_.setTrackDegrees(true);
}

std::vector<Clique> StarCliqueEnumerator::enumerateBatch(
    std::span<const EdgeId> batch,
    StarSearchStats& stats) {
    CliqueBatchOutput output(true);
    enumerateBatchInto(batch, stats, output);
    return std::move(output).takeCliques();
}

long long StarCliqueEnumerator::countBatch(
    std::span<const EdgeId> batch,
    StarSearchStats& stats) {
    CliqueBatchOutput output(false);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

long long StarCliqueEnumerator::visitBatch(
    std::span<const EdgeId> batch,
    StarSearchStats& stats,
    CliqueVisitor visitor) {
    if (!visitor) {
        throw std::invalid_argument("star-clique visitor must not be empty");
    }
    CliqueBatchOutput output(false, visitor);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

std::uint32_t StarCliqueEnumerator::beginEpoch() {
    if (epoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::ranges::fill(epochByVertex_, 0U);
        epoch_ = 0U;
    }
    return ++epoch_;
}

std::uint32_t StarCliqueEnumerator::beginUnionEpoch() {
    if (unionEpoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::ranges::fill(unionEpochByVertex_, 0U);
        unionEpoch_ = 0U;
    }
    return ++unionEpoch_;
}

const LocalVertexSet& StarCliqueEnumerator::partnerMaskFor(
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

// Does the subtree's own universe fit in fewer words than the graph it is on,
// and is it big enough to be searched long enough to repay the move.
bool StarCliqueEnumerator::reindexPays(
    const LocalSearchGraph& graph,
    const Clique& current,
    const LocalVertexSet& candidates,
    const LocalVertexSet& excluded) const noexcept {
    const std::size_t treeWords =
        (candidates.size() + excluded.size() + current.size() + 63U) / 64U;
    return candidates.size() >= configuration_.reindexMinCandidates &&
        graph.wordsPerRow() > treeWords;
}

// The child's graph is re-indexed from the bit rows of the one it leaves, so
// nothing reads the global adjacency again.
template <typename Recurse>
void StarCliqueEnumerator::descendInto(
    const LocalSearchGraph& graph,
    const Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    StarSearchStats& stats,
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
        // Every graph an unbound state runs on carries the mask of the star
        // being searched; one a bound subtree runs on may hold a stale one,
        // which nothing below reads, so it is not projected.
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
    // Taken only now: partnerMaskFor(graph) above still has to find the
    // parent's level among the live ones.
    ++deepLevel_;
    struct Release {
        std::size_t& level;
        ~Release() { --level; }
    } release{deepLevel_};
    recurse(compact, compactCandidates, compactExcluded);
}

void StarCliqueEnumerator::enumerateBatchInto(
    std::span<const EdgeId> batch,
    StarSearchStats& stats,
    CliqueBatchOutput& output) {
    if (batch.empty()) {
        return;
    }
    auto batchScope = batchIndex_.begin(batch);
    batch_ = batch;
    ++stats.batches;
    stats.batchEdges += static_cast<long long>(batch.size());

    // The partition runs once per batch, so timing it costs two clock reads
    // a batch.
    const auto assignStart = std::chrono::steady_clock::now();
    assignStars(batch, stats);
    stats.assignMilliseconds += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - assignStart).count();

    for (std::size_t starIndex = 0; starIndex < starHubs_.size(); ++starIndex) {
        const Vertex hub = starHubs_[starIndex];
        const int starRank = static_cast<int>(starIndex);
        if (starOffsets_[starIndex + 1U] - starOffsets_[starIndex] == 1U) {
            ++stats.singletonStars;
            ++stats.singletonStarEdges;
            searchSingleton(hub, starIndex, starRank, stats, output);
            continue;
        }
        buildDomain(hub, starIndex, stats);
        // One local graph per star instead of one per batch edge. The build
        // reads the adjacency of every domain vertex, and a centre's supports
        // overlap almost completely, so the edge-seeded search pays for that
        // same adjacency once per partner.
        localGraph_.rebuild(graph_, domain_, {}, {}, stats.localGraph);
        searchStar(hub, starIndex, starRank, stats, output);
    }

    for (Vertex hub : starHubs_) {
        starIndexByVertex_[static_cast<std::size_t>(hub)] = -1;
    }
    for (Vertex endpoint : endpoints_) {
        endpointRankByVertex_[static_cast<std::size_t>(endpoint)] = -1;
    }
    batch_ = {};
}

// Splits the batch into stars and ranks them.
//
// Centres. Walking the batch in order, every edge that is still uncovered puts
// both of its endpoints into the cover: a maximal matching, which is the
// Bar-Yehuda--Even local ratio scheme at unit weights and covers the batch
// edges with at most twice the fewest centres possible.
//
// Placement. An edge with one endpoint in the cover has no choice. An edge
// with both -- which a matching makes the common case -- goes to the endpoint
// carrying more batch edges, ties to the smaller vertex id.
//
// Ranks. A star ranks by the batch position of its first edge and an edge
// inside its star by batch order, so the order that decides ownership is the
// batch order the edge-seeded baseline uses.
void StarCliqueEnumerator::assignStars(
    std::span<const EdgeId> batch,
    StarSearchStats& stats) {
    starHubs_.clear();
    ownedCounts_.clear();
    endpoints_.clear();
    batchDegrees_.clear();
    ownerRankByPosition_.assign(batch.size(), -1);

    // Batch endpoints in first-appearance order, with their batch degrees.
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
                "vertex cover star partition left a batch edge uncovered");
        }
        Vertex hub = uCovered ? edge.u : edge.v;
        if (uCovered && vCovered) {
            ++stats.coverBothEndpoints;
            // Edges are stored with u < v, so a tie keeps u.
            hub = batchDegrees_[vRank] > batchDegrees_[uRank] ? edge.v : edge.u;
        }
        int& rank = starIndexByVertex_[static_cast<std::size_t>(hub)];
        if (rank < 0) {
            rank = static_cast<int>(starHubs_.size());
            starHubs_.push_back(hub);
            ownedCounts_.push_back(0U);
        }
        ownerRankByPosition_[position] = rank;
        ++ownedCounts_[static_cast<std::size_t>(rank)];
    }

    starOffsets_.assign(starHubs_.size() + 1U, 0U);
    for (std::size_t index = 0; index < ownedCounts_.size(); ++index) {
        starOffsets_[index + 1U] = starOffsets_[index] + ownedCounts_[index];
    }
    starFill_.assign(starOffsets_.begin(), starOffsets_.end() - 1);
    starEdgePositions_.assign(batch.size(), -1);
    orderByPosition_.assign(batch.size(), -1);
    for (std::size_t position = 0; position < batch.size(); ++position) {
        const auto star =
            static_cast<std::size_t>(ownerRankByPosition_[position]);
        const std::size_t slot = starFill_[star]++;
        starEdgePositions_[slot] = static_cast<int>(position);
        orderByPosition_[position] =
            static_cast<int>(slot - starOffsets_[star]);
    }

    buildIncidenceIndex(batch);

    stats.stars += static_cast<long long>(starHubs_.size());
    for (std::size_t index = 0; index < starHubs_.size(); ++index) {
        const std::size_t owned = starOffsets_[index + 1U] - starOffsets_[index];
        stats.starEdgeOccurrences += static_cast<long long>(owned);
        stats.maximumStarEdges = std::max(
            stats.maximumStarEdges, static_cast<long long>(owned));
    }
}

// Batch incidence, ordered by owning star rank and then by the edge's order
// within that star. The ownership walk wants exactly the entries ranking below
// a threshold, so this ordering lets it stop at the first entry that does not
// qualify -- usually the first one it looks at.
//
// The order comes for free: appending star by star, and within a star edge by
// edge, hands every vertex its own slice already sorted, so no vertex needs a
// sort of its own. That matters because there are about two entries per batch
// edge, and a sort call each would cost more than the walk it accelerates.
void StarCliqueEnumerator::buildIncidenceIndex(std::span<const EdgeId> batch) {
    // Each endpoint holds exactly as many entries as its batch degree, so the
    // offsets are a prefix sum of degrees already counted.
    incidenceOffsets_.assign(endpoints_.size() + 1U, 0U);
    for (std::size_t index = 0; index < batchDegrees_.size(); ++index) {
        incidenceOffsets_[index + 1U] =
            incidenceOffsets_[index] + batchDegrees_[index];
    }
    incidenceEntries_.assign(incidenceOffsets_.back(), StarIncidence{});
    incidenceFill_.assign(
        incidenceOffsets_.begin(), incidenceOffsets_.end() - 1);
    for (std::size_t star = 0; star < starHubs_.size(); ++star) {
        const std::size_t begin = starOffsets_[star];
        const std::size_t end = starOffsets_[star + 1U];
        for (std::size_t slot = begin; slot < end; ++slot) {
            const auto position =
                static_cast<std::size_t>(starEdgePositions_[slot]);
            const Edge& edge = graph_.edge(batch[position]);
            const StarIncidence base{
                0, static_cast<int>(star), static_cast<int>(slot - begin)};
            StarIncidence toV = base;
            toV.neighbor = edge.v;
            StarIncidence toU = base;
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

bool StarCliqueEnumerator::claimedEarlier(
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

Vertex StarCliqueEnumerator::partnerAt(
    std::size_t position,
    Vertex hub) const {
    const Edge& edge = graph_.edge(batch_[position]);
    return edge.u == hub ? edge.v : edge.u;
}

// The centre, its partners and every common neighbour of the centre and a
// partner. Intersection distributes over union,
//     N(u) n U_v N(v) = U_v (N(u) n N(v)),
// so the partners' neighbourhoods are marked and the centre's adjacency is
// read once per star, not once per partner.
void StarCliqueEnumerator::buildDomain(
    Vertex hub,
    std::size_t starIndex,
    StarSearchStats& stats) {
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
    mark(hub);

    const std::size_t begin = starOffsets_[starIndex];
    const std::size_t end = starOffsets_[starIndex + 1U];
    if (!configuration_.rootFilter) {
        // The ablation: the centre's whole closed neighbourhood, which holds
        // every partner, so no partner's adjacency is read here.
        for (std::size_t slot = begin; slot < end; ++slot) {
            const Vertex partner = partnerAt(
                static_cast<std::size_t>(starEdgePositions_[slot]), hub);
            partners_.push_back(partner);
            mark(partner);
        }
        stats.domainScanReads += static_cast<long long>(graph_.degree(hub));
        for (const AdjacencyEntry& entry : graph_.neighbors(hub)) {
            mark(entry.neighbor);
        }
        stats.domainVertexOccurrences +=
            static_cast<long long>(domain_.size());
        stats.maximumDomainVertices =
            std::max(stats.maximumDomainVertices, domain_.size());
        return;
    }
    const std::uint32_t unionEpoch = beginUnionEpoch();
    for (std::size_t slot = begin; slot < end; ++slot) {
        const auto position =
            static_cast<std::size_t>(starEdgePositions_[slot]);
        const Vertex partner = partnerAt(position, hub);
        partners_.push_back(partner);
        stats.domainScanReads +=
            static_cast<long long>(graph_.degree(partner));
        for (const AdjacencyEntry& entry : graph_.neighbors(partner)) {
            unionEpochByVertex_[
                static_cast<std::size_t>(entry.neighbor)] = unionEpoch;
            ++stats.domainMarkWrites;
        }
    }
    for (Vertex partner : partners_) {
        mark(partner);
    }
    stats.domainScanReads += static_cast<long long>(graph_.degree(hub));
    for (const AdjacencyEntry& entry : graph_.neighbors(hub)) {
        if (unionEpochByVertex_[
                static_cast<std::size_t>(entry.neighbor)] == unionEpoch) {
            mark(entry.neighbor);
        }
    }

    // rebuild() sorts and deduplicates the domain itself, and the partner
    // mask is set bit by bit, so neither list needs ordering here.
    stats.domainVertexOccurrences += static_cast<long long>(domain_.size());
    stats.maximumDomainVertices =
        std::max(stats.maximumDomainVertices, domain_.size());
}

void StarCliqueEnumerator::excludeClaimedNeighbors(
    const LocalSearchGraph& localGraph,
    Vertex vertex,
    int thresholdRank,
    int thresholdOrder,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    StarSearchStats& stats) const {
    const int vertexRank =
        endpointRankByVertex_[static_cast<std::size_t>(vertex)];
    if (vertexRank < 0) {
        // Not a batch endpoint, so it reaches no batch edge at all.
        return;
    }
    ++stats.ownershipWalks;
    const auto rank = static_cast<std::size_t>(vertexRank);
    const std::size_t end = incidenceOffsets_[rank + 1U];
    for (std::size_t index = incidenceOffsets_[rank]; index < end; ++index) {
        const StarIncidence& entry = incidenceEntries_[index];
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

// The one intersection a lone edge needs is its whole universe, and the
// ownership split falls out of the same pass, so the star degenerates to the
// edge-seeded search exactly instead of paying for machinery a lone edge
// cannot use.
void StarCliqueEnumerator::searchSingleton(
    Vertex hub,
    std::size_t starIndex,
    int starRank,
    StarSearchStats& stats,
    CliqueBatchOutput& output) {
    const auto position =
        static_cast<std::size_t>(starEdgePositions_[starOffsets_[starIndex]]);
    const Vertex partner = partnerAt(position, hub);

    candidateVertices_.clear();
    excludedVertices_.clear();
    graph_.forEachCommonNeighbor(
        hub,
        partner,
        [&](Vertex common, EdgeId fromHub, EdgeId fromPartner) {
            if (claimedEarlier(fromHub, starRank, earlierStarsOnly) ||
                claimedEarlier(fromPartner, starRank, earlierStarsOnly)) {
                excludedVertices_.push_back(common);
            } else {
                candidateVertices_.push_back(common);
            }
        });
    stats.domainScanReads += static_cast<long long>(
        graph_.degree(hub) + graph_.degree(partner));

    Clique current{hub, partner};
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
    branchBound(
        localGraph_, current, candidates, excluded, false, starRank, stats,
        output);
}

// The centre's row, split by whether an earlier star reaches the vertex from
// the centre by a batch edge.
void StarCliqueEnumerator::prepareRoot(
    Vertex hub,
    int starRank,
    StarSearchStats& stats) {
    allVertices_.vertices.clear();
    allVertices_.bits.assign(localGraph_.wordsPerRow(), ~std::uint64_t{0});
    localGraph_.intersectInto(allVertices_, hub, rootCandidates_);
    localGraph_.setBits({}, rootExcluded_);
    excludeClaimedNeighbors(
        localGraph_,
        hub,
        starRank,
        earlierStarsOnly,
        rootCandidates_,
        rootExcluded_,
        stats);
}

// One tree per partner v_j, in the star's order, each rooted at the edge
// (u, v_j). Its sets are the centre's cut by v_j's row, with the star's
// earlier partners moved to X -- they own every clique holding them -- and
// whatever an earlier star reaches from v_j. The centre's share was prepared
// once for the star; each root adds one intersection and the walk from v_j.
void StarCliqueEnumerator::searchPartnerRoots(
    Vertex hub,
    std::size_t starIndex,
    int starRank,
    StarSearchStats& stats,
    CliqueBatchOutput& output) {
    localGraph_.setBits({}, finishedPartners_);
    const std::span<const Vertex> localVertices = localGraph_.localVertices();
    const std::size_t begin = starOffsets_[starIndex];
    const std::size_t end = starOffsets_[starIndex + 1U];
    for (std::size_t slot = begin; slot < end; ++slot) {
        stats.perPartnerStarEdges += 1;
        const Vertex partner = partnerAt(
            static_cast<std::size_t>(starEdgePositions_[slot]), hub);
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
            starRank,
            earlierStarsOnly,
            partnerCandidates_,
            partnerExcluded_,
            stats);

        partnerClique_.assign({hub, partner});
        branchBound(
            localGraph_, partnerClique_, partnerCandidates_, partnerExcluded_,
            true, starRank, stats, output);

        const auto local =
            static_cast<std::size_t>(localGraph_.indexOf(partner));
        finishedPartners_.bits[local / 64U] |= std::uint64_t{1} << (local % 64U);
    }
}

// The centre's state is prepared once, both partitions are measured on it, and
// the chosen one is expanded from it directly: the new edges run the partner
// roots off this same state, the pivot expands it with the pivot the
// measurement already found. So nothing about the root is computed twice.
void StarCliqueEnumerator::searchStar(
    Vertex hub,
    std::size_t starIndex,
    int starRank,
    StarSearchStats& stats,
    CliqueBatchOutput& output) {
    prepareRoot(hub, starRank, stats);

    // Every partner is in P here -- the centre's walk moves only what an
    // earlier star claims, and each partner's edge to the centre is this
    // star's own -- so |W| is the partner count. The pivot is scanned into
    // scratch counters: a star that takes the new edges never has this node,
    // and one that takes the pivot has it as its root's own and is credited
    // with it below.
    struct ScanCounters {
        long long recursiveCalls = 0;
        long long pivotCalls = 0;
        long long pivotCandidatesEvaluated = 0;
        long long nodeWords = 0;
        long long pivotWords = 0;
    } scan;
    const std::size_t witnesses = partners_.size();

    // The ablation's fixed new edges need no pivot at all.
    if (configuration_.branching == StarBranching::Anchor) {
        stats.branchRootEdgeBranches += static_cast<long long>(witnesses);
        stats.branchRootChosenBranches += static_cast<long long>(witnesses);
        ++stats.branchRootEdgeStars;
        searchPartnerRoots(hub, starIndex, starRank, stats, output);
        return;
    }

    // The new edges win when the pivot's cover |P n N(q)| is at most
    // |P| - |W|. No vertex covers more of P than its degree in the domain less
    // one -- the centre, adjacent to all of them and not in P -- so when even
    // that bound keeps every vertex under the mark, the choice is made
    // without scanning a single row. Otherwise the scan runs, and it is the
    // root's own pivot scan, which the pivot side needs anyway. The
    // ablation's fixed pivot always scans and never compares.
    const bool fixedPivot = configuration_.branching == StarBranching::Pivot;
    const long long mark = static_cast<long long>(rootCandidates_.size()) -
        static_cast<long long>(witnesses);
    bool dense = fixedPivot;
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
    bool newEdges = true;
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
        newEdges = !fixedPivot && witnesses <= branches;
        ++stats.branchRootScans;
        stats.branchRootPivotBranches += static_cast<long long>(branches);
    } else {
        ++stats.branchRootDegreeDecisions;
    }
    stats.branchRootEdgeBranches += static_cast<long long>(witnesses);
    stats.branchRootChosenBranches +=
        static_cast<long long>(newEdges ? witnesses : branches);

    if (newEdges) {
        ++stats.branchRootEdgeStars;
        searchPartnerRoots(hub, starIndex, starRank, stats, output);
        return;
    }
    ++stats.branchRootPivotStars;
    stats.mergedStarEdges += static_cast<long long>(partners_.size());
    stats.pivotCalls += scan.pivotCalls;
    stats.pivotCandidatesEvaluated += scan.pivotCandidatesEvaluated;
    stats.pivotWords += scan.pivotWords;
    // Only the unbound states read the partner mask, so a star that takes the
    // new edges never builds it.
    localGraph_.setBits(partners_, partnerMask_);
    struct MaskScope {
        bool& active;
        explicit MaskScope(bool& flag) : active(flag) { active = true; }
        ~MaskScope() { active = false; }
        MaskScope(const MaskScope&) = delete;
        MaskScope& operator=(const MaskScope&) = delete;
    };
    const MaskScope maskScope(maskActive_);
    Clique current{hub};
    branchUnbound(
        localGraph_, current, rootCandidates_, rootExcluded_, true,
        BranchPartition::Pivot, pivot, starRank, stats, output);
}

// The two partitions of an unbound state. Both cover exactly the cliques the
// star still owes from here -- maximal, holding the current clique and at
// least one partner left in P -- and each covers every one of them once:
//  * new edges: a clique holds a first partner in the order taken, and every
//    earlier partner sits in X by the time that partner's child opens;
//  * pivot: a maximal clique holds a vertex of P \ N(q), or q would extend
//    it, and the usual move to X after each child keeps the children apart.
// So either may be taken at any unbound state, and the choice only changes
// how the work is cut up.
void StarCliqueEnumerator::branchUnbound(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    bool atRoot,
    BranchPartition partition,
    Vertex knownPivot,
    int starRank,
    StarSearchStats& stats,
    CliqueBatchOutput& output) {
    ++stats.recursiveCalls;
    ++stats.partnerlessNodes;
    stats.nodeWords += static_cast<long long>(graph.wordsPerRow());
    const LocalVertexSet& mask = partnerMaskFor(graph);
    std::size_t witnesses = 0U;
    for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
        witnesses += static_cast<std::size_t>(
            std::popcount(candidates.bits[word] & mask.bits[word]));
    }
    if (witnesses == 0U) {
        // The new-edge partition is empty: nothing below can hold a partner.
        ++stats.partnerPrunes;
        return;
    }

    if (unboundDepth_ == unboundFrames_.size()) {
        unboundFrames_.emplace_back();
    }
    UnboundFrame& frame = unboundFrames_[unboundDepth_];
    ++unboundDepth_;
    struct FrameRelease {
        std::size_t& depth;
        ~FrameRelease() { --depth; }
    } frameRelease{unboundDepth_};
    std::vector<Vertex>& extensions = frame.extensions;

    // The pivot is scanned only when a partition needs it: the pivot one
    // expands its branch set, and the rule compares that set's size against
    // the partners left. The rule makes the comparison on |P \ N(q)| counted
    // from the bits, so the branch list itself is built only if the pivot
    // wins.
    //
    // One partner left needs no scan. Binding it opens a single child, and
    // the rule takes that child whenever the pivot would open at least one;
    // a pivot opens none only when a vertex of X covers all of P, a state
    // that owes nothing, where the bound child costs one node and finds
    // nothing.
    Vertex pivot = knownPivot;
    bool newEdges = false;
    if (partition == BranchPartition::Rule) {
        if (witnesses == 1U) {
            newEdges = true;
            ++stats.branchUnitBindings;
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
            newEdges = witnesses <= branches;
        }
        if (newEdges) {
            ++stats.branchUnboundBindings;
        }
    }
    // Below the centre the rule is asked only by the method at the unbound
    // scope; the fixed pivot keeps the pivot at every unbound state.
    const BranchPartition below =
        configuration_.branching == StarBranching::Rule &&
                configuration_.branchScope == StarBranchScope::Unbound
        ? BranchPartition::Rule
        : BranchPartition::Pivot;

    LocalVertexSet& childCandidates = frame.childCandidates;
    LocalVertexSet& childExcluded = frame.childExcluded;
    const auto openChild = [&](Vertex vertex) {
        graph.intersectInto(candidates, vertex, childCandidates);
        graph.intersectInto(excluded, vertex, childExcluded);
        excludeClaimedNeighbors(
            graph,
            vertex,
            starRank,
            earlierStarsOnly,
            childCandidates,
            childExcluded,
            stats);
        current.push_back(vertex);
    };
    // Both partitions' branch lists are read before the first child, and
    // every child is cut from the bits, so the move need not touch the lists.
    const auto closeChild = [&](Vertex vertex) {
        current.pop_back();
        graph.moveToExcludedBits(candidates, excluded, vertex);
    };

    // A child the root's first-level moves below do not take goes through
    // descendInto, which moves it onto a narrower graph when the gate opens.
    const auto descendBound = [&] {
        descendInto(
            graph, current, childCandidates, childExcluded, stats,
            [&](const LocalSearchGraph& childGraph,
                LocalVertexSet& descendantCandidates,
                LocalVertexSet& descendantExcluded) {
                branchBound(
                    childGraph, current, descendantCandidates,
                    descendantExcluded, false, starRank, stats, output);
            });
    };
    const auto descendUnbound = [&] {
        descendInto(
            graph, current, childCandidates, childExcluded, stats,
            [&](const LocalSearchGraph& childGraph,
                LocalVertexSet& descendantCandidates,
                LocalVertexSet& descendantExcluded) {
                branchUnbound(
                    childGraph, current, descendantCandidates,
                    descendantExcluded, false, below, -1, starRank, stats,
                    output);
            });
    };

    if (newEdges) {
        // Only below the centre: the centre's own choice is searchStar's.
        // Any fixed order partitions what the state owes -- each child takes
        // the cliques whose first partner it is, the earlier ones waiting in
        // X -- so the partners still in P are read straight off the bits the
        // witness count just scanned, in vertex order.
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
            descendBound();
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
                branchBound(
                    graph, current, childCandidates, childExcluded, true,
                    starRank, stats, output);
            } else {
                descendBound();
            }
        } else if (atRoot && &graph == &localGraph_ &&
                   reindexPays(graph, current, childCandidates, childExcluded)) {
            // The centre's first-level move: the subtree only touches the
            // centre, this vertex and their common neighbours, which can be
            // many words narrower than the star.
            collectTreeVertices(
                current, childCandidates, childExcluded, treeVertices_);
            treeGraph_.rebuildFrom(localGraph_, treeVertices_, stats.localGraph);
            treeGraph_.projectSetFrom(partnerMask_, treePartnerMask_);
            ++stats.mergedReindexedSubtrees;
            LocalVertexSet compactCandidates =
                treeGraph_.makeSet(childCandidates.vertices);
            LocalVertexSet compactExcluded =
                treeGraph_.makeSet(childExcluded.vertices);
            branchUnbound(
                treeGraph_, current, compactCandidates, compactExcluded, false,
                below, -1, starRank, stats, output);
        } else {
            descendUnbound();
        }
        closeChild(vertex);
    }
}

void StarCliqueEnumerator::branchBound(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    bool mayReindex,
    int starRank,
    StarSearchStats& stats,
    CliqueBatchOutput& output) {
    const auto refineChild = [&, this](
        const LocalSearchGraph& childGraph,
        Vertex vertex,
        LocalVertexSet& childCandidates,
        LocalVertexSet& childExcluded,
        std::span<const EdgeId>) {
        excludeClaimedNeighbors(
            childGraph,
            vertex,
            starRank,
            earlierStarsOnly,
            childCandidates,
            childExcluded,
            stats);
    };
    const auto onLeaf = [&](
        std::span<const Vertex> clique,
        const EmptySearchState&) {
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
        bronKerboschLocalInPlace<false>(
            on,
            current,
            rootCandidates,
            rootExcluded,
            EmptySearchState{},
            NoStateUpdate{},
            refineChild,
            NoBranchPruning{},
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
        // Nothing in a bound subtree reads the mask, but a deeper move out of
        // treeGraph_ projects whatever mask it holds, so it must be this one.
        treeGraph_.projectSetFrom(partnerMask_, treePartnerMask_);
    }
    ++stats.reindexedTrees;
    // Both callers rebuild these sets from bits before they read them again,
    // so their vertex lists can be handed over rather than copied.
    compactCandidates_.vertices.swap(candidates.vertices);
    compactExcluded_.vertices.swap(excluded.vertices);
    treeGraph_.synchronize(compactCandidates_);
    treeGraph_.synchronize(compactExcluded_);
    search(treeGraph_, compactCandidates_, compactExcluded_);
}

} // namespace immce
