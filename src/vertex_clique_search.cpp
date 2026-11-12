#include "vertex_clique_search.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>
#include <utility>

namespace immce {
namespace {

// What a vertex-seeded search carries down: whether the clique holds a batch
// edge yet. Once it does, every maximal clique below qualifies and nothing is
// checked again.
struct VertexSearchState {
    bool holdsBatchEdge = false;
};

} // namespace

const char* vertexBranchingName(VertexBranching branching) noexcept {
    switch (branching) {
    case VertexBranching::Pivot:
        return "pivot";
    case VertexBranching::NewEdges:
        return "new-edges";
    case VertexBranching::RuleRoot:
        return "rule-root";
    case VertexBranching::RuleUnbound:
        return "rule-unbound";
    }
    return "unknown";
}

const char* vertexSeedsName(VertexSeeds seeds) noexcept {
    switch (seeds) {
    case VertexSeeds::Endpoints:
        return "endpoints";
    case VertexSeeds::Cover:
        return "cover";
    }
    return "unknown";
}

VertexCliqueEnumerator::VertexCliqueEnumerator(
    const IncrementalGraph& graph,
    VertexSearchConfiguration configuration)
    : graph_(graph),
      configuration_(configuration),
      batchIndex_(graph),
      localGraph_(graph.vertexCount()) {
    const auto vertexCount = static_cast<std::size_t>(
        graph.vertexCount() < 0 ? 0 : graph.vertexCount());
    endpointRankByVertex_.assign(vertexCount, -1);
    epochByVertex_.assign(vertexCount, 0U);
    // The rule's degree bound at the root reads the degrees the build counts.
    localGraph_.setTrackDegrees(
        configuration_.branching == VertexBranching::RuleRoot ||
        configuration_.branching == VertexBranching::RuleUnbound);
}

std::vector<Clique> VertexCliqueEnumerator::enumerateBatch(
    std::span<const EdgeId> batch,
    VertexSearchStats& stats) {
    CliqueBatchOutput output(true);
    enumerateBatchInto(batch, stats, output);
    return std::move(output).takeCliques();
}

long long VertexCliqueEnumerator::countBatch(
    std::span<const EdgeId> batch,
    VertexSearchStats& stats) {
    CliqueBatchOutput output(false);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

long long VertexCliqueEnumerator::visitBatch(
    std::span<const EdgeId> batch,
    VertexSearchStats& stats,
    CliqueVisitor visitor) {
    if (!visitor) {
        throw std::invalid_argument("vertex-clique visitor must not be empty");
    }
    CliqueBatchOutput output(false, visitor);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

std::uint32_t VertexCliqueEnumerator::beginEpoch() {
    if (epoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::ranges::fill(epochByVertex_, 0U);
        epoch_ = 0U;
    }
    return ++epoch_;
}

void VertexCliqueEnumerator::enumerateBatchInto(
    std::span<const EdgeId> batch,
    VertexSearchStats& stats,
    CliqueBatchOutput& output) {
    if (batch.empty()) {
        return;
    }
    auto batchScope = batchIndex_.begin(batch);
    ++stats.batches;
    stats.batchEdges += static_cast<long long>(batch.size());
    indexBatch(batch);
    for (std::size_t rank = 0; rank < endpoints_.size(); ++rank) {
        if (isSeed(static_cast<int>(rank))) {
            enumerateRoot(endpoints_[rank], static_cast<int>(rank), stats, output);
        }
    }
    for (Vertex endpoint : endpoints_) {
        endpointRankByVertex_[static_cast<std::size_t>(endpoint)] = -1;
    }
}

void VertexCliqueEnumerator::indexBatch(std::span<const EdgeId> batch) {
    endpoints_.clear();
    batchDegrees_.clear();
    for (const EdgeId edgeId : batch) {
        const Edge& edge = graph_.edge(edgeId);
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
    // The cover star builds: walking the batch in order, every edge that is
    // still uncovered puts both of its endpoints in.
    if (configuration_.seeds == VertexSeeds::Cover) {
        inCover_.assign(endpoints_.size(), 0U);
        for (const EdgeId edgeId : batch) {
            const Edge& edge = graph_.edge(edgeId);
            const auto uRank = static_cast<std::size_t>(
                endpointRankByVertex_[static_cast<std::size_t>(edge.u)]);
            const auto vRank = static_cast<std::size_t>(
                endpointRankByVertex_[static_cast<std::size_t>(edge.v)]);
            if (inCover_[uRank] == 0U && inCover_[vRank] == 0U) {
                inCover_[uRank] = 1U;
                inCover_[vRank] = 1U;
            }
        }
    }
    incidenceOffsets_.assign(endpoints_.size() + 1U, 0U);
    for (std::size_t index = 0; index < batchDegrees_.size(); ++index) {
        incidenceOffsets_[index + 1U] =
            incidenceOffsets_[index] + batchDegrees_[index];
    }
    incidenceNeighbors_.assign(incidenceOffsets_.back(), 0);
    incidenceFill_.assign(
        incidenceOffsets_.begin(), incidenceOffsets_.end() - 1);
    // Only the new-edge partition orders edges, so only it pays for their
    // positions.
    const bool positions = configuration_.branching != VertexBranching::Pivot;
    if (positions) {
        incidencePositions_.assign(incidenceOffsets_.back(), 0);
    }
    for (std::size_t position = 0; position < batch.size(); ++position) {
        const Edge& edge = graph_.edge(batch[position]);
        const auto uRank = static_cast<std::size_t>(
            endpointRankByVertex_[static_cast<std::size_t>(edge.u)]);
        const auto vRank = static_cast<std::size_t>(
            endpointRankByVertex_[static_cast<std::size_t>(edge.v)]);
        const std::size_t uSlot = incidenceFill_[uRank]++;
        const std::size_t vSlot = incidenceFill_[vRank]++;
        incidenceNeighbors_[uSlot] = edge.v;
        incidenceNeighbors_[vSlot] = edge.u;
        if (positions) {
            incidencePositions_[uSlot] = static_cast<int>(position);
            incidencePositions_[vSlot] = static_cast<int>(position);
        }
    }
}

bool VertexCliqueEnumerator::isSeed(int rank) const noexcept {
    return configuration_.seeds == VertexSeeds::Endpoints ||
        inCover_[static_cast<std::size_t>(rank)] != 0U;
}

std::span<const Vertex> VertexCliqueEnumerator::batchNeighbors(
    Vertex vertex) const noexcept {
    const int rank = endpointRankByVertex_[static_cast<std::size_t>(vertex)];
    if (rank < 0) {
        return {};
    }
    const auto index = static_cast<std::size_t>(rank);
    return std::span<const Vertex>(
        incidenceNeighbors_.data() + incidenceOffsets_[index],
        incidenceOffsets_[index + 1U] - incidenceOffsets_[index]);
}

std::span<const int> VertexCliqueEnumerator::batchPositions(
    Vertex vertex) const noexcept {
    const int rank = endpointRankByVertex_[static_cast<std::size_t>(vertex)];
    if (rank < 0) {
        return {};
    }
    const auto index = static_cast<std::size_t>(rank);
    return std::span<const int>(
        incidencePositions_.data() + incidenceOffsets_[index],
        incidenceOffsets_[index + 1U] - incidenceOffsets_[index]);
}

void VertexCliqueEnumerator::enumerateRoot(
    Vertex root,
    int rootRank,
    VertexSearchStats& stats,
    CliqueBatchOutput& output) {
    ++stats.roots;
    // P is the whole neighbourhood but the earlier seeds, which own every
    // clique holding them and so start in X.
    const auto startsInX = [&, this](int rank) {
        return rank >= 0 && rank < rootRank && isSeed(rank);
    };
    std::vector<Vertex> candidates;
    std::vector<Vertex> excluded;
    std::vector<Vertex> candidateEndpoints;
    const std::uint32_t epoch = beginEpoch();
    for (const AdjacencyEntry& entry : graph_.neighbors(root)) {
        const Vertex neighbor = entry.neighbor;
        const int rank =
            endpointRankByVertex_[static_cast<std::size_t>(neighbor)];
        if (startsInX(rank)) {
            excluded.push_back(neighbor);
            continue;
        }
        candidates.push_back(neighbor);
        epochByVertex_[static_cast<std::size_t>(neighbor)] = epoch;
        if (rank >= 0) {
            candidateEndpoints.push_back(neighbor);
        }
    }

    // The new-edge condition at the root, before its graph exists: one of
    // the root's own batch edges into P -- to an endpoint that is no earlier
    // seed -- or a batch edge between two vertices of P.
    const auto rootReachesBatchEdge = [&] {
        for (const Vertex partner : batchNeighbors(root)) {
            ++stats.batchEdgeChecks;
            if (!startsInX(
                    endpointRankByVertex_[static_cast<std::size_t>(partner)])) {
                return true;
            }
        }
        for (const Vertex endpoint : candidateEndpoints) {
            for (const Vertex neighbor : batchNeighbors(endpoint)) {
                ++stats.batchEdgeChecks;
                if (epochByVertex_[static_cast<std::size_t>(neighbor)] ==
                    epoch) {
                    return true;
                }
            }
        }
        return false;
    };
    if (!rootReachesBatchEdge()) {
        ++stats.prunedRoots;
        return;
    }

    Clique current{root};
    localGraph_.rebuild(
        graph_, current, candidates, excluded, stats.localGraph);
    localGraph_.setBits(candidateEndpoints, endpointMask_);

    if (configuration_.branching != VertexBranching::Pivot) {
        LocalVertexSet candidateSet = localGraph_.makeSet(std::move(candidates));
        LocalVertexSet excludedSet = localGraph_.makeSet(std::move(excluded));
        branchUnbound(
            current, candidateSet, excludedSet, true, false, stats, output);
        return;
    }

    // A vertex entering the clique completes a batch edge exactly when it is
    // the batch neighbour of a member.
    const auto updateState = [&, this](
        const LocalSearchGraph&,
        VertexSearchState parent,
        std::span<const Vertex> clique,
        Vertex vertex) {
        if (!parent.holdsBatchEdge) {
            for (const Vertex neighbor : batchNeighbors(vertex)) {
                ++stats.batchEdgeChecks;
                if (std::ranges::find(clique, neighbor) != clique.end()) {
                    return VertexSearchState{true};
                }
            }
        }
        return parent;
    };
    const auto pruneNode = [&, this](
        const LocalSearchGraph& graph,
        std::span<const Vertex> clique,
        const LocalVertexSet& candidateSet,
        const LocalVertexSet&,
        const VertexSearchState& state) {
        // The root passed the same test before its graph was built.
        if (state.holdsBatchEdge || clique.size() == 1U) {
            return false;
        }
        if (batchEdgeReachable(graph, clique, candidateSet, stats)) {
            return false;
        }
        ++stats.newEdgePrunes;
        return true;
    };
    const auto onLeaf = [&](
        std::span<const Vertex> clique,
        const VertexSearchState& state) {
        if (!state.holdsBatchEdge) {
            throw std::logic_error("a vertex-seeded leaf holds no batch edge");
        }
        ++stats.leaves;
        output.emit(clique);
        ++stats.emittedCliques;
    };
    bronKerboschLocal<false>(
        localGraph_,
        current,
        localGraph_.makeSet(std::move(candidates)),
        localGraph_.makeSet(std::move(excluded)),
        VertexSearchState{},
        updateState,
        NoChildRefinement{},
        pruneNode,
        onLeaf,
        stats);
}

bool VertexCliqueEnumerator::batchEdgeReachable(
    const LocalSearchGraph& graph,
    std::span<const Vertex> clique,
    const LocalVertexSet& candidates,
    VertexSearchStats& stats) const {
    // A batch edge from the clique into P.
    for (const Vertex member : clique) {
        for (const Vertex neighbor : batchNeighbors(member)) {
            ++stats.batchEdgeChecks;
            if (graph.contains(candidates, neighbor)) {
                return true;
            }
        }
    }
    // A batch edge inside P, which only the endpoints in the root's P can
    // carry.
    const std::span<const Vertex> locals = graph.localVertices();
    for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
        std::uint64_t bits = candidates.bits[word] & endpointMask_.bits[word];
        while (bits != 0U) {
            const Vertex endpoint = locals[word * 64U +
                static_cast<std::size_t>(std::countr_zero(bits))];
            bits &= bits - 1U;
            for (const Vertex neighbor : batchNeighbors(endpoint)) {
                ++stats.batchEdgeChecks;
                if (graph.contains(candidates, neighbor)) {
                    return true;
                }
            }
        }
    }
    return false;
}

std::size_t VertexCliqueEnumerator::collectNewEdges(
    std::span<const Vertex> clique,
    const LocalVertexSet& candidates,
    std::size_t limit,
    std::vector<NewEdge>& edges,
    VertexSearchStats& stats) const {
    const LocalSearchGraph& graph = localGraph_;
    edges.clear();
    std::size_t open = 0U;
    // Whether the edge just found takes the open ones past the limit.
    const auto add = [&](int position, Vertex first, Vertex second) {
        const bool skipped =
            joinsEarlierEdge(clique, candidates, first, second, position, stats);
        edges.push_back(NewEdge{position, first, second, skipped});
        return !skipped && ++open > limit;
    };
    // From the clique into P: R holds no batch edge of its own.
    for (const Vertex member : clique) {
        const std::span<const Vertex> neighbors = batchNeighbors(member);
        const std::span<const int> positions = batchPositions(member);
        for (std::size_t index = 0; index < neighbors.size(); ++index) {
            ++stats.batchEdgeChecks;
            if (graph.contains(candidates, neighbors[index]) &&
                add(positions[index], member, neighbors[index])) {
                return open;
            }
        }
    }
    // Inside P, which only the endpoints in the root's P can carry; each edge
    // once, from its smaller end.
    const std::span<const Vertex> locals = graph.localVertices();
    for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
        std::uint64_t bits = candidates.bits[word] & endpointMask_.bits[word];
        while (bits != 0U) {
            const Vertex endpoint = locals[word * 64U +
                static_cast<std::size_t>(std::countr_zero(bits))];
            bits &= bits - 1U;
            const std::span<const Vertex> neighbors = batchNeighbors(endpoint);
            const std::span<const int> positions = batchPositions(endpoint);
            for (std::size_t index = 0; index < neighbors.size(); ++index) {
                ++stats.batchEdgeChecks;
                const Vertex other = neighbors[index];
                if (other > endpoint && graph.contains(candidates, other) &&
                    add(positions[index], endpoint, other)) {
                    return open;
                }
            }
        }
    }
    // Complete: in batch order, the order its children are expanded in.
    std::ranges::sort(edges, {}, &NewEdge::position);
    return open;
}

bool VertexCliqueEnumerator::joinsEarlierEdge(
    std::span<const Vertex> clique,
    const LocalVertexSet& candidates,
    Vertex first,
    Vertex second,
    int position,
    VertexSearchStats& stats) const {
    // Only through an end the edge adds, as R holds no batch edge and the
    // edge is the only one between its ends; such a child owes nothing.
    for (const Vertex end : {first, second}) {
        if (!localGraph_.contains(candidates, end)) {
            continue;
        }
        const std::span<const Vertex> neighbors = batchNeighbors(end);
        const std::span<const int> positions = batchPositions(end);
        for (std::size_t index = 0;
             index < neighbors.size() && positions[index] < position; ++index) {
            ++stats.batchEdgeChecks;
            if (std::ranges::find(clique, neighbors[index]) != clique.end()) {
                return true;
            }
        }
    }
    return false;
}

void VertexCliqueEnumerator::excludeEarlierEdges(
    Vertex vertex,
    int threshold,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    VertexSearchStats& stats) const {
    const std::span<const Vertex> neighbors = batchNeighbors(vertex);
    const std::span<const int> positions = batchPositions(vertex);
    // Positions increase along the list, so the walk stops at the first edge
    // that is not earlier.
    for (std::size_t index = 0;
         index < neighbors.size() && positions[index] < threshold; ++index) {
        ++stats.batchEdgeChecks;
        const Vertex neighbor = neighbors[index];
        if (!localGraph_.contains(candidates, neighbor)) {
            continue;
        }
        localGraph_.erase(candidates, neighbor);
        localGraph_.insert(excluded, neighbor);
        ++stats.exclusionMoves;
    }
}

void VertexCliqueEnumerator::searchBound(
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    int threshold,
    VertexSearchStats& stats,
    CliqueBatchOutput& output) {
    const auto onLeaf = [&](
        std::span<const Vertex> clique,
        const EmptySearchState&) {
        ++stats.leaves;
        output.emit(clique);
        ++stats.emittedCliques;
    };
    if (threshold < 0) {
        bronKerboschLocalInPlace<false>(
            localGraph_, current, candidates, excluded, EmptySearchState{},
            NoStateUpdate{}, NoChildRefinement{}, NoBranchPruning{}, onLeaf,
            stats);
        return;
    }
    // A vertex entering the clique sends every candidate it reaches by an
    // earlier batch edge to X, as the edge-seeded baseline does.
    const auto refineChild = [&, this](
        const LocalSearchGraph&,
        Vertex vertex,
        LocalVertexSet& childCandidates,
        LocalVertexSet& childExcluded,
        std::span<const EdgeId>) {
        excludeEarlierEdges(
            vertex, threshold, childCandidates, childExcluded, stats);
    };
    bronKerboschLocalInPlace<false>(
        localGraph_, current, candidates, excluded, EmptySearchState{},
        NoStateUpdate{}, refineChild, NoBranchPruning{}, onLeaf, stats);
}

void VertexCliqueEnumerator::branchUnbound(
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    bool atRoot,
    bool pivotBelow,
    VertexSearchStats& stats,
    CliqueBatchOutput& output) {
    ++stats.recursiveCalls;
    const LocalSearchGraph& graph = localGraph_;
    if (unboundDepth_ == unboundFrames_.size()) {
        unboundFrames_.emplace_back();
    }
    UnboundFrame& frame = unboundFrames_[unboundDepth_];
    ++unboundDepth_;
    struct FrameRelease {
        std::size_t& depth;
        ~FrameRelease() { --depth; }
    } frameRelease{unboundDepth_};

    const bool decides = !pivotBelow &&
        (atRoot || configuration_.branching == VertexBranching::RuleUnbound);
    Vertex pivot = -1;
    bool newEdges = false;
    // |B| of the state, its pivot scanned into `pivot`.
    const auto pivotBranches = [&] {
        pivot = chooseLocalBkPivot(graph, candidates, excluded, stats);
        if (pivot < 0) {
            return candidates.size();
        }
        const std::span<const std::uint64_t> pivotRow = graph.rowAt(
            static_cast<std::size_t>(graph.indexOf(pivot)));
        std::size_t branches = 0U;
        for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
            branches += static_cast<std::size_t>(
                std::popcount(candidates.bits[word] & ~pivotRow[word]));
        }
        return branches;
    };
    if (!decides) {
        // As the baseline: cut when R u P holds no batch edge (a root passed
        // that test before its graph was built), and take the pivot.
        if (!atRoot && !batchEdgeReachable(graph, current, candidates, stats)) {
            ++stats.newEdgePrunes;
            return;
        }
    } else if (atRoot) {
        const std::size_t open = collectNewEdges(
            current, candidates, std::numeric_limits<std::size_t>::max(),
            frame.edges, stats);
        if (open == 0U) {
            ++stats.newEdgePrunes;
            return;
        }
        if (configuration_.branching == VertexBranching::NewEdges) {
            newEdges = true;
        } else if (open == 1U) {
            // One child either way unless a vertex of X covers all of P, a
            // state that owes nothing; as star, bound without a scan.
            newEdges = true;
            ++stats.unitBindings;
        } else {
            // No vertex covers more of P than its degree in the root's graph
            // less one, the root: when that bound keeps every vertex at or
            // under |P| - |W|, the new edges win without a scan.
            const long long mark = static_cast<long long>(candidates.size()) -
                static_cast<long long>(open);
            bool dense = false;
            for (const LocalVertexSet* set : {&excluded, &candidates}) {
                for (const Vertex vertex : set->vertices) {
                    const auto local =
                        static_cast<std::size_t>(graph.indexOf(vertex));
                    if (static_cast<long long>(graph.localDegree(local)) - 1 >
                        mark) {
                        dense = true;
                        break;
                    }
                }
                if (dense) {
                    break;
                }
            }
            if (dense) {
                newEdges = open <= pivotBranches();
            } else {
                newEdges = true;
                ++stats.degreeDecisions;
            }
        }
        ++(newEdges ? stats.rootEdgeChoices : stats.rootPivotChoices);
    } else {
        // Below a root most states take the pivot, so W is read only as far
        // as the comparison needs: up to a second open edge before the pivot
        // scan, then past |B| at most.  A list that stopped short is never
        // expanded.
        const std::size_t first =
            collectNewEdges(current, candidates, 1U, frame.edges, stats);
        if (first == 0U) {
            ++stats.newEdgePrunes;
            return;
        }
        if (first == 1U) {
            newEdges = true;
            ++stats.unitBindings;
        } else {
            const std::size_t branches = pivotBranches();
            newEdges = collectNewEdges(
                current, candidates, branches, frame.edges, stats) <= branches;
        }
        ++(newEdges ? stats.unboundEdgeChoices : stats.unboundPivotChoices);
    }

    LocalVertexSet& childCandidates = frame.childCandidates;
    LocalVertexSet& childExcluded = frame.childExcluded;
    if (newEdges) {
        // One child per edge of W in batch order, independent of one another:
        // each takes the cliques whose first batch edge it is.
        for (const NewEdge& edge : frame.edges) {
            if (edge.skipped) {
                ++stats.skippedEdgeChildren;
                continue;
            }
            ++stats.edgeChildren;
            std::size_t added = 0U;
            for (const Vertex end : {edge.first, edge.second}) {
                if (!graph.contains(candidates, end)) {
                    continue;
                }
                if (added == 0U) {
                    graph.intersectInto(candidates, end, childCandidates);
                    graph.intersectInto(excluded, end, childExcluded);
                } else {
                    graph.intersectInto(childCandidates, end, childCandidates);
                    graph.intersectInto(childExcluded, end, childExcluded);
                }
                current.push_back(end);
                ++added;
            }
            for (const Vertex member : current) {
                excludeEarlierEdges(
                    member, edge.position, childCandidates, childExcluded, stats);
            }
            searchBound(
                current, childCandidates, childExcluded, edge.position, stats,
                output);
            current.resize(current.size() - added);
        }
        return;
    }

    // The pivot partition.  Below a root-scope root that took it, every state
    // takes it again without asking.
    const bool pivotOnlyBelow = pivotBelow ||
        configuration_.branching == VertexBranching::RuleRoot;
    if (pivot < 0) {
        pivot = chooseLocalBkPivot(graph, candidates, excluded, stats);
    }
    std::vector<Vertex>& extensions = frame.extensions;
    if (pivot < 0) {
        extensions.assign(candidates.vertices.begin(), candidates.vertices.end());
    } else {
        graph.differenceVerticesInto(candidates, pivot, extensions);
    }
    for (const Vertex vertex : extensions) {
        if (!graph.contains(candidates, vertex)) {
            continue;
        }
        graph.intersectInto(candidates, vertex, childCandidates);
        graph.intersectInto(excluded, vertex, childExcluded);
        // The child holds a batch edge exactly when the vertex is the batch
        // neighbour of a member.
        bool bound = false;
        for (const Vertex neighbor : batchNeighbors(vertex)) {
            ++stats.batchEdgeChecks;
            if (std::ranges::find(current, neighbor) != current.end()) {
                bound = true;
                break;
            }
        }
        current.push_back(vertex);
        if (bound) {
            searchBound(
                current, childCandidates, childExcluded, -1, stats, output);
        } else {
            branchUnbound(
                current, childCandidates, childExcluded, false, pivotOnlyBelow,
                stats, output);
        }
        current.pop_back();
        graph.moveToExcludedBits(candidates, excluded, vertex);
    }
}

} // namespace immce
