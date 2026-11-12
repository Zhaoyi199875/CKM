#pragma once

#include "local_search_graph.hpp"

#include <concepts>
#include <span>
#include <utility>
#include <vector>

namespace immce {

struct EmptySearchState {};

struct NoStateUpdate {
    [[nodiscard]] EmptySearchState operator()(
        const LocalSearchGraph&,
        EmptySearchState parent,
        std::span<const Vertex>,
        Vertex) const noexcept {
        return parent;
    }
};

struct NoChildRefinement {
    void operator()(
        const LocalSearchGraph&,
        Vertex,
        LocalVertexSet&,
        LocalVertexSet&,
        std::span<const EdgeId>) const noexcept {}
};

struct NoBranchPruning {
    template <typename State>
    [[nodiscard]] bool operator()(
        const LocalSearchGraph&,
        std::span<const Vertex>,
        const LocalVertexSet&,
        const LocalVertexSet&,
        const State&) const noexcept {
        return false;
    }
};

template <typename Updater, typename State>
concept LocalBkStateUpdater = requires(
    const Updater& updater,
    const LocalSearchGraph& graph,
    State state,
    std::span<const Vertex> clique,
    Vertex vertex) {
    { updater(graph, std::move(state), clique, vertex) } -> std::same_as<State>;
};

template <typename Refiner>
concept LocalBkChildRefiner = requires(
    const Refiner& refiner,
    const LocalSearchGraph& graph,
    Vertex vertex,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    std::span<const EdgeId> edgeIds) {
    { refiner(graph, vertex, candidates, excluded, edgeIds) } -> std::same_as<void>;
};

template <typename Pruner, typename State>
concept LocalBkBranchPruner = requires(
    const Pruner& pruner,
    const LocalSearchGraph& graph,
    std::span<const Vertex> clique,
    const LocalVertexSet& candidates,
    const LocalVertexSet& excluded,
    State& state) {
    { pruner(graph, clique, candidates, excluded, state) } -> std::same_as<bool>;
};

template <typename Handler, typename State>
concept LocalBkLeafHandler = requires(
    const Handler& handler,
    std::span<const Vertex> clique,
    const State& state) {
    { handler(clique, state) } -> std::same_as<void>;
};

template <typename Stats>
concept LocalBkStatistics = requires(Stats& stats) {
    ++stats.recursiveCalls;
    ++stats.pivotCalls;
    ++stats.pivotCandidatesEvaluated;
};

// Statistics that also want the bitmap width the search pays for: the words
// per row of the graph under every node, and under every pivot candidate
// scored. Every set operation of a node costs one pass over that width, so
// these two sums are the search's word work, which node counts alone hide
// when two searches run on graphs of different widths.
template <typename Stats>
concept LocalBkWordStatistics = requires(Stats& stats) {
    stats.nodeWords += 1;
    stats.pivotWords += 1;
};

// The default descent into a child: the same graph, the same sets.
struct DirectDescent {
    template <typename Recurse>
    void operator()(
        const LocalSearchGraph& graph,
        const Clique&,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        const Recurse& recurse) const {
        recurse(graph, candidates, excluded);
    }
};

// The vertex of X u P covering the most of P, ties to the smaller vertex. X is
// scored first and each set in increasing vertex order, and the scan stops as
// soon as nothing left can change the answer:
//  * a vertex of X covering all of P, which no other vertex can beat;
//  * once the best score is |P| - 1 -- all a vertex of P can cover, being no
//    neighbour of itself -- every vertex of P still to come can at most tie,
//    and loses the tie by being larger than the pivot already held. So the
//    scan of P stops at the first vertex past the pivot, which is at once
//    when a vertex of P itself reaches |P| - 1.
// Either way the pivot is the one the full scan returns.
template <LocalBkStatistics Stats>
[[nodiscard]] Vertex chooseLocalBkPivot(
    const LocalSearchGraph& graph,
    const LocalVertexSet& candidates,
    const LocalVertexSet& excluded,
    Stats& stats) {
    ++stats.pivotCalls;
    Vertex pivot = -1;
    int maximumIntersection = -1;
    const auto consider = [&](Vertex vertex) {
        ++stats.pivotCandidatesEvaluated;
        if constexpr (LocalBkWordStatistics<Stats>) {
            stats.pivotWords += static_cast<long long>(graph.wordsPerRow());
        }
        const int count = graph.intersectionCount(candidates, vertex);
        if (count > maximumIntersection ||
            (count == maximumIntersection && (pivot < 0 || vertex < pivot))) {
            maximumIntersection = count;
            pivot = vertex;
        }
        return count;
    };
    const auto candidateCount = static_cast<int>(candidates.size());
    for (Vertex vertex : excluded.vertices) {
        if (consider(vertex) == candidateCount) {
            return pivot;
        }
    }
    for (Vertex vertex : candidates.vertices) {
        if (maximumIntersection == candidateCount - 1 && vertex > pivot) {
            return pivot;
        }
        consider(vertex);
    }
    return pivot;
}

// `descend` sees every child once its sets are final and decides which graph
// the child's subtree runs on: it calls the recursion it is handed with the
// graph and the sets to use, which a re-indexing descent re-expresses in a
// narrower universe of its own. DirectDescent keeps the search exactly as it
// was.
template <
    bool TrackConnectingEdgeIds,
    typename State,
    typename StateUpdater,
    typename ChildRefiner,
    typename BranchPruner,
    typename LeafHandler,
    typename Descender,
    LocalBkStatistics Stats>
requires LocalBkStateUpdater<StateUpdater, State> &&
         LocalBkChildRefiner<ChildRefiner> &&
         LocalBkBranchPruner<BranchPruner, State> &&
         LocalBkLeafHandler<LeafHandler, State>
void bronKerboschLocalInPlace(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    State state,
    const StateUpdater& updateState,
    const ChildRefiner& refineChild,
    const BranchPruner& pruneBranch,
    const LeafHandler& onLeaf,
    const Descender& descend,
    Stats& stats) {
    ++stats.recursiveCalls;
    if constexpr (LocalBkWordStatistics<Stats>) {
        stats.nodeWords += static_cast<long long>(graph.wordsPerRow());
    }
    if (pruneBranch(graph, current, candidates, excluded, state)) {
        return;
    }
    if (candidates.empty() && excluded.empty()) {
        onLeaf(current, state);
        return;
    }

    const Vertex pivot = chooseLocalBkPivot(graph, candidates, excluded, stats);
    std::vector<Vertex> extensions = pivot < 0
        ? candidates.vertices
        : graph.differenceVerticesFromNeighbors(candidates, pivot);
    LocalVertexSet childCandidates;
    LocalVertexSet childExcluded;
    std::vector<EdgeId> connectingEdgeIds;
    for (Vertex vertex : extensions) {
        if (!graph.contains(candidates, vertex)) {
            continue;
        }

        State childState = updateState(graph, state, current, vertex);
        current.push_back(vertex);
        if constexpr (TrackConnectingEdgeIds) {
            graph.intersectCandidatesWithNeighborEdges(
                candidates, vertex, childCandidates, connectingEdgeIds);
        } else {
            graph.intersectInto(candidates, vertex, childCandidates);
        }
        graph.intersectInto(excluded, vertex, childExcluded);
        refineChild(
            graph,
            vertex,
            childCandidates,
            childExcluded,
            connectingEdgeIds);
        descend(
            graph,
            current,
            childCandidates,
            childExcluded,
            [&](const LocalSearchGraph& childGraph,
                LocalVertexSet& descendantCandidates,
                LocalVertexSet& descendantExcluded) {
                bronKerboschLocalInPlace<TrackConnectingEdgeIds>(
                    childGraph,
                    current,
                    descendantCandidates,
                    descendantExcluded,
                    std::move(childState),
                    updateState,
                    refineChild,
                    pruneBranch,
                    onLeaf,
                    descend,
                    stats);
            });
        current.pop_back();

        // Only the bits: the pivot and the branch list are already read, and
        // every later child is cut from the bits. On return the two sets'
        // vertex lists are stale, which no caller reads.
        graph.moveToExcludedBits(candidates, excluded, vertex);
    }
}

template <
    bool TrackConnectingEdgeIds,
    typename State,
    typename StateUpdater,
    typename ChildRefiner,
    typename BranchPruner,
    typename LeafHandler,
    LocalBkStatistics Stats>
requires LocalBkStateUpdater<StateUpdater, State> &&
         LocalBkChildRefiner<ChildRefiner> &&
         LocalBkBranchPruner<BranchPruner, State> &&
         LocalBkLeafHandler<LeafHandler, State>
void bronKerboschLocalInPlace(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    State state,
    const StateUpdater& updateState,
    const ChildRefiner& refineChild,
    const BranchPruner& pruneBranch,
    const LeafHandler& onLeaf,
    Stats& stats) {
    bronKerboschLocalInPlace<TrackConnectingEdgeIds>(
        graph,
        current,
        candidates,
        excluded,
        std::move(state),
        updateState,
        refineChild,
        pruneBranch,
        onLeaf,
        DirectDescent{},
        stats);
}

template <
    bool TrackConnectingEdgeIds,
    typename State,
    typename StateUpdater,
    typename ChildRefiner,
    typename BranchPruner,
    typename LeafHandler,
    LocalBkStatistics Stats>
requires LocalBkStateUpdater<StateUpdater, State> &&
         LocalBkChildRefiner<ChildRefiner> &&
         LocalBkBranchPruner<BranchPruner, State> &&
         LocalBkLeafHandler<LeafHandler, State>
void bronKerboschLocal(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet candidates,
    LocalVertexSet excluded,
    State state,
    const StateUpdater& updateState,
    const ChildRefiner& refineChild,
    const BranchPruner& pruneBranch,
    const LeafHandler& onLeaf,
    Stats& stats) {
    bronKerboschLocalInPlace<TrackConnectingEdgeIds>(
        graph,
        current,
        candidates,
        excluded,
        std::move(state),
        updateState,
        refineChild,
        pruneBranch,
        onLeaf,
        stats);
}

} // namespace immce
