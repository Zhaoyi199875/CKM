#pragma once

#include "local_search_graph.hpp"

#include <concepts>
#include <span>
#include <vector>

namespace ckm {

template <typename Refiner>
concept LocalBkChildRefiner = requires(
    const Refiner& refiner,
    const LocalSearchGraph& graph,
    Vertex vertex,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded) {
    { refiner(graph, vertex, candidates, excluded) } -> std::same_as<void>;
};

template <typename Handler>
concept LocalBkLeafHandler = requires(
    const Handler& handler,
    std::span<const Vertex> clique) {
    { handler(clique) } -> std::same_as<void>;
};

template <typename Stats>
concept LocalBkStatistics = requires(Stats& stats) {
    ++stats.recursiveCalls;
    ++stats.pivotCalls;
    ++stats.pivotCandidatesEvaluated;
};

template <typename Stats>
concept LocalBkWordStatistics = requires(Stats& stats) {
    stats.nodeWords += 1;
    stats.pivotWords += 1;
};

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

template <
    typename ChildRefiner,
    typename LeafHandler,
    typename Descender,
    LocalBkStatistics Stats>
requires LocalBkChildRefiner<ChildRefiner> &&
         LocalBkLeafHandler<LeafHandler>
void bronKerboschLocalInPlace(
    const LocalSearchGraph& graph,
    Clique& current,
    LocalVertexSet& candidates,
    LocalVertexSet& excluded,
    const ChildRefiner& refineChild,
    const LeafHandler& onLeaf,
    const Descender& descend,
    Stats& stats) {
    ++stats.recursiveCalls;
    if constexpr (LocalBkWordStatistics<Stats>) {
        stats.nodeWords += static_cast<long long>(graph.wordsPerRow());
    }
    if (candidates.empty() && excluded.empty()) {
        onLeaf(current);
        return;
    }

    const Vertex pivot = chooseLocalBkPivot(graph, candidates, excluded, stats);
    std::vector<Vertex> extensions = pivot < 0
        ? candidates.vertices
        : graph.differenceVerticesFromNeighbors(candidates, pivot);
    LocalVertexSet childCandidates;
    LocalVertexSet childExcluded;
    for (Vertex vertex : extensions) {
        if (!graph.contains(candidates, vertex)) {
            continue;
        }

        current.push_back(vertex);
        graph.intersectInto(candidates, vertex, childCandidates);
        graph.intersectInto(excluded, vertex, childExcluded);
        refineChild(graph, vertex, childCandidates, childExcluded);
        descend(
            graph,
            current,
            childCandidates,
            childExcluded,
            [&](const LocalSearchGraph& childGraph,
                LocalVertexSet& descendantCandidates,
                LocalVertexSet& descendantExcluded) {
                bronKerboschLocalInPlace(
                    childGraph,
                    current,
                    descendantCandidates,
                    descendantExcluded,
                    refineChild,
                    onLeaf,
                    descend,
                    stats);
            });
        current.pop_back();

        graph.moveToExcludedBits(candidates, excluded, vertex);
    }
}

}
