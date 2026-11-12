#include "new_clique_search.hpp"

#include <stdexcept>

namespace immce {

NewCliqueEnumerator::NewCliqueEnumerator(const IncrementalGraph& graph)
    : graph_(graph),
      batchIndex_(graph),
      localGraph_(graph.vertexCount()) {}

std::vector<Clique> NewCliqueEnumerator::enumerateBatch(
    std::span<const EdgeId> batch,
    NewCliqueSearchStats& stats) {
    CliqueBatchOutput output(true);
    enumerateBatchInto(batch, stats, output);
    return std::move(output).takeCliques();
}

long long NewCliqueEnumerator::countBatch(
    std::span<const EdgeId> batch,
    NewCliqueSearchStats& stats) {
    CliqueBatchOutput output(false);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

long long NewCliqueEnumerator::visitBatch(
    std::span<const EdgeId> batch,
    NewCliqueSearchStats& stats,
    CliqueVisitor visitor) {
    if (!visitor) {
        throw std::invalid_argument("new-clique visitor must not be empty");
    }
    CliqueBatchOutput output(false, visitor);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

void NewCliqueEnumerator::enumerateBatchInto(
    std::span<const EdgeId> batch,
    NewCliqueSearchStats& stats,
    CliqueBatchOutput& output) {
    auto batchScope = batchIndex_.begin(batch);
    for (std::size_t position = 0; position < batch.size(); ++position) {
        ++stats.seeds;
        enumerateSeed(batch[position], static_cast<int>(position), stats, output);
    }
}

void NewCliqueEnumerator::enumerateSeed(
    EdgeId seedId,
    int seedPosition,
    NewCliqueSearchStats& stats,
    CliqueBatchOutput& output) {
    const Edge& seed = graph_.edge(seedId);
    Clique root{seed.u, seed.v};
    std::vector<Vertex> candidates;
    std::vector<Vertex> excluded;
    graph_.forEachCommonNeighbor(
        seed.u, seed.v,
        [&](Vertex vertex, EdgeId edgeFromU, EdgeId edgeFromV) {
            if (isEarlierBatchEdge(edgeFromU, seedPosition) ||
                isEarlierBatchEdge(edgeFromV, seedPosition)) {
                excluded.push_back(vertex);
            } else {
                candidates.push_back(vertex);
            }
        });

    localGraph_.rebuild(graph_, root, candidates, excluded, stats.localGraph);
    const auto ownershipRefiner = [&, this](
        const LocalSearchGraph& graph,
        Vertex,
        LocalVertexSet& childCandidates,
        LocalVertexSet& childExcluded,
        std::span<const EdgeId> connectingEdgeIds) {
        refineChildOwnership(
            graph,
            seedPosition,
            childCandidates,
            childExcluded,
            connectingEdgeIds,
            stats);
    };
    const auto onLeaf = [&](
        std::span<const Vertex> clique,
        const EmptySearchState&) {
        ++stats.leaves;
        output.emit(clique);
        ++stats.emittedCliques;
    };
    bronKerboschLocal<true>(
        localGraph_,
        root,
        localGraph_.makeSet(std::move(candidates)),
        localGraph_.makeSet(std::move(excluded)),
        EmptySearchState{},
        NoStateUpdate{},
        ownershipRefiner,
        NoBranchPruning{},
        onLeaf,
        stats);
}

void NewCliqueEnumerator::refineChildOwnership(
    const LocalSearchGraph& localGraph,
    int seedPosition,
    LocalVertexSet& childCandidates,
    LocalVertexSet& childExcluded,
    std::span<const EdgeId> connectingEdgeIds,
    NewCliqueSearchStats& stats) const {
    if (childCandidates.vertices.size() != connectingEdgeIds.size()) {
        throw std::logic_error("candidate/edge ownership vectors are misaligned");
    }
    std::vector<Vertex> moved;
    for (std::size_t index = 0; index < childCandidates.vertices.size(); ++index) {
        if (isEarlierBatchEdge(connectingEdgeIds[index], seedPosition)) {
            moved.push_back(childCandidates.vertices[index]);
        }
    }
    for (Vertex vertex : moved) {
        localGraph.erase(childCandidates, vertex);
        localGraph.insert(childExcluded, vertex);
        ++stats.ownershipMoves;
    }
}

bool NewCliqueEnumerator::isEarlierBatchEdge(
    EdgeId edgeIdValue,
    int seedPosition) const noexcept {
    return edgeIdValue != invalidEdgeId &&
        batchIndex_.edgeRank(edgeIdValue) >= 0 &&
        batchIndex_.edgeRank(edgeIdValue) < seedPosition;
}

} // namespace immce
