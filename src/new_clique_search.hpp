#pragma once

#include "batch_index.hpp"
#include "local_bk.hpp"

#include <span>
#include <vector>

namespace immce {

struct NewCliqueSearchStats {
    long long seeds = 0;
    long long recursiveCalls = 0;
    long long pivotCalls = 0;
    long long pivotCandidatesEvaluated = 0;
    long long ownershipMoves = 0;
    long long leaves = 0;
    long long emittedCliques = 0;
    LocalGraphStats localGraph;
};

// Edge-seeded enumeration: the baseline of the new-clique search.
//
// Every batch edge e = (u, v) seeds a search of its own. A clique holding e
// lies inside the support S_e = {u, v} u (N(u) n N(v)), and so does every
// vertex that extends it, so the maximal cliques of G[S_e] holding e are the
// graph's maximal cliques holding e. A clique holding several batch edges
// belongs to the earliest of them in the batch: a common neighbour joined to
// u or v by an earlier batch edge starts in X, and a candidate joined by one
// to a vertex entering the clique moves to X with it, so it is never branched
// on and no clique is reported twice.
class NewCliqueEnumerator {
public:
    // The ownership test needs, for every child, the edge between the vertex
    // entering the clique and each candidate; the seed graphs read it back
    // from the global adjacency.
    explicit NewCliqueEnumerator(const IncrementalGraph& graph);

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        NewCliqueSearchStats& stats);
    [[nodiscard]] long long countBatch(
        std::span<const EdgeId> batch,
        NewCliqueSearchStats& stats);
    [[nodiscard]] long long visitBatch(
        std::span<const EdgeId> batch,
        NewCliqueSearchStats& stats,
        CliqueVisitor visitor);

private:
    void enumerateBatchInto(
        std::span<const EdgeId> batch,
        NewCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    void enumerateSeed(
        EdgeId seedId,
        int seedPosition,
        NewCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    void refineChildOwnership(
        const LocalSearchGraph& localGraph,
        int seedPosition,
        LocalVertexSet& childCandidates,
        LocalVertexSet& childExcluded,
        std::span<const EdgeId> connectingEdgeIds,
        NewCliqueSearchStats& stats) const;
    [[nodiscard]] bool isEarlierBatchEdge(EdgeId edgeId, int seedPosition) const noexcept;

    const IncrementalGraph& graph_;
    BatchIndex batchIndex_;
    // One local graph for the enumerator's lifetime, as the star enumerator
    // keeps its own: its constructor fills a |V|-sized index table (7-9 MB on
    // the 2M-vertex graphs), which a graph built per batch pays every batch.
    LocalSearchGraph localGraph_;
};

} // namespace immce
