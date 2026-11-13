#pragma once

#include "batch_index.hpp"
#include "local_bk.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace ckm {

struct VertexSeedStats {
    long long batches = 0;
    long long batchEdges = 0;
    // One root per batch endpoint.  A root whose neighbourhood holds no batch
    // edge it may own is cut before its graph is built.
    long long roots = 0;
    long long prunedRoots = 0;
    long long recursiveCalls = 0;
    long long pivotCalls = 0;
    long long pivotCandidatesEvaluated = 0;
    // Search nodes cut by the live test, and the batch incidence entries
    // read deciding it and tracking whether the clique holds a batch edge
    // yet.
    long long liveTestPrunes = 0;
    long long batchEdgeChecks = 0;
    long long leaves = 0;
    long long emittedCliques = 0;
    LocalGraphStats localGraph;
};

// VertexSeed (Section III-A): the new-clique search seeded by vertices.
//
// Every endpoint s of a batch edge seeds a search of its own: the classic
// Bron--Kerbosch with Tomita's pivot, rooted at R = {s} over the closed
// neighbourhood N[s]. Every clique holding s lies in N[s], and so does every
// vertex that extends one, so the maximal cliques of G[N[s]] holding s are
// the graph's maximal cliques holding s.
//
// Ownership. Endpoints rank by their first appearance in the batch, and a new
// clique belongs to the earliest endpoint it holds: the root's earlier
// endpoints start in X, so no clique is reported twice and nothing below the
// root checks ownership again.
//
// Live test. A new maximal clique is one holding a batch edge, and most
// cliques around an endpoint hold none. The search carries whether R holds a
// batch edge; until it does, a node is cut as soon as R u P holds no batch
// edge with an end in P, since nothing below it can qualify any more. The
// root makes the same test before its graph is built.
class VertexSeed {
public:
    explicit VertexSeed(const IncrementalGraph& graph);

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        VertexSeedStats& stats);
    [[nodiscard]] long long countBatch(
        std::span<const EdgeId> batch,
        VertexSeedStats& stats);
    [[nodiscard]] long long visitBatch(
        std::span<const EdgeId> batch,
        VertexSeedStats& stats,
        CliqueVisitor visitor);

private:
    void enumerateBatchInto(
        std::span<const EdgeId> batch,
        VertexSeedStats& stats,
        CliqueBatchOutput& output);
    // The batch endpoints in first-appearance order and each one's batch
    // neighbours.
    void indexBatch(std::span<const EdgeId> batch);
    void enumerateRoot(
        Vertex root,
        int rootRank,
        VertexSeedStats& stats,
        CliqueBatchOutput& output);
    // Empty for a vertex that is not a batch endpoint.
    [[nodiscard]] std::span<const Vertex> batchNeighbors(
        Vertex vertex) const noexcept;
    // Whether R u P still holds a batch edge with an end in P: the live test.
    [[nodiscard]] bool batchEdgeReachable(
        const LocalSearchGraph& graph,
        std::span<const Vertex> clique,
        const LocalVertexSet& candidates,
        VertexSeedStats& stats) const;
    [[nodiscard]] std::uint32_t beginEpoch();

    const IncrementalGraph& graph_;
    BatchIndex batchIndex_;
    // One local graph for the enumerator's lifetime, as CoverSeed keeps its
    // own.
    LocalSearchGraph localGraph_;
    std::vector<Vertex> endpoints_;
    std::vector<int> endpointRankByVertex_;
    std::vector<std::size_t> batchDegrees_;
    // Batch incidence: endpoint rank -> its batch neighbours, in batch order.
    std::vector<std::size_t> incidenceOffsets_;
    std::vector<std::size_t> incidenceFill_;
    std::vector<Vertex> incidenceNeighbors_;
    // The endpoints in the current root's P -- its later endpoints -- the only
    // vertices of P that can carry a batch edge, as a mask over its graph.
    LocalVertexSet endpointMask_;
    // Marks the current root's P, for the test made before its graph exists.
    std::vector<std::uint32_t> epochByVertex_;
    std::uint32_t epoch_ = 0U;
};

} // namespace ckm
