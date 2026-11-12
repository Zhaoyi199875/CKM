#pragma once

#include "batch_index.hpp"
#include "local_bk.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <span>
#include <vector>

namespace immce {

// How a vertex-seeded search branches at a state whose clique holds no batch
// edge yet -- the branching ablation of 2026-10-21.  The two partitions of
// such a state are the ones star chooses between: by the pivot (one child per
// vertex of P \ N(q)) or by the new edges (one child per batch edge the state
// can own, each searched below as the edge-seeded baseline searches its
// edge).
enum class VertexBranching {
    // The pivot at every state, pruned by the new-edge condition: the
    // baseline, and the default.
    Pivot,
    // The new edges at the root, so every child already holds a batch edge.
    NewEdges,
    // The branch-count rule at the root -- the new edges when they open no
    // more children than the pivot, |W| <= |B|, the pivot otherwise and then
    // all the way down.
    RuleRoot,
    // The rule at every state whose clique holds no batch edge yet.
    RuleUnbound,
};

[[nodiscard]] const char* vertexBranchingName(VertexBranching branching) noexcept;

// Which batch endpoints seed a search -- the seed-cover ablation of
// 2026-11-01, which separates the seed cover from the edge assignment.
enum class VertexSeeds {
    // Every batch endpoint: the baseline, and the default.
    Endpoints,
    // Only the endpoints in the vertex cover of the batch edges that star
    // builds, with no edge assignment: each seed still owes every new clique
    // that holds it and no earlier seed.
    Cover,
};

[[nodiscard]] const char* vertexSeedsName(VertexSeeds seeds) noexcept;

struct VertexSearchConfiguration {
    VertexBranching branching = VertexBranching::Pivot;
    VertexSeeds seeds = VertexSeeds::Endpoints;
};

struct VertexSearchStats {
    long long batches = 0;
    long long batchEdges = 0;
    // One root per seed (every batch endpoint by default). A root whose
    // neighbourhood holds no batch edge it may own is cut before its graph is
    // built.
    long long roots = 0;
    long long prunedRoots = 0;
    long long recursiveCalls = 0;
    long long pivotCalls = 0;
    long long pivotCandidatesEvaluated = 0;
    // Search nodes cut by the new-edge condition, and the batch incidence
    // entries read deciding it and tracking whether the clique holds a batch
    // edge yet.
    long long newEdgePrunes = 0;
    long long batchEdgeChecks = 0;
    long long leaves = 0;
    long long emittedCliques = 0;
    // The branching ablation (VertexBranching).  The states that chose between
    // the two partitions, at the roots and below them, and what they took;
    // the root choices the degree bound made without a pivot scan, and the
    // states with a single new edge left, bound without one.  The children the
    // new edges opened, and those skipped because the clique with the edge
    // already holds an earlier batch edge (its cliques belong to that edge's
    // child); and the vertices moved to X for closing an earlier batch edge,
    // the edge-seeded exclusion.
    long long rootEdgeChoices = 0;
    long long rootPivotChoices = 0;
    long long unboundEdgeChoices = 0;
    long long unboundPivotChoices = 0;
    long long degreeDecisions = 0;
    long long unitBindings = 0;
    long long edgeChildren = 0;
    long long skippedEdgeChildren = 0;
    long long exclusionMoves = 0;
    LocalGraphStats localGraph;
};

// Vertex-seeded enumeration: the second baseline of the new-clique search.
//
// Every endpoint s of a batch edge seeds a search of its own: the classic
// Bron--Kerbosch with Tomita's pivot, rooted at R = {s} over the closed
// neighbourhood N[s]. Every clique holding s lies in N[s], and so does every
// vertex that extends one, so the maximal cliques of G[N[s]] holding s are
// the graph's maximal cliques holding s.
//
// Seeds (VertexSeeds). Every endpoint, or only those in the vertex cover of
// the batch edges: every new clique holds a batch edge and so a cover vertex,
// which is all the ownership below needs.
//
// Ownership. Endpoints rank by their first appearance in the batch, and a new
// clique belongs to the earliest seed it holds: the root's earlier seeds start
// in X, so no clique is reported twice and nothing below the root checks
// ownership again. An earlier endpoint outside the cover is no seed and stays
// in P.
//
// New-edge condition. A new maximal clique is one holding a batch edge, and
// most cliques around an endpoint hold none. The search carries whether R
// holds a batch edge; until it does, a node is cut as soon as R u P holds no
// batch edge with an end in P, since nothing below it can qualify any more.
// The root makes the same test before its graph is built.
//
// Branching (VertexBranching).  A state whose clique holds no batch edge yet
// owes the cliques of R u P that hold one, and the batch edges they can hold
// are W, those inside R u P (each with an end in P, as R holds none).  Every
// owed clique has a first of them in batch order, so one child per edge of W
// -- R joined by the edge, the edge-seeded baseline's search below it, every
// earlier batch edge excluded -- partitions what the state owes exactly as
// the pivot does; a child whose clique would already hold an earlier edge of
// W owes nothing and is skipped.  Each child of either partition holds at
// most |P| - |B| candidates, B = P \ N(q), so star's rule |W| <= |B| keeps
// the Tomita bound.  The root keeps its ownership by seeds either way: W
// includes batch edges that do not touch the root, which the cliques it owns
// may hold alone.
class VertexCliqueEnumerator {
public:
    explicit VertexCliqueEnumerator(
        const IncrementalGraph& graph,
        VertexSearchConfiguration configuration = {});

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        VertexSearchStats& stats);
    [[nodiscard]] long long countBatch(
        std::span<const EdgeId> batch,
        VertexSearchStats& stats);
    [[nodiscard]] long long visitBatch(
        std::span<const EdgeId> batch,
        VertexSearchStats& stats,
        CliqueVisitor visitor);

private:
    void enumerateBatchInto(
        std::span<const EdgeId> batch,
        VertexSearchStats& stats,
        CliqueBatchOutput& output);
    // The batch endpoints in first-appearance order, each one's batch
    // neighbours, and with cover seeds the cover.
    void indexBatch(std::span<const EdgeId> batch);
    // Whether the endpoint of this rank seeds a search.
    [[nodiscard]] bool isSeed(int rank) const noexcept;
    void enumerateRoot(
        Vertex root,
        int rootRank,
        VertexSearchStats& stats,
        CliqueBatchOutput& output);
    // Empty for a vertex that is not a batch endpoint.
    [[nodiscard]] std::span<const Vertex> batchNeighbors(
        Vertex vertex) const noexcept;
    // The batch positions of those edges, in the same (increasing) order.
    [[nodiscard]] std::span<const int> batchPositions(
        Vertex vertex) const noexcept;
    // Whether R u P still holds a batch edge with an end in P.
    [[nodiscard]] bool batchEdgeReachable(
        const LocalSearchGraph& graph,
        std::span<const Vertex> clique,
        const LocalVertexSet& candidates,
        VertexSearchStats& stats) const;
    [[nodiscard]] std::uint32_t beginEpoch();

    // The branching ablation.  One batch edge of W, by its ends; `skipped`
    // when R joined by it would already hold an earlier batch edge.
    struct NewEdge {
        int position = 0;
        Vertex first = 0;
        Vertex second = 0;
        bool skipped = false;
    };
    // A state whose clique holds no batch edge yet, and the root's own when
    // `atRoot`: prunes it when W is empty, otherwise expands the partition
    // the branching takes.  `pivotBelow` makes it and every state below take
    // the pivot without asking (a root-scope rule whose root took the pivot).
    void branchUnbound(
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool atRoot,
        bool pivotBelow,
        VertexSearchStats& stats,
        CliqueBatchOutput& output);
    // W of the state into `edges`, and how many of its edges are not
    // skipped.  Reading stops once more than `limit` are not, leaving the
    // list short and returning `limit` + 1: below a root, that is all the
    // rule needs to take the pivot.  A complete list is in batch order.
    std::size_t collectNewEdges(
        std::span<const Vertex> clique,
        const LocalVertexSet& candidates,
        std::size_t limit,
        std::vector<NewEdge>& edges,
        VertexSearchStats& stats) const;
    // Whether R joined by the batch edge (first, second) at `position`
    // already holds an earlier batch edge.
    [[nodiscard]] bool joinsEarlierEdge(
        std::span<const Vertex> clique,
        const LocalVertexSet& candidates,
        Vertex first,
        Vertex second,
        int position,
        VertexSearchStats& stats) const;
    // The ordinary search below a clique that holds a batch edge, every
    // batch edge before `threshold` excluded (none for a negative one).
    void searchBound(
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        int threshold,
        VertexSearchStats& stats,
        CliqueBatchOutput& output);
    // Moves to X every candidate that `vertex` reaches by a batch edge before
    // `threshold`.
    void excludeEarlierEdges(
        Vertex vertex,
        int threshold,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        VertexSearchStats& stats) const;

    const IncrementalGraph& graph_;
    VertexSearchConfiguration configuration_;
    BatchIndex batchIndex_;
    // One local graph for the enumerator's lifetime, as the edge baseline
    // keeps its own.
    LocalSearchGraph localGraph_;
    std::vector<Vertex> endpoints_;
    std::vector<int> endpointRankByVertex_;
    std::vector<std::size_t> batchDegrees_;
    // Cover seeds only: by endpoint rank, whether the endpoint is in the cover.
    std::vector<unsigned char> inCover_;
    // Batch incidence: endpoint rank -> its batch neighbours, in batch order,
    // with each edge's batch position beside it.
    std::vector<std::size_t> incidenceOffsets_;
    std::vector<std::size_t> incidenceFill_;
    std::vector<Vertex> incidenceNeighbors_;
    std::vector<int> incidencePositions_;
    // The endpoints in the current root's P -- its later endpoints, and with
    // cover seeds the earlier ones outside the cover -- the only vertices of P
    // that can carry a batch edge, as a mask over its graph.
    LocalVertexSet endpointMask_;
    // Marks the current root's P, for the test made before its graph exists.
    std::vector<std::uint32_t> epochByVertex_;
    std::uint32_t epoch_ = 0U;
    // branchUnbound's working sets, one frame per depth of states without a
    // batch edge, reused across roots; a deque, as a deeper state opens while
    // the shallower ones still use theirs.
    struct UnboundFrame {
        LocalVertexSet childCandidates;
        LocalVertexSet childExcluded;
        std::vector<Vertex> extensions;
        std::vector<NewEdge> edges;
    };
    std::deque<UnboundFrame> unboundFrames_;
    std::size_t unboundDepth_ = 0U;
};

} // namespace immce
