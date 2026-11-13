#pragma once

#include "batch_index.hpp"
#include "local_bk.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <span>
#include <vector>

namespace ckm {

struct CoverSeedConfiguration {
    // Whether a subtree moves onto a graph of its own when its universe --
    // clique, candidates and excluded set -- fits in fewer 64-bit words than
    // the graph it is on.  Off only in the tests, which hold the search with
    // and without it to be node for node the same.
    bool reindex = true;
};

// One end of a batch edge, tagged with the rank of the seed it is assigned to.
struct SeedIncidence {
    Vertex neighbor = 0;
    int ownerRank = 0;
    // Position of the edge among its seed's edges.
    int order = 0;
};

struct CoverSeedStats {
    long long batches = 0;
    long long batchEdges = 0;
    // The seeds that receive at least one edge: the subproblems.
    long long seeds = 0;
    // The seed cover, and the batch edges with both endpoints in it, which
    // the edge assignment gave to the endpoint with more new edges.
    long long coverVertices = 0;
    long long coverBothEndpoints = 0;
    long long assignedEdges = 0;
    long long maximumSeedEdges = 0;
    // Seeds with one assigned edge, which search from that edge, and the
    // batch edges each kind of root handled: those seeds, the seeds whose
    // root took anchor branching (one tree per partner) and those whose root
    // took pivot branching (one tree).
    long long singleEdgeSeeds = 0;
    long long singleEdgeSeedEdges = 0;
    long long anchorRootEdges = 0;
    long long pivotRootEdges = 0;
    // Subtrees moved onto a graph of their own: the partner trees and the
    // root's partner children (reindexedTrees), the root's other children
    // under pivot branching (rootChildReindexedSubtrees), and every move below
    // those (deepReindexedSubtrees).
    long long reindexedTrees = 0;
    long long rootChildReindexedSubtrees = 0;
    long long deepReindexedSubtrees = 0;
    // Word work: the row width of the graph under every search node, and
    // under every pivot candidate scored.  Two searches with the same node
    // count can differ several times in these when one runs on wider graphs.
    long long nodeWords = 0;
    long long pivotWords = 0;
    // Anchor-or-Pivot branching at the roots of the seeds with several
    // edges: how they chose, the children each branching would have opened
    // there summed over them (|A| and |B|) and the children the chosen one
    // did open; the roots the degree bound decided without a pivot scan, and
    // those that scanned (rootPivotBranches sums |B| over the scanned ones
    // only).
    long long anchorRootSeeds = 0;
    long long pivotRootSeeds = 0;
    long long rootAnchorBranches = 0;
    long long rootPivotBranches = 0;
    long long rootChosenBranches = 0;
    long long rootDegreeDecisions = 0;
    long long rootPivotScans = 0;
    // Unanchored states below a root that chose anchor branching, and of
    // those the states with a single partner left, anchored without a pivot
    // scan.
    long long unanchoredAnchorChoices = 0;
    long long singlePartnerAnchors = 0;
    long long domainVertexOccurrences = 0;
    std::size_t maximumDomainVertices = 0U;
    // The seed cover and the edge assignment, timed once per batch.
    double assignMilliseconds = 0.0;
    long long recursiveCalls = 0;
    long long pivotCalls = 0;
    long long pivotCandidatesEvaluated = 0;
    long long leaves = 0;
    // Earlier edges excluded: vertices moved from P to X, and the incidence
    // entries read to find them.
    long long ownershipMoves = 0;
    long long ownershipProbes = 0;
    // Unanchored states, and those the live test cut because no partner was
    // left in P.
    long long unanchoredNodes = 0;
    long long liveTestPrunes = 0;
    long long emittedCliques = 0;
    LocalGraphStats localGraph;
};

// CoverSeed: the new maximal cliques of an insertion batch.
//
// Seed cover and edge assignment.  Every new edge is assigned to one of its
// endpoints, so a batch splits into seeds s, each with its assigned edges H_s
// and their other endpoints, the partners V_s.  The seeds must form a vertex
// cover of the new edges; they come from the maximal-matching
// 2-approximation of a minimum vertex cover.  An edge with both endpoints in
// the cover goes to the one with more new edges.
//
// Domain.  A seed is searched on one local graph, built over its domain
//     Omega_s = {s} u U_{v in V_s} ({v} u (N(s) n N(v))),
// a subset of the seed's closed neighbourhood.  EdgeSeed (IMCE's search,
// Section II-C) builds one graph per new edge over that edge's domain, and
// the domains of a seed's edges overlap almost completely, so a seed reads
// the shared part once instead of once per partner.  A seed with one
// assigned edge has nothing to share and searches that edge's domain
// exactly as EdgeSeed does.
//
// Earlier edges.  A new maximal clique can lie in several seeds' domains.
// Seeds are ranked by the batch position of their first assigned edge, and a
// clique goes to the first seed whose edge it holds: a vertex joined to the
// clique by an edge of an earlier seed (H_<i) moves to X the moment its other
// endpoint enters the clique, so it is never branched on.
//
// Maximality needs no separate check: a reported clique holds a new edge
// (s, v), and any vertex extending it lies in N(s) n N(v), inside the domain.
//
// Branching.  An unanchored state owes the cliques that hold a partner still
// in P.  Anchor branching opens one child per partner in the anchor set
// A = P n V_s, each later partner yielding to the earlier ones through X;
// pivot branching opens one child per vertex of B = P \ N(q).  Both are exact,
// and Anchor-or-Pivot branching takes anchor branching when |A| <= |B|.  Every
// child of either holds at most |P| - |B| candidates, so the choice keeps the
// worst-case time complexity of BK_pivot.  An anchored state uses pivot
// branching.  Both branchings are implemented in full, and the choice between
// them is made in two places: `searchSeed` at the root and
// `branchUnanchored` below it.
//
// Representation.  A subtree whose universe fits in fewer bitmap words than
// the graph it is on is moved onto a graph of its own, read from the parent's
// bit rows: the partner trees, the root's children under pivot branching, and
// every deeper level.  The search itself is node for node the same.
class CoverSeed {
public:
    CoverSeed(
        const IncrementalGraph& graph,
        CoverSeedConfiguration configuration = {});

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        CoverSeedStats& stats);
    [[nodiscard]] long long countBatch(
        std::span<const EdgeId> batch,
        CoverSeedStats& stats);
    [[nodiscard]] long long visitBatch(
        std::span<const EdgeId> batch,
        CoverSeedStats& stats,
        CliqueVisitor visitor);

private:
    // How an unanchored state branches: by Anchor-or-Pivot, or by the pivot
    // the root already chose.
    enum class BranchPartition {
        AnchorOrPivot,
        Pivot,
    };

    void enumerateBatchInto(
        std::span<const EdgeId> batch,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    void assignEdges(std::span<const EdgeId> batch, CoverSeedStats& stats);
    void buildIncidenceIndex(std::span<const EdgeId> batch);
    void buildDomain(
        Vertex seed,
        std::size_t seedIndex,
        CoverSeedStats& stats);
    // A seed with one assigned edge: EdgeSeed's search on that edge's domain.
    void searchSingleEdgeSeed(
        Vertex seed,
        std::size_t seedIndex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    // A seed with several assigned edges: Anchor-or-Pivot branching at its
    // root.
    void searchSeed(
        Vertex seed,
        std::size_t seedIndex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    // The root state, which both branchings start from: the seed's row,
    // split into rootCandidates_ / rootExcluded_ by the edges of earlier
    // seeds at the seed.
    void prepareRoot(
        Vertex seed,
        int seedRank,
        CoverSeedStats& stats);
    // Anchor branching at the root: one tree per partner, its sets cut from
    // the root state.
    void searchPartnerRoots(
        Vertex seed,
        std::size_t seedIndex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    // An unanchored state: the live test cuts it when no partner is left in
    // P, otherwise it branches as `partition` says.  `atRoot` marks the root
    // state, whose children may move onto a graph of their own.  A
    // `knownPivot` of 0 or more is this state's pivot, already scanned by the
    // caller.  Under Anchor-or-Pivot a state with one partner left anchors it
    // without a pivot scan, and anchor branching below the root reads the
    // partners off the partner mask's bits.
    void branchUnanchored(
        const LocalSearchGraph& graph,
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool atRoot,
        BranchPartition partition,
        Vertex knownPivot,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    // An anchored state, whose clique already holds a partner: pivot
    // branching.  With `mayReindex` a tree on the seed's own graph first
    // moves onto a graph of its own when that is narrower.
    void branchAnchored(
        const LocalSearchGraph& graph,
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool mayReindex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    // Hands a child's subtree to `recurse`, on a graph of its own when that
    // is narrower and on `graph` otherwise.
    template <typename Recurse>
    void descendInto(
        const LocalSearchGraph& graph,
        const Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        CoverSeedStats& stats,
        const Recurse& recurse);
    // Whether a subtree about to be searched on `graph` fits in fewer words
    // on a graph of its own.
    [[nodiscard]] bool reindexPays(
        const LocalSearchGraph& graph,
        const Clique& current,
        const LocalVertexSet& candidates,
        const LocalVertexSet& excluded) const noexcept;
    // The partner mask in `graph`'s own index space.
    [[nodiscard]] const LocalVertexSet& partnerMaskFor(
        const LocalSearchGraph& graph) const noexcept;
    // Moves every candidate that `vertex` reaches by a batch edge ranking
    // below `thresholdRank`/`thresholdOrder` into the excluded set, so it is
    // never branched on.  Entries are visited in ranking order and the walk
    // stops at the first one that no longer qualifies, which is usually the
    // first entry.
    void excludeClaimedNeighbors(
        const LocalSearchGraph& localGraph,
        Vertex vertex,
        int thresholdRank,
        int thresholdOrder,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        CoverSeedStats& stats) const;
    [[nodiscard]] Vertex partnerAt(std::size_t position, Vertex seed) const;
    // True when `edgeId` is a batch edge assigned to an earlier seed.
    [[nodiscard]] bool claimedEarlier(
        EdgeId edgeId,
        int thresholdRank,
        int thresholdOrder) const noexcept;
    [[nodiscard]] std::uint32_t beginEpoch();
    [[nodiscard]] std::uint32_t beginUnionEpoch();

    const IncrementalGraph& graph_;
    CoverSeedConfiguration configuration_;
    BatchIndex batchIndex_;
    // The seed's graph, over its domain (or a single edge's domain).
    LocalSearchGraph localGraph_;
    // The narrow universe of a single tree, re-indexed out of localGraph_.
    LocalSearchGraph treeGraph_;
    std::vector<Vertex> treeVertices_;
    // One graph per nesting level below treeGraph_, each re-indexed out of
    // the one above it, with the partner mask projected into it and the
    // scratch its universe is assembled in.  Grown on demand; deepLevel_ is
    // the next free level.  A level is added while the levels above it are
    // still being searched, and those searches hold references to their graph
    // and mask, so all three must keep their addresses as they grow: hence
    // unique_ptr and deque, not plain vectors.
    std::vector<std::unique_ptr<LocalSearchGraph>> deepGraphs_;
    std::deque<LocalVertexSet> deepMasks_;
    std::deque<std::vector<Vertex>> deepVertices_;
    std::size_t deepLevel_ = 0U;
    // Whether partnerMask_ belongs to the seed being searched, so that a
    // nested graph has to carry it: set while a root that took pivot
    // branching is searched, whose unanchored states read it.
    bool maskActive_ = false;
    // branchUnanchored's working sets, one frame per depth of unanchored
    // states and reused across seeds, so an unanchored state allocates
    // nothing once the frames have grown.  A deeper state opens while the
    // shallower ones still use theirs, so the frames live in a deque, whose
    // growth moves nothing.
    struct UnanchoredFrame {
        LocalVertexSet childCandidates;
        LocalVertexSet childExcluded;
        std::vector<Vertex> extensions;
    };
    std::deque<UnanchoredFrame> unanchoredFrames_;
    std::size_t unanchoredDepth_ = 0U;

    // The batch currently being enumerated, so a seed can recover the
    // partner behind one of its assigned batch positions.
    std::span<const EdgeId> batch_;
    // Seeds, in rank order.  Their assigned edge positions live in a CSR
    // indexed by seed, so a batch needs no per-seed allocation.
    std::vector<Vertex> seeds_;
    std::vector<std::size_t> seedOffsets_;
    std::vector<std::size_t> seedFill_;
    std::vector<std::size_t> ownedCounts_;
    std::vector<int> seedEdgePositions_;
    // Rank of the seed each batch edge is assigned to, by batch position, and
    // the edge's order among that seed's edges.
    std::vector<int> ownerRankByPosition_;
    std::vector<int> orderByPosition_;
    // Seed rank by vertex, -1 for every vertex that is no seed.
    std::vector<int> seedIndexByVertex_;
    // Cover membership by endpoint rank.
    std::vector<unsigned char> inCover_;
    // Batch endpoints in first-appearance order, with their rank and batch
    // degree.
    std::vector<Vertex> endpoints_;
    std::vector<int> endpointRankByVertex_;
    std::vector<std::size_t> batchDegrees_;
    // Batch incidence carrying the rank of the edge's seed next to the
    // neighbour, ordered so the walk over earlier edges can stop early.
    std::vector<std::size_t> incidenceOffsets_;
    std::vector<SeedIncidence> incidenceEntries_;
    std::vector<std::size_t> incidenceFill_;

    std::vector<Vertex> domain_;
    std::vector<Vertex> partners_;
    std::vector<Vertex> candidateVertices_;
    std::vector<Vertex> excludedVertices_;
    std::vector<std::uint32_t> epochByVertex_;
    std::uint32_t epoch_ = 0U;
    std::vector<std::uint32_t> unionEpochByVertex_;
    std::uint32_t unionEpoch_ = 0U;

    LocalVertexSet partnerMask_;
    // Every vertex of the seed's graph: the mask the seed's row is read
    // through.
    LocalVertexSet allVertices_;
    // The root state: the seed's row split into the part no earlier seed
    // reaches by its edges and the part one does.
    LocalVertexSet rootCandidates_;
    LocalVertexSet rootExcluded_;
    // searchPartnerRoots: the partners whose trees have finished, and the
    // current tree's root, all reused across trees.
    LocalVertexSet finishedPartners_;
    LocalVertexSet partnerCandidates_;
    LocalVertexSet partnerExcluded_;
    Clique partnerClique_;
    // A tree's sets once moved onto treeGraph_, and the partner mask in that
    // graph's index space.
    LocalVertexSet compactCandidates_;
    LocalVertexSet compactExcluded_;
    LocalVertexSet treePartnerMask_;
};

} // namespace ckm
