#pragma once

#include "batch_index.hpp"
#include "local_bk.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <span>
#include <vector>

namespace immce {

// Where the branch-count rule chooses between the two partitions of a
// multi-edge star's search (see StarCliqueEnumerator). The two settings are
// the two star methods.
enum class StarBranchScope {
    // At the centre only. A star whose root takes the pivot runs the pivot
    // partition all the way down; one whose root takes the new edges runs one
    // tree per partner.
    Root,
    // At every state whose clique holds the centre and no partner yet: the
    // centre, and below it every child the pivot partition opens on a vertex
    // that is not a partner. Each such state is still a search for the cliques
    // holding at least one of the partners left in P, so either partition is
    // exact there too. The default.
    Unbound,
};

[[nodiscard]] const char* starBranchScopeName(StarBranchScope scope) noexcept;

// Which partition a multi-edge star takes at a state whose clique holds the
// centre and no partner yet. The rule is the method; the two fixed settings
// are its ablation (2026-10-31), and the scope only matters under the rule.
enum class StarBranching {
    // The partition with fewer children, where StarBranchScope says. The
    // default.
    Rule,
    // Always the new edges: one child per partner left, as separate searches
    // of the star's edges would open. The centre takes it, so every state
    // below holds a partner and the rule is never asked.
    Anchor,
    // Always the pivot, at the centre and at every state below it that holds
    // no partner yet, as the vertex baseline branches.
    Pivot,
};

[[nodiscard]] const char* starBranchingName(StarBranching branching) noexcept;

struct StarSearchConfiguration {
    StarBranchScope branchScope = StarBranchScope::Unbound;
    StarBranching branching = StarBranching::Rule;
    // The root filter: a multi-edge star's graph is built over its domain, the
    // centre, its partners and the centre's neighbours adjacent to a partner,
    // so its root's P and X hold no vertex that no new edge of the star
    // reaches. Off is the ablation: the graph covers the centre's whole
    // closed neighbourhood and its root's P and X keep every neighbour.
    // Singleton stars are not affected: their root already holds the edge.
    bool rootFilter = true;
    // The re-indexing gate. A subtree moves onto a graph of its own when its
    // universe -- clique, candidates and excluded set -- fits in fewer 64-bit
    // words than the graph it is on, and it still has at least this many
    // candidates. The word test asks whether the move narrows the rows at
    // all; keying it on domain width instead re-indexes trees that are
    // already one word wide. The candidate count predicts whether the tree is
    // then searched long enough to repay the move: trees that average two
    // nodes never do, however wide their domains. The tests lower it so that
    // their small graphs move subtrees too.
    std::size_t reindexMinCandidates = 8U;
};

// One end of a batch edge, tagged with the rank of the star that owns it.
struct StarIncidence {
    Vertex neighbor = 0;
    int ownerRank = 0;
    // Position of the edge inside its own star, for tie-breaks within a star.
    int order = 0;
};

struct StarSearchStats {
    long long batches = 0;
    long long batchEdges = 0;
    // One star per centre that owns at least one edge.
    long long stars = 0;
    // The vertex cover the centres come from, and the batch edges with both
    // endpoints in it, which the tie-break placed rather than the cover.
    long long coverVertices = 0;
    long long coverBothEndpoints = 0;
    long long starEdgeOccurrences = 0;
    long long maximumStarEdges = 0;
    // Stars owning one edge, which search that edge's support alone, and the
    // batch edges each kind of root handled: the singletons, the multi-edge
    // stars whose root took the new edges (one tree per partner) and those
    // whose root took the pivot (one merged tree).
    long long singletonStars = 0;
    long long singletonStarEdges = 0;
    long long perPartnerStarEdges = 0;
    long long mergedStarEdges = 0;
    // Subtrees moved onto a graph of their own: per-partner trees and the
    // centre's partner children (reindexedTrees), the pivot partition's
    // first-level subtrees on non-partners (mergedReindexedSubtrees), and every
    // move below those (deepReindexedSubtrees).
    long long reindexedTrees = 0;
    long long mergedReindexedSubtrees = 0;
    long long deepReindexedSubtrees = 0;
    // Word work: the row width of the graph under every search node, and
    // under every pivot candidate scored. Two searches with the same node
    // count can differ several times in these when one runs on wider graphs.
    long long nodeWords = 0;
    long long pivotWords = 0;
    // The branch-count rule at the roots of the multi-edge stars: how they
    // chose, the children each partition would have opened there summed over
    // them (|W| and |B|) and the children the chosen one did open; the roots
    // the degree bound decided without a pivot scan, and those that scanned
    // (branchRootPivotBranches sums |B| over the scanned ones only).
    long long branchRootEdgeStars = 0;
    long long branchRootPivotStars = 0;
    long long branchRootEdgeBranches = 0;
    long long branchRootPivotBranches = 0;
    long long branchRootChosenBranches = 0;
    long long branchRootDegreeDecisions = 0;
    long long branchRootScans = 0;
    // Unbound scope: states below a root that switched to the new edges, and
    // of those the states with a single partner left, bound without a scan.
    long long branchUnboundBindings = 0;
    long long branchUnitBindings = 0;
    long long domainVertexOccurrences = 0;
    std::size_t maximumDomainVertices = 0U;
    // Adjacency entries read, and union marks written, assembling the domains.
    long long domainScanReads = 0;
    long long domainMarkWrites = 0;
    // The partition, timed once per batch.
    double assignMilliseconds = 0.0;
    long long recursiveCalls = 0;
    long long pivotCalls = 0;
    long long pivotCandidatesEvaluated = 0;
    long long leaves = 0;
    long long ownershipMoves = 0;
    long long ownershipProbes = 0;
    long long ownershipWalks = 0;
    // States whose clique holds no partner yet, and those cut because no
    // partner was left in P.
    long long partnerlessNodes = 0;
    long long partnerPrunes = 0;
    long long emittedCliques = 0;
    LocalGraphStats localGraph;
};

// Star-seeded enumeration: the method of the new-clique search.
//
// Partition. Every batch edge goes to one of its endpoints, so a batch splits
// into stars: a centre u together with the partners P_u it was given. A centre
// set has to be a vertex cover of the batch edges, so the fewest stars a batch
// can have is its minimum vertex cover; the centres are taken from the
// maximal-matching 2-approximation of it (the Bar-Yehuda--Even local ratio
// scheme at unit weights).
//
// Domain. A star is searched on one local graph, built over its domain
//     D_u = {u} u U_{v in P_u} ({v} u (N(u) n N(v))),
// a subset of u's closed neighbourhood. The edge-seeded baseline builds one
// graph per batch edge at the cost of the degrees summed over that edge's
// support, and a centre's supports overlap almost completely, so a star reads
// the shared part once instead of once per partner. A star owning one edge has
// nothing to share and searches that edge's support exactly as the baseline
// does.
//
// Ownership. A new maximal clique can lie in several stars' domains. Stars are
// ranked by the batch position of their first edge and the edges of a star by
// batch order, and a clique goes to the first (star, edge) it holds. That is
// enforced the way the baseline enforces its own order: a vertex joined to the
// clique by an earlier batch edge moves to the excluded set the moment its
// counterpart enters the clique, so it is never branched on, and no leaf is
// re-searched to find out it was already reported.
//
// Maximality needs no separate check: a qualifying clique holds a batch edge
// (u, v), and any vertex extending it lies in N(u) n N(v), inside the domain.
//
// Search. Below the centre of a multi-edge star two partitions of the cliques
// the star owes are exact: one child per partner still in P -- the new edges,
// each later partner yielding to the earlier ones through X -- or one child
// per vertex of the pivot's branch set B = P \ N(q). The branch-count rule
// takes the new edges when they open no more children, |W| <= |B|, and the
// pivot otherwise; StarBranchScope says at which states it chooses. Every
// child of either partition holds at most |P| - |B| candidates, so the rule
// keeps the Tomita bound and picks the partition with fewer children under it.
//
// Representation. A subtree whose universe fits in fewer bitmap words than the
// graph it is on is moved onto a graph of its own, read from the parent's bit
// rows: the per-partner trees, the first-level subtrees of the pivot
// partition, and every deeper level, under the gate of
// StarSearchConfiguration. The search itself is node for node the same.
class StarCliqueEnumerator {
public:
    StarCliqueEnumerator(
        const IncrementalGraph& graph,
        StarSearchConfiguration configuration = {});

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        StarSearchStats& stats);
    [[nodiscard]] long long countBatch(
        std::span<const EdgeId> batch,
        StarSearchStats& stats);
    [[nodiscard]] long long visitBatch(
        std::span<const EdgeId> batch,
        StarSearchStats& stats,
        CliqueVisitor visitor);

private:
    // Which partition an unbound state expands: the rule's choice, or the
    // pivot's.
    enum class BranchPartition {
        Rule,
        Pivot,
    };

    void enumerateBatchInto(
        std::span<const EdgeId> batch,
        StarSearchStats& stats,
        CliqueBatchOutput& output);
    void assignStars(std::span<const EdgeId> batch, StarSearchStats& stats);
    void buildIncidenceIndex(std::span<const EdgeId> batch);
    void buildDomain(
        Vertex hub,
        std::size_t starIndex,
        StarSearchStats& stats);
    // A star owning exactly one edge: the edge-seeded search on its support.
    void searchSingleton(
        Vertex hub,
        std::size_t starIndex,
        int starRank,
        StarSearchStats& stats,
        CliqueBatchOutput& output);
    // A star owning several edges: the branch-count rule at its centre.
    void searchStar(
        Vertex hub,
        std::size_t starIndex,
        int starRank,
        StarSearchStats& stats,
        CliqueBatchOutput& output);
    // The centre's state, which both partitions start from: its row, split
    // into rootCandidates_ / rootExcluded_ by what earlier stars claim from it.
    void prepareRoot(
        Vertex hub,
        int starRank,
        StarSearchStats& stats);
    // The new-edge partition at the centre: one tree per partner, its sets cut
    // from the centre's state.
    void searchPartnerRoots(
        Vertex hub,
        std::size_t starIndex,
        int starRank,
        StarSearchStats& stats,
        CliqueBatchOutput& output);
    // A state whose clique holds the centre and no partner yet: prunes it
    // when no partner is left in P, otherwise expands `partition`. `atRoot`
    // marks the centre's own state, whose children may move onto a graph of
    // their own. A `knownPivot` of 0 or more is this state's pivot, already
    // scanned by the caller. Under the rule a state with one partner left
    // binds it without a pivot scan, and the new edges below the centre are
    // read off the partner mask's bits.
    void branchUnbound(
        const LocalSearchGraph& graph,
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool atRoot,
        BranchPartition partition,
        Vertex knownPivot,
        int starRank,
        StarSearchStats& stats,
        CliqueBatchOutput& output);
    // A state whose clique already holds a partner: the ordinary search. With
    // `mayReindex` a tree on the star's own graph first moves onto a graph of
    // its own when the gate opens.
    void branchBound(
        const LocalSearchGraph& graph,
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool mayReindex,
        int starRank,
        StarSearchStats& stats,
        CliqueBatchOutput& output);
    // Hands a child's subtree to `recurse`, on a graph of its own when the
    // gate opens and on `graph` otherwise.
    template <typename Recurse>
    void descendInto(
        const LocalSearchGraph& graph,
        const Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        StarSearchStats& stats,
        const Recurse& recurse);
    // The re-indexing gate for a subtree about to be searched on `graph`.
    [[nodiscard]] bool reindexPays(
        const LocalSearchGraph& graph,
        const Clique& current,
        const LocalVertexSet& candidates,
        const LocalVertexSet& excluded) const noexcept;
    // The partner mask in `graph`'s own index space.
    [[nodiscard]] const LocalVertexSet& partnerMaskFor(
        const LocalSearchGraph& graph) const noexcept;
    // Moves every candidate that `vertex` reaches by a batch edge ranking
    // below `thresholdRank`/`thresholdOrder` into the excluded set, so it
    // is never branched on. Entries are visited in ranking order and the walk
    // stops at the first one that no longer qualifies, which is usually the
    // first entry.
    void excludeClaimedNeighbors(
        const LocalSearchGraph& localGraph,
        Vertex vertex,
        int thresholdRank,
        int thresholdOrder,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        StarSearchStats& stats) const;
    [[nodiscard]] Vertex partnerAt(std::size_t position, Vertex hub) const;
    // True when `edgeId` is a batch edge that a lower-ranked owner claims.
    [[nodiscard]] bool claimedEarlier(
        EdgeId edgeId,
        int thresholdRank,
        int thresholdOrder) const noexcept;
    [[nodiscard]] std::uint32_t beginEpoch();
    [[nodiscard]] std::uint32_t beginUnionEpoch();

    const IncrementalGraph& graph_;
    StarSearchConfiguration configuration_;
    BatchIndex batchIndex_;
    // The star's graph, over its domain (or a singleton's support).
    LocalSearchGraph localGraph_;
    // The narrow universe of a single tree, re-indexed out of localGraph_.
    LocalSearchGraph treeGraph_;
    std::vector<Vertex> treeVertices_;
    // One graph per nesting level below treeGraph_, each re-indexed out of
    // the one above it, with the partner mask projected into it and the
    // scratch its universe is assembled in. Grown on demand; deepLevel_ is
    // the next free level. A level is added while the levels above it are
    // still being searched, and those searches hold references to their graph
    // and mask, so all three must keep their addresses as they grow: hence
    // unique_ptr and deque, not plain vectors.
    std::vector<std::unique_ptr<LocalSearchGraph>> deepGraphs_;
    std::deque<LocalVertexSet> deepMasks_;
    std::deque<std::vector<Vertex>> deepVertices_;
    std::size_t deepLevel_ = 0U;
    // Whether partnerMask_ belongs to the star being searched, so that a
    // nested graph has to carry it: set while the pivot partition runs, whose
    // unbound states read it.
    bool maskActive_ = false;
    // branchUnbound's working sets, one frame per depth of unbound states and
    // reused across stars, so an unbound state allocates nothing once the
    // frames have grown. A deeper state opens while the shallower ones still
    // use theirs, so the frames live in a deque, whose growth moves nothing.
    struct UnboundFrame {
        LocalVertexSet childCandidates;
        LocalVertexSet childExcluded;
        std::vector<Vertex> extensions;
    };
    std::deque<UnboundFrame> unboundFrames_;
    std::size_t unboundDepth_ = 0U;

    // The batch currently being enumerated, so a star can recover the
    // partner behind one of its owned batch positions.
    std::span<const EdgeId> batch_;
    // Stars, in rank order. Owned edge positions live in a CSR indexed by
    // star, so a batch needs no per-star allocation.
    std::vector<Vertex> starHubs_;
    std::vector<std::size_t> starOffsets_;
    std::vector<std::size_t> starFill_;
    std::vector<std::size_t> ownedCounts_;
    std::vector<int> starEdgePositions_;
    // Star rank of the owner of each batch edge, by batch position, and the
    // edge's order inside that star.
    std::vector<int> ownerRankByPosition_;
    std::vector<int> orderByPosition_;
    // Star rank by centre vertex, -1 for every other vertex.
    std::vector<int> starIndexByVertex_;
    // Cover membership by endpoint rank.
    std::vector<unsigned char> inCover_;
    // Batch endpoints in first-appearance order, with their rank and batch
    // degree.
    std::vector<Vertex> endpoints_;
    std::vector<int> endpointRankByVertex_;
    std::vector<std::size_t> batchDegrees_;
    // Batch incidence carrying the owning star's rank next to the neighbour,
    // ordered so the ownership walk can stop early.
    std::vector<std::size_t> incidenceOffsets_;
    std::vector<StarIncidence> incidenceEntries_;
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
    // Every vertex of the star's graph: the mask the centre's row is read
    // through.
    LocalVertexSet allVertices_;
    // The centre's state: its row split into the part no earlier star claims
    // and the part one does.
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

} // namespace immce
