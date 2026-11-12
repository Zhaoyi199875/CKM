#pragma once

#include "incremental_graph.hpp"
#include "types.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

namespace immce {

// Third baseline: the procedure-oriented IMCE of
//
//   T. Yu, T. Jiang, M. J. Bah, C. Zhao, H. Huang, M. Liu, S. Zhou, Z. Li,
//   J. Zhang, "Incremental Maximal Clique Enumeration for Hybrid Edge Changes
//   in Large Dynamic Graphs", IEEE TKDE 36(4):1650-1665, 2024,
//
// which keeps every v-header maximal clique (the candidate clique part of its
// index SOMEi) next to an inverted clique tree whose root-to-leaf branches are
// the maximal cliques (the ICTree), and answers a batch by rolling every
// header the batch can touch back to the last snapshot it shares with G + H
// and redoing that header's construction from there.  The construction step
// is the paper's reference [42]: T. Yu and M. Liu, "A linear time algorithm
// for maximal clique enumeration in large sparse graphs", IPL 125:35-40, 2017
// (CM-Constructor, Algorithm 2 UpdateClique).  The paper publishes no code;
// this is a direct implementation of its pseudo-code:
//
//   Algorithm 1  ICTreeDelete             treeDelete
//   Algorithm 2  ICTreeAdd                treeAdd
//   Algorithm 3  UpdateSOMEi              applyBatch
//   Algorithm 4  GetUpdateScope           computeUpdateScope, pruning off
//   Algorithm 5  RollBackTimestamp        rollBack
//   Algorithm 6  UpdateCliqueBodyCluster  redo + updateClique
//   Algorithm 7  GetUpdateScope_Pruning   computeUpdateScope, pruning on
//
// Vertices are ordered by their dense id, which keeps the order of the
// original ids -- the paper's default ("the alphabetical order of vertex
// ID").  A timestamp T_r is the vertex r it follows, so a body's creation
// timestamp BT is stored as that vertex.
//
// Where the pseudo-code is ambiguous, the reading below is the one the
// paper's own lemmas and worked example force:
//
//   * Algorithm 3 line 5 builds the redo scope as Gamma_r(v_j, N^b) u {v_j}:
//     the big neighbours from v_j on.  Line 5 of Algorithm 6 adds {v_r} as a
//     body, so v_r has to be a neighbour of the header; the scope is taken as
//     the big neighbours >= v_j, which for an insertion always holds v_j.
//   * Algorithm 5 line 7 keeps "Gamma_l(v_j, l) u {v_j}": the vertices of l
//     up to the rollback timestamp, which Algorithm 3 sets to T_{j-1}.  A body
//     therefore keeps the vertices smaller than v_j.
//   * Algorithm 7 tests updateNeighborMap[v_j] against N^b_{v_k} only, which
//     on its own would also prune the small endpoint v_x itself (v_x is never
//     in its own N^b).  Lemma 8 prunes only headers below v_x, and the
//     paper's example keeps v_3 for the changes at v_4 and v_5 (BNI = {v_2:
//     v_4, v_3: v_4}), so a header that is itself a changed small endpoint is
//     kept.
//
// The paper handles insertions and deletions by the same procedure.  For an
// insertion the batch's edges are already in `graph` when `applyBatch` runs,
// and the update scope is read off G + H, whose neighbour lists contain G's.
// A deletion needs the union of the two graphs as well -- its small endpoint
// has left N^s_{v_y} in G - H -- so it reads the scope before the edges go and
// redoes the headers after (`beginBatch` / `endBatch`).
//
// Reporting changes.  SOMEi maintains the maximal cliques, it does not list
// how they changed: a header rebuilt from its rollback point deletes its
// cliques from the ICTree and adds the redone ones back, most of them
// unchanged.  To produce the same change stream as the other frameworks,
// every leaf the ICTree gains or loses during a batch is recorded with the
// clique it stands for, and the batch's net change is read off those records
// afterwards: a clique that lost its leaf once more than it gained one is a
// deleted maximal clique, the reverse a new one.  This is extra work the paper
// does not do; it is timed apart (`change-ms`) and `trackChanges = false`
// switches it off, leaving the paper's algorithm alone.
//
// The untimed build (2026-11-06).  A deletion or mixed run first builds its
// starting graph untimed, and the paper's construction of it is CM-Constructor
// itself: hours on youtube, wikipedia and flickr, where a static enumeration
// takes a minute (experiments/somei_build_cost_20261106).  Nearly all of it is
// lines 7-22 of Algorithm 6, which find the maximal cliques of G[R]: the
// maximal lists among the l_j n R, found by comparing R with every body of
// the header and then every such list with every other.  A body that lies in
// R is one of those cliques, and is extended; every other one, K, becomes the
// body K u {v_r}.  So between setBuilding(true) and setBuilding(false) a step
// may instead enumerate the maximal cliques of G[R] by Bron-Kerbosch with
// pivoting and look each one up among the header's bodies: the same bodies,
// created in the same order with the same timestamps, and the same ICTree, so
// the measured batches start from the very state the paper's construction
// leaves.  Those batches always run the paper's steps.
struct SomeiConfiguration {
    // Algorithm 7 (true) or Algorithm 4 (false) for the update scope.
    bool pruning = true;
    // Record the ICTree's leaf changes and net them into the batch's new and
    // deleted maximal cliques.
    bool trackChanges = true;
    // Clock reads around each header's rollback, redo and ICTreeAdd -- the
    // phases of the paper's Fig. 12.  Off, only the batch-level phases are
    // timed.
    bool profilePhases = false;
    // The untimed build's steps by Bron-Kerbosch on G[R] (true) or by the
    // paper's lines 7-22 (false); see above.  The same state either way.
    bool enumerateBuild = true;
    // A build step with R larger than this runs the paper's lines instead:
    // G[R] is held as |R|^2 bits.  The tests lower it.
    std::size_t enumerateBuildLimit = 4096U;
};

struct SomeiStats {
    long long batches = 0;
    long long batchEdges = 0;
    // GetUpdateScope: big endpoints (|SNI|), small neighbours v_k <=
    // SNI[v_j] visited, visits Lemma 8 rejected, and headers to rebuild
    // (|BNI|).
    long long bigVertices = 0;
    long long scopeVisits = 0;
    long long prunedVisits = 0;
    long long rebuiltHeaders = 0;
    // RollBackTimestamp: bodies visited, removed (created after the rollback
    // timestamp) and cut back to it.
    long long rolledBackBodies = 0;
    long long removedBodies = 0;
    long long truncatedBodies = 0;
    // UpdateCliqueBodyCluster: timestamps replayed, those with B or R empty,
    // l_j-against-R comparisons, bodies extended in place and created, and
    // the partial intersections (pSet) left after duplicates.
    long long redoSteps = 0;
    long long emptyRenovations = 0;
    long long bodyComparisons = 0;
    long long extendedBodies = 0;
    long long createdBodies = 0;
    long long partialIntersections = 0;
    // ICTree.
    long long treeAdds = 0;
    long long treeDeletes = 0;
    long long treeNodesCreated = 0;
    long long treeNodesDeleted = 0;
    // Leaf gains and losses recorded, and the net change they gave.
    long long leafEvents = 0;
    long long newCliques = 0;
    long long deletedCliques = 0;
    double scopeMilliseconds = 0.0;
    double rebuildMilliseconds = 0.0;
    double changeMilliseconds = 0.0;
    // profilePhases only: the rebuild split into its three parts.
    double rollbackMilliseconds = 0.0;
    double redoMilliseconds = 0.0;
    double treeAddMilliseconds = 0.0;
};

class SomeiIndex {
public:
    explicit SomeiIndex(
        const IncrementalGraph& graph,
        SomeiConfiguration configuration = {});

    // Algorithm 3 for one batch of insertions.  The batch's edges must already
    // be in the graph.  With change tracking on, `newCliques` and
    // `deletedCliques` receive the batch's new and deleted maximal cliques
    // (sorted vertices, no order between cliques); off, they are left
    // untouched.
    void applyBatch(
        std::span<const EdgeId> batch,
        SomeiStats& stats,
        std::vector<Clique>& newCliques,
        std::vector<Clique>& deletedCliques);
    // Algorithm 3 in its two halves, which is how a batch of deletions runs:
    // the update scope is read off the union of the graphs before and after
    // the batch -- for a deletion that is the graph before it, where a
    // deleted edge's small endpoint is still a small neighbour of its big one
    // -- and the headers are rolled back and redone on the graph after it.
    // So the caller calls `beginBatch` with the batch's edges still in the
    // graph, removes them, then calls `endBatch`.  (For insertions the graph
    // after the batch is the union, and `applyBatch` runs both halves on it.)
    void beginBatch(std::span<const EdgeId> batch, SomeiStats& stats);
    void endBatch(
        SomeiStats& stats,
        std::vector<Clique>& newCliques,
        std::vector<Clique>& deletedCliques);

    [[nodiscard]] const SomeiConfiguration& configuration() const noexcept {
        return configuration_;
    }
    // Change tracking off for an untimed build, whose changes nobody reads,
    // and back on after it.  The ICTree is the same either way.
    void setChangeTracking(bool on) noexcept { configuration_.trackChanges = on; }
    // The batches between setBuilding(true) and setBuilding(false) are the
    // untimed build: their steps enumerate G[R] if `enumerateBuild` is set.
    void setBuilding(bool on) noexcept { building_ = on; }
    // Leaves of the ICTree: |MC(G)|.
    [[nodiscard]] std::size_t maximalCliqueCount() const noexcept {
        return leafCount_;
    }
    // Bodies in the candidate clique part: |HMC(G)|.
    [[nodiscard]] std::size_t candidateCliqueCount() const noexcept {
        return bodyCount_;
    }
    // ICTree nodes, root excluded.
    [[nodiscard]] std::size_t treeNodeCount() const noexcept;

    // Inspection, for the tests.
    struct BodyView {
        std::span<const Vertex> vertices;
        Vertex createdAt = 0;
    };
    [[nodiscard]] std::vector<BodyView> bodies(Vertex header) const;
    // The ICTree's root-to-leaf branches, each as a sorted clique.
    [[nodiscard]] std::vector<Clique> maximalCliques() const;
    // The last batch's BNI: (header, v_j) pairs in header order.
    [[nodiscard]] std::span<const std::pair<Vertex, Vertex>>
    lastUpdateScope() const noexcept {
        return lastScope_;
    }

private:
    using NodeId = std::uint32_t;
    static constexpr NodeId noNode = ~NodeId{0};
    static constexpr NodeId rootNode = 0;

    struct Body {
        std::vector<Vertex> vertices;   // ascending, all above the header
        Vertex createdAt = 0;           // BT: the vertex of the timestamp
        NodeId link = noNode;           // the ICTree node of its header
    };

    struct TreeNode {
        Vertex vertex = -1;             // -1: the root or a free slot
        NodeId parent = noNode;
        std::uint32_t children = 0;
        // Marked nodes only: the linked body, as its index among the bodies
        // of `vertex` (the header).  A linked body is never moved or changed,
        // so the index stays valid while the mark does.
        std::uint32_t body = 0;
        bool marked = false;
    };

    // node.children.find / add / delete of the ICTree: one open-addressing
    // table keyed by (parent, vertex) for every node's children.
    class ChildTable {
    public:
        ChildTable();
        [[nodiscard]] NodeId find(NodeId parent, Vertex vertex) const noexcept;
        void insert(NodeId parent, Vertex vertex, NodeId child);
        void erase(NodeId parent, Vertex vertex);

    private:
        struct Slot {
            std::uint64_t key = emptyKey;
            NodeId child = noNode;
        };
        static constexpr std::uint64_t emptyKey = ~std::uint64_t{0};

        [[nodiscard]] static std::uint64_t key(
            NodeId parent, Vertex vertex) noexcept;
        [[nodiscard]] std::size_t home(std::uint64_t value) const noexcept;
        void grow();

        std::vector<Slot> slots_;
        std::size_t size_ = 0;
        std::size_t mask_ = 0;
    };

    struct LeafEvent {
        std::size_t offset = 0;         // into eventVertices_
        std::uint32_t length = 0;
        int sign = 0;                   // +1 leaf gained, -1 leaf lost
    };

    // Algorithms 4 / 7: fills scope_ with BNI in header order.
    void computeUpdateScope(std::span<const EdgeId> batch, SomeiStats& stats);
    // Lemma 8 (see above): can the changed small neighbours `changed` of a big
    // endpoint alter the cliques of `header`?
    [[nodiscard]] bool touchesHeader(
        Vertex header, std::span<const Vertex> changed) const;
    // Algorithm 5 with timestamp T_{from-1}.
    void rollBack(Vertex header, Vertex from, SomeiStats& stats);
    // Algorithm 6 lines 1-22 over the big neighbours from `from` on.
    void redo(Vertex header, Vertex from, SomeiStats& stats);
    // Algorithm 6 lines 4-22 (UpdateClique of [42]) for one timestamp T_r.
    void updateClique(Vertex header, Vertex renovator, SomeiStats& stats);
    // The same step of the untimed build, lines 7-22 by Bron-Kerbosch on
    // G[R] (see above).
    void updateCliqueByEnumeration(
        Vertex header, Vertex renovator, SomeiStats& stats);
    // Bron-Kerbosch with pivoting on renovationRows_, the P / X rows of
    // `depth` in enumerationStack_; every maximal clique into foundRanges_.
    void enumerateRenovation(std::size_t depth);
    // bodyByHash_ for the header about to be redone.
    void indexBodies(Vertex header);
    void releaseBodyIndex();
    [[nodiscard]] static std::uint64_t extendHash(
        std::uint64_t hash, Vertex vertex) noexcept;
    // Algorithm 6 lines 23-24.
    void addAllToTree(Vertex header, SomeiStats& stats);

    void treeAdd(Vertex header, std::uint32_t bodyIndex, SomeiStats& stats);
    void treeDelete(Vertex header, std::uint32_t bodyIndex, SomeiStats& stats);
    [[nodiscard]] NodeId newChild(NodeId parent, Vertex vertex, SomeiStats& stats);
    void deleteNode(NodeId node, SomeiStats& stats);

    void recordLeaf(NodeId node, int sign);
    void recordLeaf(Vertex header, const Body& body, int sign);
    void netChanges(
        std::vector<Clique>& newCliques,
        std::vector<Clique>& deletedCliques,
        SomeiStats& stats);

    const IncrementalGraph& graph_;
    SomeiConfiguration configuration_;

    // Candidate clique part: the bodies of every header, in the order they
    // were created, so also in creation-timestamp order.
    std::vector<std::vector<Body>> headers_;
    std::size_t bodyCount_ = 0;

    // ICTree.
    std::vector<TreeNode> nodes_;
    std::vector<NodeId> freeNodes_;
    ChildTable children_;
    std::size_t leafCount_ = 0;

    // Per-batch scope: SNI / updateNeighborMap as the batch's edges grouped
    // by big endpoint (changedSmall_ holds one group's small endpoints), BNI
    // as a vertex-indexed array.
    std::vector<std::pair<Vertex, Vertex>> changes_;
    std::vector<Vertex> changedSmall_;
    std::vector<Vertex> rebuildFrom_;
    std::vector<Vertex> scopeHeaders_;
    std::vector<std::pair<Vertex, Vertex>> lastScope_;

    // UpdateClique scratch: R, l_j n R, pSet and fSet as ranges into arenas.
    std::vector<Vertex> renovation_;
    std::vector<Vertex> intersection_;
    std::vector<Vertex> partialVertices_;
    std::vector<std::pair<std::size_t, std::size_t>> partialRanges_;
    std::vector<Vertex> fullVertices_;
    std::vector<std::pair<std::size_t, std::size_t>> fullRanges_;

    // Leaf changes of the current batch.
    std::vector<Vertex> eventVertices_;
    std::vector<LeafEvent> events_;
    std::vector<std::uint32_t> eventOrder_;

    // The untimed build by enumeration: whether the batches are the build,
    // the bodies of the header being redone by the hash of their vertices,
    // G[R] as bit rows over R's positions, the recursion's P and X rows (one
    // pair per depth), its current clique, and the maximal cliques of G[R]
    // (positions, ranges into an arena) and those to create.
    bool building_ = false;
    std::unordered_multimap<std::uint64_t, std::uint32_t> bodyByHash_;
    std::size_t renovationWords_ = 0;
    std::vector<std::uint64_t> renovationRows_;
    std::vector<std::uint64_t> enumerationStack_;
    std::vector<std::uint32_t> enumerationClique_;
    std::vector<std::uint32_t> foundPositions_;
    std::vector<std::pair<std::size_t, std::size_t>> foundRanges_;
    std::vector<std::uint32_t> toCreate_;
};

} // namespace immce
