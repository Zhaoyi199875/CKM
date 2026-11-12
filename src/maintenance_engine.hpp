#pragma once

#include "incremental_graph.hpp"
#include "insertion_stream.hpp"
#include "maximal_clique_index.hpp"
#include "new_clique_search.hpp"
#include "somei.hpp"
#include "star_clique_search.hpp"
#include "subsumed_clique_search.hpp"
#include "vertex_clique_search.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <ostream>
#include <span>
#include <vector>

namespace immce {

// A batch that deletes some edges and inserts others is one change: the
// graph before it becomes the graph after it, and the maximal cliques change
// by M(after) \ M(before) and M(before) \ M(after).  Two ways to get there
// (designs/mixed_semantics_plan_20261010.md):
enum class MixedStrategy {
    // Over the two real graphs.  The cliques the batch destroys are found on
    // the graph before it, the new cliques holding an inserted edge on the
    // graph after it, and the candidates the deletions free are checked for
    // maximality on the graph after it -- so a clique that would exist only
    // between the deletions and the insertions is never formed.
    SinglePass,
    // Das et al., Sect. 4.5: every deletion first, then every insertion, each
    // maintained as a batch of its own through the graph between them; the
    // cliques the deletions create and the insertions subsume ("transient"
    // cliques) are then taken out of both change lists.
    TwoPhase,
};

[[nodiscard]] const char* mixedStrategyName(MixedStrategy strategy) noexcept;

struct MaintenanceConfiguration {
    bool runStar = false;
    bool runEdge = false;
    bool runVertex = false;
    bool runSomei = false;
    bool fullMaintenance = false;
    // Whether `build` is used: it inserts the state the measured batches
    // start from, untimed, with star and residual-kernel whatever the
    // measured framework is (any exact pair leaves the same maximal cliques;
    // SOMEi builds its own structure).
    bool builds = false;
    StarBranchScope starBranchScope = StarBranchScope::Unbound;
    // Anything but the rule is the ablation; see StarBranching.
    StarBranching starBranching = StarBranching::Rule;
    // Off only for the ablation; see StarSearchConfiguration::rootFilter.
    bool starRootFilter = true;
    // The vertex baseline's branching; anything but the pivot is the
    // ablation (see VertexBranching).
    VertexBranching vertexBranching = VertexBranching::Pivot;
    // Which endpoints seed it; only the cover for the seed-cover ablation
    // (see VertexSeeds).
    VertexSeeds vertexSeeds = VertexSeeds::Endpoints;
    SubsumedCliqueAlgorithm subsumedAlgorithm =
        SubsumedCliqueAlgorithm::Optimized;
    // The index's live clique-size histogram: the residual-kernel size gate.
    bool sizeHistogram = false;
    SomeiConfiguration somei;
    MixedStrategy mixedStrategy = MixedStrategy::SinglePass;
};

// Milliseconds over the measured batches.
struct MaintenanceTimes {
    double edge = 0.0;
    double vertex = 0.0;
    double star = 0.0;
    // The subsumed search, and for deletions its run backwards (recovery).
    double subsumed = 0.0;
    // SOMEi: the whole of Algorithm 3 per batch, change listing included.
    double somei = 0.0;
    double cliqueIndexUpdate = 0.0;
    double graphUpdate = 0.0;
    // Streaming only: the change lines written, and the whole batch from
    // `prepareBatch` through the closing `finishBatch`.
    double streamOutput = 0.0;
    double streamOuter = 0.0;
    // Splits of the above for batches that delete: the producers on the
    // deleted edges and on the inserted ones, the recovery's share of
    // `subsumed`, and (two-phase) taking the transient cliques out.
    double destroyedSearch = 0.0;
    double newSearch = 0.0;
    double recovery = 0.0;
    double transientNetting = 0.0;
    // The single-vertex cliques: reading the batch endpoints' degrees.
    double singleVertex = 0.0;
};

// The parts of the measured batches' changes, by where they come from.
struct UpdateCounts {
    long long batches = 0;
    long long insertedEdges = 0;
    long long deletedEdges = 0;
    // Old cliques holding a deleted edge.
    long long destroyed = 0;
    // Cliques the deletions create that the batch leaves maximal.
    long long recovered = 0;
    // New cliques holding an inserted edge.
    long long holdingNew = 0;
    // Old cliques the new ones subsume.
    long long subsumed = 0;
    // Two-phase: cliques the deletions created and the insertions subsumed.
    long long transient = 0;
    // The producers' recursive calls on the deleted and the inserted edges.
    long long destroyedSideCalls = 0;
    long long newSideCalls = 0;
    // Single-vertex cliques the batches created and removed, and (two-phase)
    // those the deletions created and the insertions removed again.  The
    // parts above count cliques of two or more vertices only.
    long long singleVertexCreated = 0;
    long long singleVertexRemoved = 0;
    long long transientSingleVertex = 0;
};

// One batch's changes as found: canonical cliques, in no particular order.
struct BatchChanges {
    std::vector<Clique> created;
    std::vector<Clique> removed;
};

// Where the streaming pipeline writes its change lines as it finds them.
struct StreamTarget {
    std::ostream& output;
    const InsertionStream& stream;
    bool writeNew = false;
    bool writeSubsumed = false;
};

// Keeps the maximal cliques of an IncrementalGraph across its batches: the
// new-clique searches the configuration selects (compared batch by batch when
// there are several), the subsumed search and its index, or SOMEi in their
// place.  Every graph update goes through it.
//
// Single-vertex cliques.  The graph's vertices are those that have appeared
// in the stream: a vertex appears with its first edge, as an isolated vertex
// of the graph before that batch, and stays when it loses its edges.  An
// isolated vertex v is the maximal clique {v}, as in the baselines (Das et
// al. count {a}, {b}, {c} as subsumed when three isolated vertices become a
// triangle).  Such a clique changes only when a batch endpoint's degree
// crosses zero: {v} is subsumed when v had no edge before the batch and has
// one after it, and new in the reverse case.  The graph's updates find exactly
// those vertices anyway -- the ones an insertion gives a first edge and the
// ones a deletion leaves without -- to keep their list of active vertices, so
// the engine reads the cliques off those two lists, the same way for every
// framework, and they never enter the searches or the clique index, which
// hold the cliques of two or more vertices.  Only full maintenance lists
// them, and SOMEi only with its change tracking on; `currentCliqueCount`
// counts them whenever it counts cliques.
class MaintenanceEngine {
public:
    MaintenanceEngine(
        IncrementalGraph& graph,
        const MaintenanceConfiguration& configuration);

    // A warm-up batch that is not searched: the graph takes it, nothing else.
    void applyOnly(std::span<const EdgeId> edges);
    // Untimed: inserts edges into the state the measured batches start from.
    void build(std::span<const EdgeId> edges);
    // An insertion, its new cliques counted: one producer, no maintenance,
    // no output.
    void countBatch(std::span<const EdgeId> batch);
    // An insertion streamed through full maintenance, its change lines
    // written as they are found.
    void streamBatch(const EdgeBatch& batch, const StreamTarget& target);
    // A materialized batch: insertions (new cliques only, or with full
    // maintenance), deletions, or both at once (full maintenance only).
    void updateBatch(
        int timestamp,
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    // The end of the warm-up: every timer and counter starts again.
    void resetMeasurements();

    [[nodiscard]] const MaintenanceTimes& times() const noexcept {
        return times_;
    }
    [[nodiscard]] const UpdateCounts& updateCounts() const noexcept {
        return counts_;
    }
    [[nodiscard]] long long totalCreated() const noexcept { return totalCreated_; }
    [[nodiscard]] long long totalRemoved() const noexcept { return totalRemoved_; }
    [[nodiscard]] const NewCliqueSearchStats& edgeStats() const noexcept {
        return edgeStats_;
    }
    [[nodiscard]] const VertexSearchStats& vertexStats() const noexcept {
        return vertexStats_;
    }
    [[nodiscard]] const StarSearchStats& starStats() const noexcept {
        return starStats_;
    }
    [[nodiscard]] const SubsumedCliqueSearchStats& subsumedStats() const noexcept {
        return subsumedStats_;
    }
    [[nodiscard]] const SomeiStats& someiStats() const noexcept {
        return someiStats_;
    }
    [[nodiscard]] const MaximalCliqueIndex& index() const noexcept {
        return index_;
    }
    [[nodiscard]] const SomeiIndex* somei() const noexcept {
        return someiIndex_.get();
    }
    // The maximal cliques held now, the single-vertex ones included (none
    // when only new cliques are searched).
    [[nodiscard]] std::size_t currentCliqueCount() const noexcept;
    // The vertices that have appeared and have no edge now: the single-vertex
    // maximal cliques.
    [[nodiscard]] std::size_t isolatedVertexCount() const noexcept;

private:
    // Every selected producer on `batch`, compared when there are several.
    [[nodiscard]] std::vector<Clique> searchCliques(
        std::span<const EdgeId> batch,
        int timestamp);
    // The same, timed into one side of a batch that deletes.
    [[nodiscard]] std::vector<Clique> searchSide(
        std::span<const EdgeId> batch,
        int timestamp,
        bool deletedSide);
    [[nodiscard]] long long producerCalls() const noexcept;
    [[nodiscard]] std::vector<Clique> searchSubsumed(
        std::span<const EdgeId> inserted,
        const std::vector<Clique>& fresh);
    [[nodiscard]] std::vector<Clique> recover(
        std::span<const EdgeId> deleted,
        const std::vector<Clique>& destroyed);
    void insertFresh(const std::vector<Clique>& fresh);
    void eraseDestroyed(const std::vector<Clique>& destroyed);
    void updateSinglePass(
        int timestamp,
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    void updateTwoPhase(
        int timestamp,
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    void updateSomei(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    // Whether the batches' single-vertex cliques are listed: full
    // maintenance, and for SOMEi its change tracking.
    [[nodiscard]] bool listsSingleVertexCliques() const noexcept;
    // The graph updates of a maintained batch, timed as graph updates.  Each
    // notes that the graph's list of the vertices it activated (or
    // isolated) belongs to this batch.
    void insertIntoGraph(std::span<const EdgeId> edges);
    void removeFromGraph(std::span<const EdgeId> edges);
    // Marks the vertices `activated` as appeared.
    void markAppeared(std::span<const Vertex> activated);
    // The batch's single-vertex cliques from the vertices its insertions
    // activated (A) and its deletions isolated (R), both ascending: {v} is
    // subsumed for v in A \ R and new for v in R \ A, whichever order the
    // two updates ran in -- removing first, A also holds the vertices R
    // regained, whose cliques (two-phase: transient) did not change; inserting
    // first, A and R are disjoint.  Fills singleVertexRemoved_ and
    // singleVertexCreated_, counts them, and marks A as appeared.
    void settleSingleVertexCliques(bool twoPhase);

    IncrementalGraph& graph_;
    MaintenanceConfiguration configuration_;
    std::unique_ptr<NewCliqueEnumerator> edgeEnumerator_;
    std::unique_ptr<VertexCliqueEnumerator> vertexEnumerator_;
    std::unique_ptr<StarCliqueEnumerator> starEnumerator_;
    MaximalCliqueIndex index_;
    std::unique_ptr<SomeiIndex> someiIndex_;
    std::unique_ptr<SubsumedCliqueEnumerator> subsumedEnumerator_;
    std::unique_ptr<StarCliqueEnumerator> buildEnumerator_;

    NewCliqueSearchStats edgeStats_;
    VertexSearchStats vertexStats_;
    StarSearchStats starStats_;
    SubsumedCliqueSearchStats subsumedStats_;
    SomeiStats someiStats_;
    StarSearchStats buildStarStats_;
    SubsumedCliqueSearchStats buildSubsumedStats_;
    SomeiStats buildSomeiStats_;
    MaintenanceTimes times_;
    UpdateCounts counts_;
    long long totalCreated_ = 0;
    long long totalRemoved_ = 0;
    // A SOMEi batch that deletes and inserts: its changed edges together.
    std::vector<EdgeId> changedEdges_;
    // Single-vertex cliques: which vertices have appeared, and how many;
    // whether this batch inserted and deleted edges, so that the graph's
    // lists of activated and isolated vertices are this batch's; and the
    // batch's subsumed and new single-vertex cliques.
    std::vector<unsigned char> appeared_;
    std::size_t appearedCount_ = 0U;
    bool insertedThisBatch_ = false;
    bool removedThisBatch_ = false;
    std::vector<Vertex> singleVertexRemoved_;
    std::vector<Vertex> singleVertexCreated_;
};

} // namespace immce
