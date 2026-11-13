#pragma once

#include "cover_seed.hpp"
#include "incremental_graph.hpp"
#include "maximal_clique_index.hpp"
#include "subsumed_search.hpp"
#include "vertex_seed.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace ckm {

// The search that finds the new maximal cliques of an insertion batch.
enum class NewCliqueSearch {
    CoverSeed,
    VertexSeed,
};

[[nodiscard]] const char* newCliqueSearchName(NewCliqueSearch search) noexcept;

struct MaintenanceConfiguration {
    NewCliqueSearch newCliqueSearch = NewCliqueSearch::CoverSeed;
    bool fullMaintenance = false;
    // Whether `build` is used: it inserts the state the measured batches
    // start from, untimed, with CoverSeed and KernelEnum whichever new-clique
    // search is measured (both leave the same maximal cliques).
    bool builds = false;
    bool sizeHistogram = false;
};

// Milliseconds over the measured batches.
struct MaintenanceTimes {
    double vertexSeed = 0.0;
    double coverSeed = 0.0;
    // The subsumed search, and for deletions its run backwards (recovery).
    double subsumed = 0.0;
    double cliqueIndexUpdate = 0.0;
    double graphUpdate = 0.0;
    // Splits of the above for batches that delete: the new-clique search on
    // the deleted edges and on the inserted ones, the recovery's share of
    // `subsumed`, and (hybrid batches) taking the transient cliques out.
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
    // Hybrid batches: cliques the deletions created and the insertions
    // subsumed.
    long long transient = 0;
    // The new-clique search's recursive calls on the deleted and the inserted
    // edges.
    long long destroyedSideCalls = 0;
    long long newSideCalls = 0;
    // Single-vertex cliques the batches created and removed, and (hybrid
    // batches) those the deletions created and the insertions removed again.
    // The parts above count cliques of two or more vertices only.
    long long singleVertexCreated = 0;
    long long singleVertexRemoved = 0;
    long long transientSingleVertex = 0;
};

// One batch's changes as found: canonical cliques, in no particular order.
struct BatchChanges {
    std::vector<Clique> created;
    std::vector<Clique> removed;
};

// Keeps the maximal cliques of an IncrementalGraph across its batches: the
// new-clique search the configuration selects, then KernelEnum and the
// maximal-clique index.  Every graph update goes through it.
//
// Single-vertex cliques.  The graph's vertices are those that have appeared
// in the stream: a vertex appears with its first edge, as an isolated vertex
// of the graph before that batch, and stays when it loses its edges.  An
// isolated vertex v is the maximal clique {v}, as in Das et al., who count
// {a}, {b}, {c} as subsumed when three isolated vertices become a triangle.
// Such a clique changes only when a batch endpoint's degree crosses zero:
// {v} is subsumed when v had no edge before the batch and has one after it,
// and new in the reverse case.  The graph's updates find exactly those
// vertices anyway -- the ones an insertion gives a first edge and the ones a
// deletion leaves without -- to keep their list of active vertices, so the
// engine reads the cliques off those two lists, and they never enter the
// searches or the clique index, which hold the cliques of two or more
// vertices.  Only full maintenance lists them; `currentCliqueCount` counts
// them whenever it counts cliques.
class MaintenanceEngine {
public:
    MaintenanceEngine(
        IncrementalGraph& graph,
        const MaintenanceConfiguration& configuration);

    // Untimed: inserts edges into the state the measured batches start from.
    void build(std::span<const EdgeId> edges);
    // An insertion, its new cliques counted: no maintenance, no output.
    void countBatch(std::span<const EdgeId> batch);
    // A materialized batch: insertions (new cliques only, or with full
    // maintenance), deletions, or both at once (full maintenance only).
    void updateBatch(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    // After the build: every timer and counter starts again.
    void resetMeasurements();

    [[nodiscard]] const MaintenanceTimes& times() const noexcept {
        return times_;
    }
    [[nodiscard]] const UpdateCounts& updateCounts() const noexcept {
        return counts_;
    }
    [[nodiscard]] long long totalCreated() const noexcept { return totalCreated_; }
    [[nodiscard]] long long totalRemoved() const noexcept { return totalRemoved_; }
    [[nodiscard]] const VertexSeedStats& vertexSeedStats() const noexcept {
        return vertexSeedStats_;
    }
    [[nodiscard]] const CoverSeedStats& coverSeedStats() const noexcept {
        return coverSeedStats_;
    }
    [[nodiscard]] const SubsumedSearchStats& subsumedStats() const noexcept {
        return subsumedStats_;
    }
    [[nodiscard]] const MaximalCliqueIndex& index() const noexcept {
        return index_;
    }
    // The maximal cliques held now, the single-vertex ones included (none
    // when only new cliques are searched).
    [[nodiscard]] std::size_t currentCliqueCount() const noexcept;
    // The vertices that have appeared and have no edge now: the single-vertex
    // maximal cliques.
    [[nodiscard]] std::size_t isolatedVertexCount() const noexcept;

private:
    // The selected new-clique search on `batch`: canonical cliques, sorted.
    [[nodiscard]] std::vector<Clique> searchCliques(std::span<const EdgeId> batch);
    // The same, timed into one side of a batch that deletes.
    [[nodiscard]] std::vector<Clique> searchSide(
        std::span<const EdgeId> batch,
        bool deletedSide);
    [[nodiscard]] long long searchCalls() const noexcept;
    [[nodiscard]] std::vector<Clique> searchSubsumed(
        std::span<const EdgeId> inserted,
        const std::vector<Clique>& fresh);
    [[nodiscard]] std::vector<Clique> recover(
        std::span<const EdgeId> deleted,
        const std::vector<Clique>& destroyed);
    void insertFresh(const std::vector<Clique>& fresh);
    void eraseDestroyed(const std::vector<Clique>& destroyed);
    // A batch that only inserts or only deletes.
    void updateOneWay(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    // A hybrid batch: its deletions first, then its insertions.
    void updateHybrid(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
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
    // regained, whose cliques (transient in a hybrid batch) did not change;
    // inserting first, A and R are disjoint.  Fills singleVertexRemoved_ and
    // singleVertexCreated_, counts them, and marks A as appeared.
    void settleSingleVertexCliques(bool hybrid);

    IncrementalGraph& graph_;
    MaintenanceConfiguration configuration_;
    std::unique_ptr<VertexSeed> vertexSeed_;
    std::unique_ptr<CoverSeed> coverSeed_;
    MaximalCliqueIndex index_;
    std::unique_ptr<SubsumedSearch> subsumedSearch_;
    std::unique_ptr<CoverSeed> buildCoverSeed_;

    VertexSeedStats vertexSeedStats_;
    CoverSeedStats coverSeedStats_;
    SubsumedSearchStats subsumedStats_;
    CoverSeedStats buildCoverSeedStats_;
    SubsumedSearchStats buildSubsumedStats_;
    MaintenanceTimes times_;
    UpdateCounts counts_;
    long long totalCreated_ = 0;
    long long totalRemoved_ = 0;
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

} // namespace ckm
