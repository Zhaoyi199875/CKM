#include "maintenance_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace immce {
namespace {

void normalizeResults(std::vector<Clique>& cliques) {
    for (Clique& clique : cliques) {
        clique = canonicalClique(std::move(clique));
    }
    std::ranges::sort(cliques, CliqueLexicographicOrder{});
}

void ensureNoDuplicates(
    const std::vector<Clique>& cliques,
    std::string_view label) {
    if (std::adjacent_find(cliques.begin(), cliques.end()) != cliques.end()) {
        throw std::logic_error(
            std::string(label) + " emitted a duplicate clique");
    }
}

void acceptOrCompareResults(
    std::vector<Clique> candidate,
    std::string_view label,
    int timestamp,
    std::vector<Clique>& accepted,
    std::string& acceptedLabel,
    bool& hasAccepted) {
    // Sorted so that the order is deterministic and a second framework can be
    // compared against the first.
    normalizeResults(candidate);
    ensureNoDuplicates(candidate, label);
    if (!hasAccepted) {
        accepted = std::move(candidate);
        acceptedLabel = label;
        hasAccepted = true;
        return;
    }
    if (accepted != candidate) {
        throw std::logic_error(
            acceptedLabel + "/" + std::string(label) +
            " mismatch at timestamp " + std::to_string(timestamp));
    }
}

template <typename Callable>
[[nodiscard]] auto measureEnumeration(
    Callable&& enumerate,
    double& elapsedMilliseconds) {
    const auto begin = std::chrono::steady_clock::now();
    auto result = enumerate();
    elapsedMilliseconds += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    return result;
}

template <typename Callable>
void measureAccumulatedPhase(
    Callable&& operation,
    double& elapsedMilliseconds) {
    const auto begin = std::chrono::steady_clock::now();
    operation();
    elapsedMilliseconds += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
}

// Drives full maintenance from inside the Lambda_new search: every new
// maximal clique is handed over the moment it is emitted, so the batch is
// never collected into a vector.  A new clique enters the index only after
// its own subsumed search has run, and no other new clique of the batch can
// be a proper subset of it, so the interleaving cannot erase one of them.
// That argument is about maximality in G+H alone, so it holds for whichever
// framework produced the clique.
class StreamingMaintenanceConsumer {
public:
    StreamingMaintenanceConsumer(
        int timestamp,
        const StreamTarget& target,
        SubsumedCliqueEnumerator& subsumedEnumerator,
        MaximalCliqueIndex& maximalCliques,
        SubsumedCliqueSearchStats& subsumedStats,
        SubsumedCliqueAlgorithm subsumedAlgorithm,
        double& outputMilliseconds)
        : timestamp_(timestamp),
          stream_(target.stream),
          output_(target.output),
          writeNew_(target.writeNew),
          writeSubsumed_(target.writeSubsumed),
          subsumedEnumerator_(subsumedEnumerator),
          maximalCliques_(maximalCliques),
          subsumedStats_(subsumedStats),
          subsumedAlgorithm_(subsumedAlgorithm),
          outputMilliseconds_(outputMilliseconds),
          subsumedOutput_(
              false,
              CliqueVisitor{
                  this, &StreamingMaintenanceConsumer::visitSubsumed}) {}

    [[nodiscard]] CliqueVisitor newCliqueVisitor() noexcept {
        return CliqueVisitor{this, &StreamingMaintenanceConsumer::visitNew};
    }

    [[nodiscard]] long long subsumedCount() const noexcept {
        return subsumedOutput_.count();
    }

    void finishBatch() {
        subsumedEnumerator_.finishBatch(
            maximalCliques_, subsumedStats_, subsumedOutput_);
    }

    // The batch's subsumed single-vertex cliques, after its other changes.
    void writeSubsumedSingleVertex(std::span<const Vertex> vertices) {
        if (!writeSubsumed_) {
            return;
        }
        for (const Vertex& vertex : vertices) {
            writeTimedClique("del", std::span<const Vertex>(&vertex, 1U));
        }
    }

private:
    static void visitNew(void* context, std::span<const Vertex> clique) {
        static_cast<StreamingMaintenanceConsumer*>(context)->consumeNew(clique);
    }

    static void visitSubsumed(
        void* context,
        std::span<const Vertex> clique) {
        auto& self = *static_cast<StreamingMaintenanceConsumer*>(context);
        if (self.writeSubsumed_) {
            self.writeTimedClique("del", clique);
        }
    }

    void consumeNew(std::span<const Vertex> rawClique) {
        Clique clique =
            canonicalClique(Clique(rawClique.begin(), rawClique.end()));
        if (writeNew_) {
            writeTimedClique("new", clique);
        }
        // Deliberately untimed: splitting search from maintenance here would
        // cost two clock reads per new clique, which measured more than the
        // batch materialization this pipeline removes.  The whole visit is
        // timed once per batch instead, and reported as
        // `search-and-maintenance-ms`.
        subsumedEnumerator_.consumeNewClique(
            clique,
            maximalCliques_,
            subsumedAlgorithm_,
            subsumedStats_,
            subsumedOutput_);
        // The subsumed search hashed the clique on its way through, and the
        // insertion reuses that hash.
        if (!maximalCliques_.insertHashed(
                clique, subsumedEnumerator_.newCliqueHashes().back())) {
            throw std::logic_error(
                "new clique already exists in the current maximal-clique index");
        }
    }

    void writeTimedClique(
        std::string_view change,
        std::span<const Vertex> clique) {
        const auto begin = std::chrono::steady_clock::now();
        output_ << timestamp_ << '\t' << change << '\t';
        for (std::size_t index = 0; index < clique.size(); ++index) {
            if (index != 0U) {
                output_ << ' ';
            }
            output_ << stream_.originalVertexId(clique[index]);
        }
        output_ << '\n';
        outputMilliseconds_ +=
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - begin).count();
    }

    int timestamp_ = 0;
    const InsertionStream& stream_;
    std::ostream& output_;
    bool writeNew_ = false;
    bool writeSubsumed_ = false;
    SubsumedCliqueEnumerator& subsumedEnumerator_;
    MaximalCliqueIndex& maximalCliques_;
    SubsumedCliqueSearchStats& subsumedStats_;
    SubsumedCliqueAlgorithm subsumedAlgorithm_;
    double& outputMilliseconds_;
    CliqueBatchOutput subsumedOutput_;
};

void append(std::vector<Clique>& target, std::vector<Clique>&& source) {
    target.insert(target.end(), std::make_move_iterator(source.begin()),
                  std::make_move_iterator(source.end()));
}

} // namespace

const char* mixedStrategyName(MixedStrategy strategy) noexcept {
    return strategy == MixedStrategy::TwoPhase ? "two-phase" : "single-pass";
}

MaintenanceEngine::MaintenanceEngine(
    IncrementalGraph& graph,
    const MaintenanceConfiguration& configuration)
    : graph_(graph),
      configuration_(configuration),
      // Only the framework that reads the histogram pays to maintain it.
      // SOMEi keeps its own maximal cliques and leaves this index empty.
      index_(configuration.fullMaintenance && !configuration.runSomei &&
             configuration.sizeHistogram),
      appeared_(static_cast<std::size_t>(graph.vertexCount()), 0U) {
    if (configuration.runEdge) {
        edgeEnumerator_ = std::make_unique<NewCliqueEnumerator>(graph);
    }
    if (configuration.runVertex) {
        VertexSearchConfiguration vertexConfiguration;
        vertexConfiguration.branching = configuration.vertexBranching;
        vertexConfiguration.seeds = configuration.vertexSeeds;
        vertexEnumerator_ = std::make_unique<VertexCliqueEnumerator>(
            graph, vertexConfiguration);
    }
    if (configuration.runStar) {
        StarSearchConfiguration starConfiguration;
        starConfiguration.branchScope = configuration.starBranchScope;
        starConfiguration.branching = configuration.starBranching;
        starConfiguration.rootFilter = configuration.starRootFilter;
        starEnumerator_ =
            std::make_unique<StarCliqueEnumerator>(graph, starConfiguration);
    }
    if (configuration.runSomei) {
        someiIndex_ = std::make_unique<SomeiIndex>(graph, configuration.somei);
    }
    if (configuration.fullMaintenance && !configuration.runSomei) {
        subsumedEnumerator_ = std::make_unique<SubsumedCliqueEnumerator>(graph);
        if (configuration.builds) {
            buildEnumerator_ = std::make_unique<StarCliqueEnumerator>(graph);
        }
    }
}

std::size_t MaintenanceEngine::isolatedVertexCount() const noexcept {
    // Only full maintenance tracks which vertices have appeared.
    const std::size_t active = graph_.activeVertices().size();
    return configuration_.fullMaintenance && appearedCount_ > active
        ? appearedCount_ - active
        : 0U;
}

std::size_t MaintenanceEngine::currentCliqueCount() const noexcept {
    if (someiIndex_ != nullptr) {
        return someiIndex_->maximalCliqueCount() + isolatedVertexCount();
    }
    return configuration_.fullMaintenance
        ? index_.size() + isolatedVertexCount()
        : 0U;
}

bool MaintenanceEngine::listsSingleVertexCliques() const noexcept {
    return configuration_.fullMaintenance &&
        (someiIndex_ == nullptr || someiIndex_->configuration().trackChanges);
}

void MaintenanceEngine::insertIntoGraph(std::span<const EdgeId> edges) {
    measureAccumulatedPhase(
        [&] { graph_.applyBatch(edges); }, times_.graphUpdate);
    insertedThisBatch_ = true;
}

void MaintenanceEngine::removeFromGraph(std::span<const EdgeId> edges) {
    measureAccumulatedPhase(
        [&] { graph_.removeBatch(edges); }, times_.graphUpdate);
    removedThisBatch_ = true;
}

void MaintenanceEngine::markAppeared(std::span<const Vertex> activated) {
    for (const Vertex vertex : activated) {
        unsigned char& appeared = appeared_[static_cast<std::size_t>(vertex)];
        if (appeared == 0U) {
            appeared = 1U;
            ++appearedCount_;
        }
    }
}

void MaintenanceEngine::settleSingleVertexCliques(bool twoPhase) {
    const std::span<const Vertex> activated = insertedThisBatch_
        ? graph_.lastActivatedVertices()
        : std::span<const Vertex>{};
    const std::span<const Vertex> isolated = removedThisBatch_
        ? graph_.lastIsolatedVertices()
        : std::span<const Vertex>{};
    insertedThisBatch_ = false;
    removedThisBatch_ = false;
    markAppeared(activated);
    singleVertexRemoved_.clear();
    singleVertexCreated_.clear();
    if (!listsSingleVertexCliques()) {
        return;
    }
    std::ranges::set_difference(
        activated, isolated, std::back_inserter(singleVertexRemoved_));
    std::ranges::set_difference(
        isolated, activated, std::back_inserter(singleVertexCreated_));
    counts_.singleVertexRemoved +=
        static_cast<long long>(singleVertexRemoved_.size());
    counts_.singleVertexCreated +=
        static_cast<long long>(singleVertexCreated_.size());
    if (twoPhase) {
        // A n R: left without an edge between the halves, given one again.
        counts_.transientSingleVertex += static_cast<long long>(
            activated.size() - singleVertexRemoved_.size());
    }
}

void MaintenanceEngine::applyOnly(std::span<const EdgeId> edges) {
    measureAccumulatedPhase(
        [&] { graph_.applyBatch(edges); }, times_.graphUpdate);
}

void MaintenanceEngine::build(std::span<const EdgeId> edges) {
    graph_.applyBatch(edges);
    markAppeared(graph_.lastActivatedVertices());
    if (someiIndex_ != nullptr) {
        // A build lists no changes, so SOMEi records none: in one batch
        // (--somei-build=one-batch) every maximal clique would be one.  Its
        // steps may enumerate G[R] (--somei-build-update=bk, somei.hpp): the
        // same state, so the measured batches, which never do, start from
        // what the paper's construction leaves.
        std::vector<Clique> added;
        std::vector<Clique> removed;
        const bool tracking = someiIndex_->configuration().trackChanges;
        someiIndex_->setChangeTracking(false);
        someiIndex_->setBuilding(true);
        someiIndex_->applyBatch(edges, buildSomeiStats_, added, removed);
        someiIndex_->setBuilding(false);
        someiIndex_->setChangeTracking(tracking);
        return;
    }
    if (buildEnumerator_ == nullptr) {
        throw std::logic_error("the engine was not configured to build");
    }
    std::vector<Clique> fresh =
        buildEnumerator_->enumerateBatch(edges, buildStarStats_);
    normalizeResults(fresh);
    static_cast<void>(subsumedEnumerator_->enumerateBatch(
        edges, fresh, index_, buildSubsumedStats_,
        SubsumedCliqueAlgorithm::ResidualKernel));
    const std::span<const std::uint64_t> hashes =
        subsumedEnumerator_->newCliqueHashes();
    for (std::size_t position = 0; position < fresh.size(); ++position) {
        if (!index_.insertHashed(fresh[position], hashes[position])) {
            throw std::logic_error(
                "new clique already exists in the current maximal-clique index");
        }
    }
}

void MaintenanceEngine::countBatch(std::span<const EdgeId> batch) {
    measureAccumulatedPhase(
        [&] { graph_.applyBatch(batch); }, times_.graphUpdate);
    long long batchCliques = 0;
    if (starEnumerator_ != nullptr) {
        batchCliques = measureEnumeration(
            [&] { return starEnumerator_->countBatch(batch, starStats_); },
            times_.star);
    } else if (edgeEnumerator_ != nullptr) {
        batchCliques = measureEnumeration(
            [&] { return edgeEnumerator_->countBatch(batch, edgeStats_); },
            times_.edge);
    } else {
        batchCliques = measureEnumeration(
            [&] { return vertexEnumerator_->countBatch(batch, vertexStats_); },
            times_.vertex);
    }
    totalCreated_ += batchCliques;
}

void MaintenanceEngine::streamBatch(
    const EdgeBatch& batch,
    const StreamTarget& target) {
    insertIntoGraph(batch.edgeIds);
    StreamingMaintenanceConsumer consumer(
        batch.timestamp,
        target,
        *subsumedEnumerator_,
        index_,
        subsumedStats_,
        configuration_.subsumedAlgorithm,
        times_.streamOutput);
    const double outputBefore = times_.streamOutput;
    long long batchCliques = 0;
    double* searchMilliseconds = nullptr;
    // One timer over the whole batch, not just the visit.
    //
    // `prepareBatch` builds this batch's edge index and `finishBatch` issues
    // every probe the batch held back -- which on a large index is most of
    // them, and is where a framework that produces fewer probes actually
    // spends less.  Timing only the visit would leave exactly the difference
    // being measured outside the measurement, so the boundary is drawn here
    // and is the same for every framework.
    const auto outerBegin = std::chrono::steady_clock::now();
    subsumedEnumerator_->prepareBatch(
        batch.edgeIds, configuration_.subsumedAlgorithm, subsumedStats_);
    if (starEnumerator_ != nullptr) {
        searchMilliseconds = &times_.star;
        batchCliques = measureEnumeration(
            [&] {
                return starEnumerator_->visitBatch(
                    batch.edgeIds, starStats_, consumer.newCliqueVisitor());
            },
            times_.star);
    } else if (edgeEnumerator_ != nullptr) {
        searchMilliseconds = &times_.edge;
        batchCliques = measureEnumeration(
            [&] {
                return edgeEnumerator_->visitBatch(
                    batch.edgeIds, edgeStats_, consumer.newCliqueVisitor());
            },
            times_.edge);
    } else {
        searchMilliseconds = &times_.vertex;
        batchCliques = measureEnumeration(
            [&] {
                return vertexEnumerator_->visitBatch(
                    batch.edgeIds, vertexStats_, consumer.newCliqueVisitor());
            },
            times_.vertex);
    }
    consumer.finishBatch();
    times_.streamOuter += std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - outerBegin).count();
    // Writing change lines is untimed in the batched pipeline, so take it
    // back out here too; under --no-output it is zero and costs no clock
    // reads at all.
    const double outputDelta = times_.streamOutput - outputBefore;
    *searchMilliseconds -= outputDelta;
    times_.streamOuter -= outputDelta;
    // The vertices the batch gave their first edge: their single-vertex
    // cliques are subsumed.
    measureAccumulatedPhase(
        [&] { settleSingleVertexCliques(false); }, times_.singleVertex);
    consumer.writeSubsumedSingleVertex(singleVertexRemoved_);
    totalCreated_ += batchCliques;
    totalRemoved_ += consumer.subsumedCount() +
        static_cast<long long>(singleVertexRemoved_.size());
}

void MaintenanceEngine::updateBatch(
    int timestamp,
    std::span<const EdgeId> inserted,
    std::span<const EdgeId> deleted,
    BatchChanges& changes) {
    changes.created.clear();
    changes.removed.clear();
    if (!configuration_.fullMaintenance) {
        if (!deleted.empty()) {
            throw std::logic_error("deleting edges needs full maintenance");
        }
        measureAccumulatedPhase(
            [&] { graph_.applyBatch(inserted); }, times_.graphUpdate);
        changes.created = searchCliques(inserted, timestamp);
        totalCreated_ += static_cast<long long>(changes.created.size());
        return;
    }
    ++counts_.batches;
    counts_.insertedEdges += static_cast<long long>(inserted.size());
    counts_.deletedEdges += static_cast<long long>(deleted.size());
    const bool twoPhase = someiIndex_ == nullptr &&
        configuration_.mixedStrategy == MixedStrategy::TwoPhase &&
        !inserted.empty() && !deleted.empty();
    if (someiIndex_ != nullptr) {
        updateSomei(inserted, deleted, changes);
    } else if (twoPhase) {
        updateTwoPhase(timestamp, inserted, deleted, changes);
    } else {
        updateSinglePass(timestamp, inserted, deleted, changes);
    }
    // The single-vertex cliques, from the vertices the graph updates
    // activated and isolated.
    measureAccumulatedPhase(
        [&] {
            settleSingleVertexCliques(twoPhase);
            for (const Vertex vertex : singleVertexRemoved_) {
                changes.removed.push_back(Clique{vertex});
            }
            for (const Vertex vertex : singleVertexCreated_) {
                changes.created.push_back(Clique{vertex});
            }
        },
        times_.singleVertex);
    totalCreated_ += static_cast<long long>(changes.created.size());
    totalRemoved_ += static_cast<long long>(changes.removed.size());
}

// The batch at once.  The cliques it destroys are the new-clique search's
// answer for the deleted edges on the graph that still holds them; the new
// cliques holding an inserted edge are its answer for the inserted edges on
// the graph after the batch, and the subsumed search finds the old cliques
// they subsume; the recovery checks the candidates the deletions free against
// that same graph, so a candidate an inserted edge extends is dropped there
// (in two phases it would be a transient clique).  None of the four lists can
// meet another: a destroyed clique holds a deleted edge, a recovered one no
// inserted edge, and the subsumed and recovered candidates are proper subsets
// of new and destroyed cliques.  So the index doubles as the recovery's table
// of what it already accepted, as when deleting alone.
void MaintenanceEngine::updateSinglePass(
    int timestamp,
    std::span<const EdgeId> inserted,
    std::span<const EdgeId> deleted,
    BatchChanges& changes) {
    std::vector<Clique> destroyed;
    if (!deleted.empty()) {
        destroyed = searchSide(deleted, timestamp, true);
        eraseDestroyed(destroyed);
    }
    if (!deleted.empty()) {
        removeFromGraph(deleted);
    }
    if (!inserted.empty()) {
        insertIntoGraph(inserted);
    }
    std::vector<Clique> fresh;
    std::vector<Clique> subsumed;
    if (!inserted.empty()) {
        fresh = searchSide(inserted, timestamp, false);
        subsumed = searchSubsumed(inserted, fresh);
        // Before the recovery: its `prepareBatch` clears the hashes the
        // subsumed search left for these insertions.
        insertFresh(fresh);
    }
    std::vector<Clique> recovered;
    if (!deleted.empty()) {
        recovered = recover(deleted, destroyed);
    }
    counts_.destroyed += static_cast<long long>(destroyed.size());
    counts_.recovered += static_cast<long long>(recovered.size());
    counts_.holdingNew += static_cast<long long>(fresh.size());
    counts_.subsumed += static_cast<long long>(subsumed.size());
    changes.created = std::move(fresh);
    append(changes.created, std::move(recovered));
    changes.removed = std::move(destroyed);
    append(changes.removed, std::move(subsumed));
}

// Das et al.'s order: the deletions maintained on their own through the graph
// between the two halves of the batch, where the recovery checks maximality
// and inserts what it accepts; then the insertions, whose subsumed search
// finds those cliques again if an inserted edge extends them.  Such a clique
// was created and subsumed inside the batch and is in neither list of its
// net change.
void MaintenanceEngine::updateTwoPhase(
    int timestamp,
    std::span<const EdgeId> inserted,
    std::span<const EdgeId> deleted,
    BatchChanges& changes) {
    std::vector<Clique> destroyed = searchSide(deleted, timestamp, true);
    removeFromGraph(deleted);
    std::vector<Clique> recovered = recover(deleted, destroyed);
    eraseDestroyed(destroyed);
    insertIntoGraph(inserted);
    std::vector<Clique> fresh = searchSide(inserted, timestamp, false);
    std::vector<Clique> subsumed = searchSubsumed(inserted, fresh);
    insertFresh(fresh);
    measureAccumulatedPhase(
        [&] {
            std::ranges::sort(recovered, CliqueLexicographicOrder{});
            std::ranges::sort(subsumed, CliqueLexicographicOrder{});
            std::vector<Clique> transient;
            std::ranges::set_intersection(
                recovered, subsumed, std::back_inserter(transient),
                CliqueLexicographicOrder{});
            if (transient.empty()) {
                return;
            }
            counts_.transient += static_cast<long long>(transient.size());
            const auto withoutTransient = [&](std::vector<Clique>& cliques) {
                std::vector<Clique> kept;
                kept.reserve(cliques.size() - transient.size());
                std::ranges::set_difference(
                    cliques, transient, std::back_inserter(kept),
                    CliqueLexicographicOrder{});
                cliques = std::move(kept);
            };
            withoutTransient(recovered);
            withoutTransient(subsumed);
        },
        times_.transientNetting);
    counts_.destroyed += static_cast<long long>(destroyed.size());
    counts_.recovered += static_cast<long long>(recovered.size());
    counts_.holdingNew += static_cast<long long>(fresh.size());
    counts_.subsumed += static_cast<long long>(subsumed.size());
    changes.created = std::move(fresh);
    append(changes.created, std::move(recovered));
    changes.removed = std::move(destroyed);
    append(changes.removed, std::move(subsumed));
}

// SOMEi reads its update scope off the union of the graphs before and after
// the batch -- the graph holding both the inserted and the deleted edges --
// and rolls back and redoes its headers on the graph after it.
void MaintenanceEngine::updateSomei(
    std::span<const EdgeId> inserted,
    std::span<const EdgeId> deleted,
    BatchChanges& changes) {
    if (!inserted.empty()) {
        insertIntoGraph(inserted);
    }
    if (deleted.empty()) {
        measureAccumulatedPhase(
            [&] {
                someiIndex_->applyBatch(
                    inserted, someiStats_, changes.created, changes.removed);
            },
            times_.somei);
        return;
    }
    changedEdges_.assign(inserted.begin(), inserted.end());
    changedEdges_.insert(changedEdges_.end(), deleted.begin(), deleted.end());
    measureAccumulatedPhase(
        [&] { someiIndex_->beginBatch(changedEdges_, someiStats_); },
        times_.somei);
    removeFromGraph(deleted);
    measureAccumulatedPhase(
        [&] {
            someiIndex_->endBatch(someiStats_, changes.created, changes.removed);
        },
        times_.somei);
}

std::vector<Clique> MaintenanceEngine::searchCliques(
    std::span<const EdgeId> batch,
    int timestamp) {
    std::vector<Clique> results;
    std::string resultLabel;
    bool hasResults = false;
    if (starEnumerator_ != nullptr) {
        acceptOrCompareResults(
            measureEnumeration(
                [&] { return starEnumerator_->enumerateBatch(batch, starStats_); },
                times_.star),
            "star", timestamp, results, resultLabel, hasResults);
    }
    if (edgeEnumerator_ != nullptr) {
        acceptOrCompareResults(
            measureEnumeration(
                [&] { return edgeEnumerator_->enumerateBatch(batch, edgeStats_); },
                times_.edge),
            "edge", timestamp, results, resultLabel, hasResults);
    }
    if (vertexEnumerator_ != nullptr) {
        acceptOrCompareResults(
            measureEnumeration(
                [&] {
                    return vertexEnumerator_->enumerateBatch(batch, vertexStats_);
                },
                times_.vertex),
            "vertex", timestamp, results, resultLabel, hasResults);
    }
    if (!hasResults) {
        throw std::logic_error("no enumeration framework was selected");
    }
    return results;
}

std::vector<Clique> MaintenanceEngine::searchSide(
    std::span<const EdgeId> batch,
    int timestamp,
    bool deletedSide) {
    const double timeBefore = times_.star + times_.edge + times_.vertex;
    const long long callsBefore = producerCalls();
    std::vector<Clique> cliques = searchCliques(batch, timestamp);
    const double elapsed =
        times_.star + times_.edge + times_.vertex - timeBefore;
    const long long calls = producerCalls() - callsBefore;
    if (deletedSide) {
        times_.destroyedSearch += elapsed;
        counts_.destroyedSideCalls += calls;
    } else {
        times_.newSearch += elapsed;
        counts_.newSideCalls += calls;
    }
    return cliques;
}

long long MaintenanceEngine::producerCalls() const noexcept {
    return starStats_.recursiveCalls + edgeStats_.recursiveCalls +
        vertexStats_.recursiveCalls;
}

std::vector<Clique> MaintenanceEngine::searchSubsumed(
    std::span<const EdgeId> inserted,
    const std::vector<Clique>& fresh) {
    return measureEnumeration(
        [&] {
            return subsumedEnumerator_->enumerateBatch(
                inserted, fresh, index_, subsumedStats_,
                configuration_.subsumedAlgorithm);
        },
        times_.subsumed);
}

std::vector<Clique> MaintenanceEngine::recover(
    std::span<const EdgeId> deleted,
    const std::vector<Clique>& destroyed) {
    const auto begin = std::chrono::steady_clock::now();
    std::vector<Clique> created = subsumedEnumerator_->recoverBatch(
        deleted, destroyed, index_, subsumedStats_,
        configuration_.subsumedAlgorithm);
    const double elapsed = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    times_.subsumed += elapsed;
    times_.recovery += elapsed;
    return created;
}

void MaintenanceEngine::insertFresh(const std::vector<Clique>& fresh) {
    measureAccumulatedPhase(
        [&] {
            // The subsumed search hashed every new clique, in order, and the
            // insertions reuse those hashes.
            const std::span<const std::uint64_t> hashes =
                subsumedEnumerator_->newCliqueHashes();
            if (hashes.size() != fresh.size()) {
                throw std::logic_error(
                    "the subsumed search did not hash every new clique");
            }
            // Each insertion misses twice: on the table slot and on the
            // clique's own heap buffer, which the batch vector holds one
            // allocation per clique.  Both are known `insertLookahead`
            // cliques ahead of time.
            constexpr std::size_t insertLookahead = 24U;
            const auto prefetchInsert = [&](std::size_t index) {
                index_.prefetchForInsert(hashes[index]);
                __builtin_prefetch(fresh[index].data(), 0, 1);
            };
            for (std::size_t index = 0;
                 index < std::min(insertLookahead, fresh.size()); ++index) {
                prefetchInsert(index);
            }
            for (std::size_t index = 0; index < fresh.size(); ++index) {
                if (index + insertLookahead < fresh.size()) {
                    prefetchInsert(index + insertLookahead);
                }
                if (!index_.insertHashed(fresh[index], hashes[index])) {
                    throw std::logic_error(
                        "new clique already exists in the current "
                        "maximal-clique index");
                }
            }
        },
        times_.cliqueIndexUpdate);
}

void MaintenanceEngine::eraseDestroyed(const std::vector<Clique>& destroyed) {
    measureAccumulatedPhase(
        [&] {
            for (const Clique& clique : destroyed) {
                if (!index_.eraseCanonical(clique)) {
                    throw std::logic_error(
                        "a destroyed clique is missing from the maximal-clique "
                        "index");
                }
            }
        },
        times_.cliqueIndexUpdate);
}

void MaintenanceEngine::resetMeasurements() {
    edgeStats_ = NewCliqueSearchStats{};
    vertexStats_ = VertexSearchStats{};
    starStats_ = StarSearchStats{};
    subsumedStats_ = SubsumedCliqueSearchStats{};
    someiStats_ = SomeiStats{};
    times_ = MaintenanceTimes{};
    counts_ = UpdateCounts{};
    totalCreated_ = 0;
    totalRemoved_ = 0;
}

} // namespace immce
