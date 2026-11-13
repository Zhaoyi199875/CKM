#include "maintenance_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

namespace ckm {
namespace {

// Canonical cliques in a fixed order, so that the subsumed search and the
// change output see the same sequence whichever search produced them.
void normalizeResults(std::vector<Clique>& cliques) {
    for (Clique& clique : cliques) {
        clique = canonicalClique(std::move(clique));
    }
    std::ranges::sort(cliques, CliqueLexicographicOrder{});
    if (std::adjacent_find(cliques.begin(), cliques.end()) != cliques.end()) {
        throw std::logic_error("the new-clique search emitted a duplicate clique");
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

void append(std::vector<Clique>& target, std::vector<Clique>&& source) {
    target.insert(target.end(), std::make_move_iterator(source.begin()),
                  std::make_move_iterator(source.end()));
}

} // namespace

const char* newCliqueSearchName(NewCliqueSearch search) noexcept {
    switch (search) {
    case NewCliqueSearch::VertexSeed:
        return "vertexseed";
    case NewCliqueSearch::CoverSeed:
        break;
    }
    return "coverseed";
}

MaintenanceEngine::MaintenanceEngine(
    IncrementalGraph& graph,
    const MaintenanceConfiguration& configuration)
    : graph_(graph),
      configuration_(configuration),
      index_(configuration.fullMaintenance && configuration.sizeHistogram),
      appeared_(static_cast<std::size_t>(graph.vertexCount()), 0U) {
    switch (configuration.newCliqueSearch) {
    case NewCliqueSearch::VertexSeed:
        vertexSeed_ = std::make_unique<VertexSeed>(graph);
        break;
    case NewCliqueSearch::CoverSeed:
        coverSeed_ = std::make_unique<CoverSeed>(graph);
        break;
    }
    if (configuration.fullMaintenance) {
        subsumedSearch_ = std::make_unique<SubsumedSearch>(graph);
        if (configuration.builds) {
            buildCoverSeed_ = std::make_unique<CoverSeed>(graph);
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
    return configuration_.fullMaintenance
        ? index_.size() + isolatedVertexCount()
        : 0U;
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

void MaintenanceEngine::settleSingleVertexCliques(bool hybrid) {
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
    std::ranges::set_difference(
        activated, isolated, std::back_inserter(singleVertexRemoved_));
    std::ranges::set_difference(
        isolated, activated, std::back_inserter(singleVertexCreated_));
    counts_.singleVertexRemoved +=
        static_cast<long long>(singleVertexRemoved_.size());
    counts_.singleVertexCreated +=
        static_cast<long long>(singleVertexCreated_.size());
    if (hybrid) {
        // A n R: left without an edge between the halves, given one again.
        counts_.transientSingleVertex += static_cast<long long>(
            activated.size() - singleVertexRemoved_.size());
    }
}

void MaintenanceEngine::build(std::span<const EdgeId> edges) {
    graph_.applyBatch(edges);
    markAppeared(graph_.lastActivatedVertices());
    if (buildCoverSeed_ == nullptr) {
        throw std::logic_error("the engine was not configured to build");
    }
    std::vector<Clique> fresh =
        buildCoverSeed_->enumerateBatch(edges, buildCoverSeedStats_);
    normalizeResults(fresh);
    static_cast<void>(subsumedSearch_->enumerateBatch(
        edges, fresh, index_, buildSubsumedStats_));
    const std::span<const std::uint64_t> hashes =
        subsumedSearch_->newCliqueHashes();
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
    if (coverSeed_ != nullptr) {
        batchCliques = measureEnumeration(
            [&] { return coverSeed_->countBatch(batch, coverSeedStats_); },
            times_.coverSeed);
    } else {
        batchCliques = measureEnumeration(
            [&] { return vertexSeed_->countBatch(batch, vertexSeedStats_); },
            times_.vertexSeed);
    }
    totalCreated_ += batchCliques;
}

void MaintenanceEngine::updateBatch(
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
        changes.created = searchCliques(inserted);
        totalCreated_ += static_cast<long long>(changes.created.size());
        return;
    }
    ++counts_.batches;
    counts_.insertedEdges += static_cast<long long>(inserted.size());
    counts_.deletedEdges += static_cast<long long>(deleted.size());
    const bool hybrid = !inserted.empty() && !deleted.empty();
    if (hybrid) {
        updateHybrid(inserted, deleted, changes);
    } else {
        updateOneWay(inserted, deleted, changes);
    }
    // The single-vertex cliques, from the vertices the graph updates
    // activated and isolated.
    measureAccumulatedPhase(
        [&] {
            settleSingleVertexCliques(hybrid);
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

// Insertions: the new cliques holding an inserted edge are the new-clique
// search's answer on the graph after the batch, and the subsumed search finds
// the old cliques they subsume.  Deletions: the cliques the batch destroys are
// the new-clique search's answer for the deleted edges on the graph that still
// holds them, and the recovery checks the candidates the deletions free
// against the graph without them.  The index doubles as the recovery's table
// of what it already accepted.
void MaintenanceEngine::updateOneWay(
    std::span<const EdgeId> inserted,
    std::span<const EdgeId> deleted,
    BatchChanges& changes) {
    std::vector<Clique> destroyed;
    if (!deleted.empty()) {
        destroyed = searchSide(deleted, true);
        eraseDestroyed(destroyed);
        removeFromGraph(deleted);
    }
    std::vector<Clique> fresh;
    std::vector<Clique> subsumed;
    if (!inserted.empty()) {
        insertIntoGraph(inserted);
        fresh = searchSide(inserted, false);
        subsumed = searchSubsumed(inserted, fresh);
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
// was created and subsumed inside the batch (a transient clique) and is in
// neither list of its net change.
void MaintenanceEngine::updateHybrid(
    std::span<const EdgeId> inserted,
    std::span<const EdgeId> deleted,
    BatchChanges& changes) {
    std::vector<Clique> destroyed = searchSide(deleted, true);
    removeFromGraph(deleted);
    std::vector<Clique> recovered = recover(deleted, destroyed);
    eraseDestroyed(destroyed);
    insertIntoGraph(inserted);
    std::vector<Clique> fresh = searchSide(inserted, false);
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

std::vector<Clique> MaintenanceEngine::searchCliques(
    std::span<const EdgeId> batch) {
    std::vector<Clique> results;
    if (coverSeed_ != nullptr) {
        results = measureEnumeration(
            [&] { return coverSeed_->enumerateBatch(batch, coverSeedStats_); },
            times_.coverSeed);
    } else {
        results = measureEnumeration(
            [&] { return vertexSeed_->enumerateBatch(batch, vertexSeedStats_); },
            times_.vertexSeed);
    }
    normalizeResults(results);
    return results;
}

std::vector<Clique> MaintenanceEngine::searchSide(
    std::span<const EdgeId> batch,
    bool deletedSide) {
    const double timeBefore = times_.coverSeed + times_.vertexSeed;
    const long long callsBefore = searchCalls();
    std::vector<Clique> cliques = searchCliques(batch);
    const double elapsed = times_.coverSeed + times_.vertexSeed - timeBefore;
    const long long calls = searchCalls() - callsBefore;
    if (deletedSide) {
        times_.destroyedSearch += elapsed;
        counts_.destroyedSideCalls += calls;
    } else {
        times_.newSearch += elapsed;
        counts_.newSideCalls += calls;
    }
    return cliques;
}

long long MaintenanceEngine::searchCalls() const noexcept {
    return coverSeedStats_.recursiveCalls + vertexSeedStats_.recursiveCalls;
}

std::vector<Clique> MaintenanceEngine::searchSubsumed(
    std::span<const EdgeId> inserted,
    const std::vector<Clique>& fresh) {
    return measureEnumeration(
        [&] {
            return subsumedSearch_->enumerateBatch(
                inserted, fresh, index_, subsumedStats_);
        },
        times_.subsumed);
}

std::vector<Clique> MaintenanceEngine::recover(
    std::span<const EdgeId> deleted,
    const std::vector<Clique>& destroyed) {
    const auto begin = std::chrono::steady_clock::now();
    std::vector<Clique> created = subsumedSearch_->recoverBatch(
        deleted, destroyed, index_, subsumedStats_);
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
                subsumedSearch_->newCliqueHashes();
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
    vertexSeedStats_ = VertexSeedStats{};
    coverSeedStats_ = CoverSeedStats{};
    subsumedStats_ = SubsumedSearchStats{};
    times_ = MaintenanceTimes{};
    counts_ = UpdateCounts{};
    totalCreated_ = 0;
    totalRemoved_ = 0;
}

} // namespace ckm
