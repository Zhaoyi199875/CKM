#pragma once

#include "cover_seed.hpp"
#include "incremental_graph.hpp"
#include "maximal_clique_index.hpp"
#include "subsumed_search.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace ckm {

struct MaintenanceConfiguration {
    bool builds = false;
    bool sizeHistogram = false;
};

struct MaintenanceTimes {
    double coverSeed = 0.0;
    double subsumed = 0.0;
    double cliqueIndexUpdate = 0.0;
    double graphUpdate = 0.0;
    double destroyedSearch = 0.0;
    double newSearch = 0.0;
    double recovery = 0.0;
    double transientNetting = 0.0;
    double singleVertex = 0.0;
};

struct UpdateCounts {
    long long batches = 0;
    long long insertedEdges = 0;
    long long deletedEdges = 0;
    long long destroyed = 0;
    long long recovered = 0;
    long long holdingNew = 0;
    long long subsumed = 0;
    long long transient = 0;
    long long destroyedSideCalls = 0;
    long long newSideCalls = 0;
    long long singleVertexCreated = 0;
    long long singleVertexRemoved = 0;
    long long transientSingleVertex = 0;
};

struct BatchChanges {
    std::vector<Clique> created;
    std::vector<Clique> removed;
};

class MaintenanceEngine {
public:
    MaintenanceEngine(
        IncrementalGraph& graph,
        const MaintenanceConfiguration& configuration);

    void build(std::span<const EdgeId> edges);
    void updateBatch(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    void resetMeasurements();

    [[nodiscard]] const MaintenanceTimes& times() const noexcept {
        return times_;
    }
    [[nodiscard]] const UpdateCounts& updateCounts() const noexcept {
        return counts_;
    }
    [[nodiscard]] long long totalCreated() const noexcept { return totalCreated_; }
    [[nodiscard]] long long totalRemoved() const noexcept { return totalRemoved_; }
    [[nodiscard]] const CoverSeedStats& coverSeedStats() const noexcept {
        return coverSeedStats_;
    }
    [[nodiscard]] const SubsumedSearchStats& subsumedStats() const noexcept {
        return subsumedStats_;
    }
    [[nodiscard]] const MaximalCliqueIndex& index() const noexcept {
        return index_;
    }
    [[nodiscard]] std::size_t currentCliqueCount() const noexcept;
    [[nodiscard]] std::size_t isolatedVertexCount() const noexcept;

private:
    [[nodiscard]] std::vector<Clique> searchCliques(std::span<const EdgeId> batch);
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
    void updateOneWay(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    void updateHybrid(
        std::span<const EdgeId> inserted,
        std::span<const EdgeId> deleted,
        BatchChanges& changes);
    void insertIntoGraph(std::span<const EdgeId> edges);
    void removeFromGraph(std::span<const EdgeId> edges);
    void markAppeared(std::span<const Vertex> activated);
    void settleSingleVertexCliques(bool hybrid);

    IncrementalGraph& graph_;
    std::unique_ptr<CoverSeed> coverSeed_;
    MaximalCliqueIndex index_;
    std::unique_ptr<SubsumedSearch> subsumedSearch_;
    std::unique_ptr<CoverSeed> buildCoverSeed_;

    CoverSeedStats coverSeedStats_;
    SubsumedSearchStats subsumedStats_;
    CoverSeedStats buildCoverSeedStats_;
    SubsumedSearchStats buildSubsumedStats_;
    MaintenanceTimes times_;
    UpdateCounts counts_;
    long long totalCreated_ = 0;
    long long totalRemoved_ = 0;
    std::vector<unsigned char> appeared_;
    std::size_t appearedCount_ = 0U;
    bool insertedThisBatch_ = false;
    bool removedThisBatch_ = false;
    std::vector<Vertex> singleVertexRemoved_;
    std::vector<Vertex> singleVertexCreated_;
};

}
