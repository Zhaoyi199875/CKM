#pragma once

#include "incremental_graph.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace immce {

// Each batch edge's position in the batch, looked up by edge id.  Both
// new-clique searches order ownership by it.  `begin` marks one batch and
// hands back a scope that unmarks exactly the entries it set, so a batch pays
// for its own edges only, not for the |E|-sized table.
class BatchIndex {
public:
    class Scope {
    public:
        explicit Scope(BatchIndex& owner) noexcept : owner_(&owner) {}

        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;

        Scope(Scope&& other) noexcept
            : owner_(std::exchange(other.owner_, nullptr)) {}

        Scope& operator=(Scope&& other) noexcept {
            if (this != &other) {
                reset();
                owner_ = std::exchange(other.owner_, nullptr);
            }
            return *this;
        }

        ~Scope() { reset(); }

    private:
        void reset() noexcept {
            if (owner_ != nullptr) {
                owner_->clear();
                owner_ = nullptr;
            }
        }

        BatchIndex* owner_ = nullptr;
    };

    explicit BatchIndex(const IncrementalGraph& graph)
        : graph_(graph),
          edgeRankById_(static_cast<std::size_t>(graph.edgeCount()), -1) {}

    BatchIndex(const BatchIndex&) = delete;
    BatchIndex& operator=(const BatchIndex&) = delete;
    BatchIndex(BatchIndex&&) = delete;
    BatchIndex& operator=(BatchIndex&&) = delete;

    [[nodiscard]] Scope begin(std::span<const EdgeId> batch) {
        if (marked_) {
            throw std::logic_error("a batch index scope is already active");
        }
        marked_ = true;
        try {
            if (batch.size() > static_cast<std::size_t>(graph_.edgeCount())) {
                throw std::invalid_argument("batch contains too many edge IDs");
            }
            for (std::size_t position = 0; position < batch.size(); ++position) {
                const EdgeId edgeId = batch[position];
                if (edgeId < 0 || edgeId >= graph_.edgeCount() ||
                    !graph_.isActive(edgeId)) {
                    throw std::invalid_argument(
                        "batch edge must already be active in the graph");
                }
                int& rank = edgeRankById_[static_cast<std::size_t>(edgeId)];
                if (rank >= 0) {
                    throw std::invalid_argument("batch contains a duplicate edge ID");
                }
                rank = static_cast<int>(position);
                touchedEdges_.push_back(edgeId);
            }
        } catch (...) {
            clear();
            throw;
        }
        return Scope(*this);
    }

    // The edge's position in the current batch, or -1 when it is not a batch
    // edge (or is `invalidEdgeId`).
    [[nodiscard]] int edgeRank(EdgeId edgeId) const noexcept {
        return edgeId == invalidEdgeId
            ? -1
            : edgeRankById_[static_cast<std::size_t>(edgeId)];
    }

private:
    void clear() noexcept {
        for (EdgeId edgeId : touchedEdges_) {
            edgeRankById_[static_cast<std::size_t>(edgeId)] = -1;
        }
        touchedEdges_.clear();
        marked_ = false;
    }

    const IncrementalGraph& graph_;
    std::vector<int> edgeRankById_;
    std::vector<EdgeId> touchedEdges_;
    bool marked_ = false;
};

} // namespace immce
