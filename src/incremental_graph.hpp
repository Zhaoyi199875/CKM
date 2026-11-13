#pragma once

#include "types.hpp"

#include <span>
#include <type_traits>
#include <vector>

namespace ckm {

struct AdjacencyEntry {
    Vertex neighbor = 0;
    EdgeId edgeId = invalidEdgeId;
};

class IncrementalGraph {
public:
    IncrementalGraph(int vertexCount, std::span<const Edge> edgeCatalog);
    IncrementalGraph(int vertexCount, std::vector<Edge>&& edgeCatalog);

    void applyBatch(std::span<const EdgeId> edgeIds);
    void removeBatch(std::span<const EdgeId> edgeIds);

    [[nodiscard]] int vertexCount() const noexcept;
    [[nodiscard]] int edgeCount() const noexcept;
    [[nodiscard]] int activeEdgeCount() const noexcept;
    [[nodiscard]] const Edge& edge(EdgeId edgeId) const;
    [[nodiscard]] bool isActive(EdgeId edgeId) const;
    [[nodiscard]] std::span<const Vertex> activeVertices() const noexcept;
    [[nodiscard]] std::span<const Vertex> lastActivatedVertices() const noexcept {
        return lastActivated_;
    }
    [[nodiscard]] std::span<const Vertex> lastIsolatedVertices() const noexcept {
        return lastIsolated_;
    }
    [[nodiscard]] std::span<const AdjacencyEntry> neighbors(Vertex vertex) const;
    [[nodiscard]] std::span<const AdjacencyEntry> neighborsUnchecked(
        Vertex vertex) const noexcept {
        return adjacency_[static_cast<std::size_t>(vertex)];
    }
    void prefetchListHeader(Vertex vertex) const noexcept {
        __builtin_prefetch(adjacency_.data() + static_cast<std::size_t>(vertex), 0, 1);
    }
    void prefetchListStart(Vertex vertex) const noexcept {
        __builtin_prefetch(adjacency_[static_cast<std::size_t>(vertex)].data(), 0, 1);
    }
    [[nodiscard]] int degree(Vertex vertex) const;

    template <typename Visitor>
    void forEachCommonNeighbor(Vertex u, Vertex v, Visitor&& visitor) const {
        const auto left = neighbors(u);
        const auto right = neighbors(v);
        std::size_t i = 0;
        std::size_t j = 0;
        while (i < left.size() && j < right.size()) {
            if (left[i].neighbor < right[j].neighbor) {
                ++i;
            } else if (right[j].neighbor < left[i].neighbor) {
                ++j;
            } else {
                if constexpr (std::is_same_v<
                        decltype(visitor(
                            left[i].neighbor, left[i].edgeId, right[j].edgeId)),
                        bool>) {
                    if (!visitor(
                            left[i].neighbor, left[i].edgeId, right[j].edgeId)) {
                        return;
                    }
                } else {
                    visitor(left[i].neighbor, left[i].edgeId, right[j].edgeId);
                }
                ++i;
                ++j;
            }
        }
    }

private:
    void checkVertex(Vertex vertex) const;

    std::vector<Edge> edges_;
    std::vector<std::vector<AdjacencyEntry>> adjacency_;
    std::vector<unsigned char> active_;
    std::vector<Vertex> activeVertices_;
    std::vector<Vertex> lastActivated_;
    std::vector<Vertex> lastIsolated_;
    int activeEdgeCount_ = 0;
};

}
