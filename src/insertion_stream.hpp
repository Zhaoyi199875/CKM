#pragma once

#include "types.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace ckm {

class InsertionStream {
public:
    [[nodiscard]] static InsertionStream readEdgeList(const std::string& path);

    [[nodiscard]] int vertexCount() const noexcept;
    [[nodiscard]] int edgeCount() const noexcept;
    [[nodiscard]] const Edge& edge(EdgeId edgeId) const;
    [[nodiscard]] std::span<const Edge> edges() const noexcept;
    [[nodiscard]] std::vector<Edge> takeEdges() noexcept;
    [[nodiscard]] int originalVertexId(Vertex denseVertex) const;
    [[nodiscard]] std::span<const int> insertionTimestamps() const noexcept;

    [[nodiscard]] std::vector<EdgeBatch> makeBatches(
        std::optional<std::size_t> batchSize,
        std::optional<std::size_t> maxBatches = std::nullopt,
        std::size_t firstEdge = 0U,
        std::size_t endEdge = static_cast<std::size_t>(-1)) const;
    [[nodiscard]] std::size_t batchCount(
        std::optional<std::size_t> batchSize) const;

private:
    std::vector<int> originalVertexIds_;
    std::vector<Edge> edges_;
    std::vector<int> insertionTimestamps_;
};

}
