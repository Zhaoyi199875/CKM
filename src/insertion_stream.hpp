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
    // The stream is in increasing timestamp order, and edges that share a
    // timestamp keep the order the input file lists them in (first occurrence
    // of each edge): a stable sort by timestamp, which is all "increasing
    // timestamp order" asks for. On datasets whose timestamps are days or
    // years nearly every edge has such ties, so this order decides which edges
    // a fixed-size batch holds; the datasets list one event's edges together
    // -- a paper's author clique, a page's links -- and it keeps them
    // together.
    [[nodiscard]] static InsertionStream readEdgeList(const std::string& path);

    [[nodiscard]] int vertexCount() const noexcept;
    [[nodiscard]] int edgeCount() const noexcept;
    [[nodiscard]] const Edge& edge(EdgeId edgeId) const;
    [[nodiscard]] std::span<const Edge> edges() const noexcept;
    // Hands the catalog to the graph instead of copying it.  At these graph
    // sizes the duplicate is hundreds of megabytes, and nothing reads the
    // stream's copy afterwards.
    [[nodiscard]] std::vector<Edge> takeEdges() noexcept;
    [[nodiscard]] int originalVertexId(Vertex denseVertex) const;
    // Edge ids are assigned in insertion order, so an event is just its
    // timestamp; only the timestamps are stored.
    [[nodiscard]] std::span<const int> insertionTimestamps() const noexcept;

    // A positive value chunks the globally time-ordered insertion stream by
    // edge count, matching the paper's batch model. nullopt is the explicit
    // compatibility mode that groups equal creation timestamps.  `maxBatches`
    // stops early so a truncated run never materializes the whole stream.
    // `firstEdge` and `endEdge` chunk only the edges [firstEdge, endEdge) --
    // what a sliding window's first window and the measured run after it
    // need; they require a
    // fixed batch size.
    [[nodiscard]] std::vector<EdgeBatch> makeBatches(
        std::optional<std::size_t> batchSize,
        std::optional<std::size_t> maxBatches = std::nullopt,
        std::size_t firstEdge = 0U,
        std::size_t endEdge = static_cast<std::size_t>(-1)) const;
    // How many batches the stream would yield without truncation.
    [[nodiscard]] std::size_t batchCount(
        std::optional<std::size_t> batchSize) const;

private:
    std::vector<int> originalVertexIds_;
    std::vector<Edge> edges_;
    std::vector<int> insertionTimestamps_;
};

} // namespace ckm
