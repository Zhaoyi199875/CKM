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
    // Takes ownership instead of copying: at these graph sizes the duplicate
    // catalog is hundreds of megabytes.
    IncrementalGraph(int vertexCount, std::vector<Edge>&& edgeCatalog);

    void applyBatch(std::span<const EdgeId> edgeIds);
    // The reverse: every edge must be active and listed once.  An edge keeps
    // its id and may be applied again later; a vertex left without edges
    // drops out of `activeVertices`.
    void removeBatch(std::span<const EdgeId> edgeIds);

    [[nodiscard]] int vertexCount() const noexcept;
    [[nodiscard]] int edgeCount() const noexcept;
    [[nodiscard]] int activeEdgeCount() const noexcept;
    [[nodiscard]] const Edge& edge(EdgeId edgeId) const;
    [[nodiscard]] bool isActive(EdgeId edgeId) const;
    [[nodiscard]] std::span<const Vertex> activeVertices() const noexcept;
    // By-products of keeping `activeVertices`, each ascending: the vertices
    // the last applyBatch gave their first edge (none before it), and those
    // the last removeBatch left without an edge.
    [[nodiscard]] std::span<const Vertex> lastActivatedVertices() const noexcept {
        return lastActivated_;
    }
    [[nodiscard]] std::span<const Vertex> lastIsolatedVertices() const noexcept {
        return lastIsolated_;
    }
    [[nodiscard]] std::span<const AdjacencyEntry> neighbors(Vertex vertex) const;
    // `neighbors` without the range check, for a loop over vertices that
    // came out of the graph itself.
    [[nodiscard]] std::span<const AdjacencyEntry> neighborsUnchecked(
        Vertex vertex) const noexcept {
        return adjacency_[static_cast<std::size_t>(vertex)];
    }
    // For a loop that knows its vertices ahead: starts fetching a vertex's
    // list header, and, once that is in, the list's first entries.
    void prefetchListHeader(Vertex vertex) const noexcept {
        __builtin_prefetch(adjacency_.data() + static_cast<std::size_t>(vertex), 0, 1);
    }
    void prefetchListStart(Vertex vertex) const noexcept {
        __builtin_prefetch(adjacency_[static_cast<std::size_t>(vertex)].data(), 0, 1);
    }
    [[nodiscard]] int degree(Vertex vertex) const;

    // Both adjacency lists are sorted, so a straight two-pointer merge over
    // them is optimal for this access pattern: purely sequential, prefetcher
    // friendly. Common neighbours are visited in increasing vertex order,
    // which callers rely on to keep their candidate lists sorted.
    //
    // A visitor returning void sees every common neighbour. One returning
    // bool stops the merge the first time it returns false, which lets a
    // caller that only needs to know whether the count exceeds a small bound
    // quit as soon as it does instead of walking both lists to the end. The
    // `if constexpr` compiles the unused half away, so void visitors keep the
    // exact loop they had before.
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

} // namespace ckm
