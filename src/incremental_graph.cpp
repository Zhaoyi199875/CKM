#include "incremental_graph.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace immce {

IncrementalGraph::IncrementalGraph(
    int vertexCount,
    std::span<const Edge> edgeCatalog)
    : IncrementalGraph(
          vertexCount, std::vector<Edge>(edgeCatalog.begin(), edgeCatalog.end())) {}

IncrementalGraph::IncrementalGraph(
    int vertexCount,
    std::vector<Edge>&& edgeCatalog)
    : edges_(std::move(edgeCatalog)),
      adjacency_(vertexCount < 0 ? 0U : static_cast<std::size_t>(vertexCount)),
      active_(edges_.size(), 0) {
    if (vertexCount < 0) {
        throw std::invalid_argument("vertex count must not be negative");
    }
    for (const Edge& value : edges_) {
        checkVertex(value.u);
        checkVertex(value.v);
        if (value.u >= value.v) {
            throw std::invalid_argument("edge catalog must be normalized and loop-free");
        }
    }
    // Order a permutation rather than a copy of the catalog: on the large
    // graphs the copy alone is hundreds of megabytes, and this check runs
    // while the catalog itself is still live.
    std::vector<EdgeId> order(edges_.size());
    std::iota(order.begin(), order.end(), 0);
    std::ranges::sort(order, [this](EdgeId lhs, EdgeId rhs) {
        return edges_[static_cast<std::size_t>(lhs)] <
            edges_[static_cast<std::size_t>(rhs)];
    });
    for (std::size_t index = 1; index < order.size(); ++index) {
        if (edges_[static_cast<std::size_t>(order[index - 1U])] ==
            edges_[static_cast<std::size_t>(order[index])]) {
            throw std::invalid_argument(
                "edge catalog contains duplicate undirected edges");
        }
    }
}

void IncrementalGraph::applyBatch(std::span<const EdgeId> edgeIds) {
    std::vector<EdgeId> sortedIds(edgeIds.begin(), edgeIds.end());
    std::ranges::sort(sortedIds);
    if (std::adjacent_find(sortedIds.begin(), sortedIds.end()) != sortedIds.end()) {
        throw std::logic_error("an insertion edge must occur exactly once");
    }
    for (EdgeId id : sortedIds) {
        if (id < 0 || static_cast<std::size_t>(id) >= edges_.size()) {
            throw std::out_of_range("batch contains an invalid edge ID");
        }
        const auto index = static_cast<std::size_t>(id);
        if (active_[index] != 0) {
            throw std::logic_error("an insertion edge must occur exactly once");
        }
    }

    struct PendingEntry {
        Vertex vertex = 0;
        AdjacencyEntry adjacency;
    };
    if (edgeIds.size() > std::numeric_limits<std::size_t>::max() / 2U) {
        throw std::length_error("batch is too large");
    }
    std::vector<PendingEntry> pending;
    pending.reserve(edgeIds.size() * 2U);
    for (EdgeId id : edgeIds) {
        const Edge& value = edge(id);
        pending.push_back({value.u, {value.v, id}});
        pending.push_back({value.v, {value.u, id}});
    }
    std::ranges::sort(pending, [](const PendingEntry& lhs, const PendingEntry& rhs) {
        if (lhs.vertex != rhs.vertex) {
            return lhs.vertex < rhs.vertex;
        }
        return lhs.adjacency.neighbor < rhs.adjacency.neighbor;
    });

    std::size_t begin = 0;
    std::vector<Vertex>& newlyActiveVertices = lastActivated_;
    newlyActiveVertices.clear();
    while (begin < pending.size()) {
        std::size_t end = begin + 1U;
        while (end < pending.size() && pending[end].vertex == pending[begin].vertex) {
            ++end;
        }
        auto& old = adjacency_[static_cast<std::size_t>(pending[begin].vertex)];
        if (old.empty()) {
            newlyActiveVertices.push_back(pending[begin].vertex);
        }
        std::vector<AdjacencyEntry> merged;
        merged.reserve(old.size() + end - begin);
        std::size_t oldPosition = 0;
        std::size_t newPosition = begin;
        while (oldPosition < old.size() && newPosition < end) {
            if (old[oldPosition].neighbor < pending[newPosition].adjacency.neighbor) {
                merged.push_back(old[oldPosition++]);
            } else {
                merged.push_back(pending[newPosition++].adjacency);
            }
        }
        merged.insert(merged.end(), old.begin() + static_cast<std::ptrdiff_t>(oldPosition), old.end());
        while (newPosition < end) {
            merged.push_back(pending[newPosition++].adjacency);
        }
        old.swap(merged);
        begin = end;
    }
    if (!newlyActiveVertices.empty()) {
        std::vector<Vertex> mergedVertices;
        mergedVertices.reserve(activeVertices_.size() + newlyActiveVertices.size());
        std::ranges::merge(
            activeVertices_,
            newlyActiveVertices,
            std::back_inserter(mergedVertices));
        activeVertices_.swap(mergedVertices);
    }
    for (EdgeId id : edgeIds) {
        active_[static_cast<std::size_t>(id)] = 1;
    }
    activeEdgeCount_ += static_cast<int>(edgeIds.size());
}

void IncrementalGraph::removeBatch(std::span<const EdgeId> edgeIds) {
    std::vector<EdgeId> sortedIds(edgeIds.begin(), edgeIds.end());
    std::ranges::sort(sortedIds);
    if (std::adjacent_find(sortedIds.begin(), sortedIds.end()) != sortedIds.end()) {
        throw std::logic_error("a deleted edge must occur exactly once");
    }
    for (EdgeId id : sortedIds) {
        if (id < 0 || static_cast<std::size_t>(id) >= edges_.size()) {
            throw std::out_of_range("batch contains an invalid edge ID");
        }
        if (active_[static_cast<std::size_t>(id)] == 0) {
            throw std::logic_error("a deleted edge must be active");
        }
    }

    // Both directions of every edge, grouped by the vertex that loses it, so
    // each affected list is compacted in one pass however many of its edges
    // go.
    std::vector<std::pair<Vertex, Vertex>> pending;
    pending.reserve(edgeIds.size() * 2U);
    for (EdgeId id : edgeIds) {
        const Edge& value = edge(id);
        pending.emplace_back(value.u, value.v);
        pending.emplace_back(value.v, value.u);
    }
    std::ranges::sort(pending);

    std::vector<Vertex>& isolatedVertices = lastIsolated_;
    isolatedVertices.clear();
    std::size_t begin = 0;
    while (begin < pending.size()) {
        const Vertex vertex = pending[begin].first;
        std::size_t end = begin + 1U;
        while (end < pending.size() && pending[end].first == vertex) {
            ++end;
        }
        auto& list = adjacency_[static_cast<std::size_t>(vertex)];
        std::size_t removal = begin;
        std::size_t kept = 0;
        for (std::size_t position = 0; position < list.size(); ++position) {
            if (removal < end && list[position].neighbor == pending[removal].second) {
                ++removal;
                continue;
            }
            list[kept++] = list[position];
        }
        if (removal != end) {
            throw std::logic_error("adjacency is missing a deleted edge");
        }
        list.resize(kept);
        if (list.empty()) {
            isolatedVertices.push_back(vertex);
        }
        begin = end;
    }
    if (!isolatedVertices.empty()) {
        std::vector<Vertex> remaining;
        remaining.reserve(activeVertices_.size() - isolatedVertices.size());
        std::ranges::set_difference(
            activeVertices_, isolatedVertices, std::back_inserter(remaining));
        activeVertices_.swap(remaining);
    }
    for (EdgeId id : edgeIds) {
        active_[static_cast<std::size_t>(id)] = 0;
    }
    activeEdgeCount_ -= static_cast<int>(edgeIds.size());
}

int IncrementalGraph::vertexCount() const noexcept {
    return static_cast<int>(adjacency_.size());
}

int IncrementalGraph::edgeCount() const noexcept {
    return static_cast<int>(edges_.size());
}

int IncrementalGraph::activeEdgeCount() const noexcept {
    return activeEdgeCount_;
}

const Edge& IncrementalGraph::edge(EdgeId edgeIdValue) const {
    return edges_.at(static_cast<std::size_t>(edgeIdValue));
}

bool IncrementalGraph::isActive(EdgeId edgeIdValue) const {
    return active_.at(static_cast<std::size_t>(edgeIdValue)) != 0;
}

std::span<const Vertex> IncrementalGraph::activeVertices() const noexcept {
    return activeVertices_;
}

std::span<const AdjacencyEntry> IncrementalGraph::neighbors(Vertex vertex) const {
    checkVertex(vertex);
    return adjacency_[static_cast<std::size_t>(vertex)];
}

int IncrementalGraph::degree(Vertex vertex) const {
    return static_cast<int>(neighbors(vertex).size());
}

EdgeId IncrementalGraph::edgeId(Vertex u, Vertex v) const {
    checkVertex(u);
    checkVertex(v);
    const auto list = neighbors(u);
    const auto position = std::lower_bound(
        list.begin(), list.end(), v,
        [](const AdjacencyEntry& entry, Vertex value) {
            return entry.neighbor < value;
        });
    return position != list.end() && position->neighbor == v
        ? position->edgeId
        : invalidEdgeId;
}

bool IncrementalGraph::hasEdge(Vertex u, Vertex v) const {
    return edgeId(u, v) != invalidEdgeId;
}

void IncrementalGraph::checkVertex(Vertex vertex) const {
    if (vertex < 0 || static_cast<std::size_t>(vertex) >= adjacency_.size()) {
        throw std::out_of_range("vertex is outside the graph");
    }
}

} // namespace immce
