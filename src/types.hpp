#pragma once

#include <algorithm>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace ckm {

using Vertex = int;
using EdgeId = int;
using Clique = std::vector<Vertex>;

constexpr EdgeId invalidEdgeId = -1;

struct Edge {
    Vertex u = 0;
    Vertex v = 0;

    auto operator<=>(const Edge&) const = default;
};

[[nodiscard]] inline Edge normalizedEdge(Vertex u, Vertex v) noexcept {
    if (v < u) {
        std::swap(u, v);
    }
    return Edge{u, v};
}

struct InsertionEvent {
    EdgeId edgeId = invalidEdgeId;
    int timestamp = 0;
};

struct EdgeBatch {
    int timestamp = 0;
    std::vector<EdgeId> edgeIds;
};

[[nodiscard]] inline std::vector<Vertex> uniqueSorted(
    std::vector<Vertex> values) {
    std::ranges::sort(values);
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

inline void insertSortedUnique(std::vector<Vertex>& values, Vertex value) {
    const auto position = std::lower_bound(values.begin(), values.end(), value);
    if (position == values.end() || *position != value) {
        values.insert(position, value);
    }
}

[[nodiscard]] inline Clique canonicalClique(Clique clique) {
    return uniqueSorted(std::move(clique));
}

class CliqueBatchOutput {
public:
    void emit(std::span<const Vertex> clique) {
        cliques_.push_back(
            canonicalClique(Clique(clique.begin(), clique.end())));
    }

    void emitCanonical(Clique clique) {
        cliques_.push_back(std::move(clique));
    }

    [[nodiscard]] std::vector<Clique> takeCliques() && noexcept {
        return std::move(cliques_);
    }

private:
    std::vector<Clique> cliques_;
};

struct CliqueOutputOrder {
    [[nodiscard]] bool operator()(const Clique& lhs, const Clique& rhs) const noexcept {
        if (lhs.size() != rhs.size()) {
            return lhs.size() > rhs.size();
        }
        return lhs < rhs;
    }
};

struct CliqueLexicographicOrder {
    [[nodiscard]] bool operator()(const Clique& lhs, const Clique& rhs) const noexcept {
        return lhs < rhs;
    }
};

}
