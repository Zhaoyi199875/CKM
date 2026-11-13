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

// Receives each clique at the moment the search emits it, so a consumer can
// process it without the producer materializing the whole batch first.  A
// default-constructed visitor is inert.
struct CliqueVisitor {
    using Function = void (*)(void*, std::span<const Vertex>);

    void* context = nullptr;
    Function function = nullptr;

    void operator()(std::span<const Vertex> clique) const {
        if (function != nullptr) {
            function(context, clique);
        }
    }

    [[nodiscard]] explicit operator bool() const noexcept {
        return function != nullptr;
    }
};

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
    explicit CliqueBatchOutput(
        bool materialize,
        CliqueVisitor visitor = {}) noexcept
        : materialize_(materialize), visitor_(visitor) {}

    void emit(std::span<const Vertex> clique) {
        ++count_;
        visitor_(clique);
        if (materialize_) {
            cliques_.push_back(
                canonicalClique(Clique(clique.begin(), clique.end())));
        }
    }

    // For producers that already hold a canonical clique: skips the
    // re-normalization `emit` has to perform on an arbitrary span.
    void emitCanonical(Clique clique) {
        ++count_;
        visitor_(clique);
        if (materialize_) {
            cliques_.push_back(std::move(clique));
        }
    }

    [[nodiscard]] long long count() const noexcept { return count_; }

    [[nodiscard]] std::vector<Clique> takeCliques() && noexcept {
        return std::move(cliques_);
    }

private:
    bool materialize_ = false;
    CliqueVisitor visitor_;
    long long count_ = 0;
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

} // namespace ckm
