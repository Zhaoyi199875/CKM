#pragma once

#include "incremental_graph.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace ckm {

struct LocalGraphStats {
    long long builds = 0;
    long long vertexOccurrences = 0;
    long long edgeOccurrences = 0;
    std::size_t maximumVertices = 0;
    std::size_t peakBytes = 0;
};

struct LocalVertexSet {
    std::vector<Vertex> vertices;
    std::vector<std::uint64_t> bits;

    [[nodiscard]] bool empty() const noexcept { return vertices.empty(); }
    [[nodiscard]] std::size_t size() const noexcept { return vertices.size(); }
};

class LocalSearchGraph {
public:
    explicit LocalSearchGraph(int globalVertexCount)
        : localIndexByVertex_(
              globalVertexCount < 0 ? 0U : static_cast<std::size_t>(globalVertexCount),
              -1) {
        if (globalVertexCount < 0) {
            throw std::invalid_argument("global vertex count must not be negative");
        }
    }

    void rebuild(
        const IncrementalGraph& graph,
        std::span<const Vertex> root,
        std::span<const Vertex> candidates,
        std::span<const Vertex> excluded,
        LocalGraphStats& stats) {
        for (Vertex vertex : localToGlobal_) {
            localIndexByVertex_[static_cast<std::size_t>(vertex)] = -1;
        }
        localToGlobal_.clear();
        localToGlobal_.reserve(root.size() + candidates.size() + excluded.size());
        localToGlobal_.insert(localToGlobal_.end(), root.begin(), root.end());
        localToGlobal_.insert(localToGlobal_.end(), candidates.begin(), candidates.end());
        localToGlobal_.insert(localToGlobal_.end(), excluded.begin(), excluded.end());
        std::ranges::sort(localToGlobal_);
        localToGlobal_.erase(
            std::unique(localToGlobal_.begin(), localToGlobal_.end()),
            localToGlobal_.end());

        const std::size_t n = localToGlobal_.size();
        for (std::size_t local = 0; local < n; ++local) {
            localIndexByVertex_[static_cast<std::size_t>(localToGlobal_[local])] =
                static_cast<int>(local);
        }
        wordsPerRow_ = (n + 63U) / 64U;
        if (wordsPerRow_ != 0 && n > std::numeric_limits<std::size_t>::max() / wordsPerRow_) {
            throw std::length_error("local adjacency matrix is too large");
        }
        adjacencyRows_.assign(n * wordsPerRow_, 0);

        if (trackDegrees_) {
            localDegrees_.assign(n, 0);
        }
        long long localEdges = 0;
        for (std::size_t local = 0; local < n; ++local) {
            const Vertex vertex = localToGlobal_[local];
            for (const AdjacencyEntry& entry : graph.neighbors(vertex)) {
                const int neighborLocalValue = localIndex(entry.neighbor);
                if (neighborLocalValue < 0) {
                    continue;
                }
                const auto neighborLocal = static_cast<std::size_t>(neighborLocalValue);
                adjacencyRows_[local * wordsPerRow_ + neighborLocal / 64U] |=
                    std::uint64_t{1} << (neighborLocal % 64U);
                if (trackDegrees_) {
                    ++localDegrees_[local];
                }
                if (local < neighborLocal) {
                    ++localEdges;
                }
            }
        }

        ++stats.builds;
        stats.vertexOccurrences += static_cast<long long>(n);
        stats.edgeOccurrences += localEdges;
        stats.maximumVertices = std::max(stats.maximumVertices, n);
        const std::size_t bytes = adjacencyRows_.size() * sizeof(std::uint64_t) +
            localToGlobal_.size() * sizeof(Vertex);
        stats.peakBytes = std::max(stats.peakBytes, bytes);
    }

    void setTrackDegrees(bool track) noexcept { trackDegrees_ = track; }
    [[nodiscard]] int localDegree(std::size_t local) const noexcept {
        return localDegrees_[local];
    }

    void rebuildFrom(
        const LocalSearchGraph& parent,
        std::span<const Vertex> vertices,
        LocalGraphStats& stats) {
        if (!std::is_sorted(vertices.begin(), vertices.end()) ||
            std::adjacent_find(vertices.begin(), vertices.end()) !=
                vertices.end()) {
            throw std::invalid_argument(
                "re-indexed vertices must be sorted and unique");
        }
        for (Vertex vertex : localToGlobal_) {
            localIndexByVertex_[static_cast<std::size_t>(vertex)] = -1;
        }
        localToGlobal_.assign(vertices.begin(), vertices.end());

        const std::size_t n = localToGlobal_.size();
        for (std::size_t local = 0; local < n; ++local) {
            localIndexByVertex_[static_cast<std::size_t>(localToGlobal_[local])] =
                static_cast<int>(local);
        }
        wordsPerRow_ = (n + 63U) / 64U;
        if (wordsPerRow_ != 0 && n > std::numeric_limits<std::size_t>::max() / wordsPerRow_) {
            throw std::length_error("local adjacency matrix is too large");
        }
        adjacencyRows_.assign(n * wordsPerRow_, 0);

        const std::size_t parentWords = parent.wordsPerRow_;
        parentLocalToCompact_.resize(parent.localToGlobal_.size(), -1);
        parentMask_.assign(parentWords, 0);
        parentLocals_.clear();
        parentLocals_.reserve(n);
        for (std::size_t local = 0; local < n; ++local) {
            const int parentLocalValue = parent.localIndex(localToGlobal_[local]);
            if (parentLocalValue < 0) {
                throw std::logic_error(
                    "re-indexed vertex is outside the parent search graph");
            }
            const auto parentLocal = static_cast<std::size_t>(parentLocalValue);
            parentLocals_.push_back(parentLocal);
            parentLocalToCompact_[parentLocal] = static_cast<int>(local);
            parentMask_[parentLocal / 64U] |=
                std::uint64_t{1} << (parentLocal % 64U);
        }

        long long localEdges = 0;
        for (std::size_t local = 0; local < n; ++local) {
            const auto parentRow = parent.rowAt(parentLocals_[local]);
            std::uint64_t* const row = adjacencyRows_.data() + local * wordsPerRow_;
            for (std::size_t wordIndex = 0; wordIndex < parentWords; ++wordIndex) {
                std::uint64_t word = parentRow[wordIndex] & parentMask_[wordIndex];
                while (word != 0U) {
                    const std::size_t neighborParentLocal = wordIndex * 64U +
                        static_cast<std::size_t>(std::countr_zero(word));
                    word &= word - 1U;
                    const auto neighborLocal = static_cast<std::size_t>(
                        parentLocalToCompact_[neighborParentLocal]);
                    row[neighborLocal / 64U] |=
                        std::uint64_t{1} << (neighborLocal % 64U);
                    if (local < neighborLocal) {
                        ++localEdges;
                    }
                }
            }
        }
        for (std::size_t parentLocal : parentLocals_) {
            parentLocalToCompact_[parentLocal] = -1;
        }

        ++stats.builds;
        stats.vertexOccurrences += static_cast<long long>(n);
        stats.edgeOccurrences += localEdges;
        stats.maximumVertices = std::max(stats.maximumVertices, n);
        const std::size_t bytes = adjacencyRows_.size() * sizeof(std::uint64_t) +
            localToGlobal_.size() * sizeof(Vertex);
        stats.peakBytes = std::max(stats.peakBytes, bytes);
    }

    void projectSetFrom(
        const LocalVertexSet& parentSet,
        LocalVertexSet& result) const {
        result.vertices.clear();
        result.bits.assign(wordsPerRow_, 0);
        for (std::size_t local = 0; local < parentLocals_.size(); ++local) {
            const std::size_t parentLocal = parentLocals_[local];
            if ((parentSet.bits[parentLocal / 64U] &
                 (std::uint64_t{1} << (parentLocal % 64U))) != 0U) {
                result.bits[local / 64U] |= std::uint64_t{1} << (local % 64U);
            }
        }
    }

    [[nodiscard]] std::size_t wordsPerRow() const noexcept { return wordsPerRow_; }

    [[nodiscard]] std::span<const Vertex> localVertices() const noexcept {
        return localToGlobal_;
    }

    [[nodiscard]] int indexOf(Vertex vertex) const noexcept {
        return localIndex(vertex);
    }

    [[nodiscard]] std::span<const std::uint64_t> rowAt(
        std::size_t local) const noexcept {
        return std::span<const std::uint64_t>(
            adjacencyRows_.data() + local * wordsPerRow_, wordsPerRow_);
    }

    void setBits(std::span<const Vertex> vertices, LocalVertexSet& set) const {
        set.vertices.clear();
        set.bits.assign(wordsPerRow_, 0);
        for (Vertex vertex : vertices) {
            const int localValue = localIndex(vertex);
            if (localValue < 0) {
                throw std::logic_error("vertex set is outside the local search graph");
            }
            const auto local = static_cast<std::size_t>(localValue);
            set.bits[local / 64U] |= std::uint64_t{1} << (local % 64U);
        }
    }

    [[nodiscard]] LocalVertexSet makeSet(std::vector<Vertex> vertices) const {
        LocalVertexSet result{std::move(vertices), {}};
        synchronize(result);
        return result;
    }

    void intersectInto(
        const LocalVertexSet& input,
        Vertex vertex,
        LocalVertexSet& result) const {
        result.bits.resize(wordsPerRow_);
        const auto adjacent = row(vertex);
        for (std::size_t word = 0; word < wordsPerRow_; ++word) {
            result.bits[word] = input.bits[word] & adjacent[word];
        }
        rebuildVertices(result);
    }

    [[nodiscard]] std::vector<Vertex> differenceVerticesFromNeighbors(
        const LocalVertexSet& input,
        Vertex vertex) const {
        std::vector<Vertex> result;
        result.reserve(std::min(input.size(), std::size_t{8U}));
        differenceVerticesInto(input, vertex, result);
        return result;
    }

    void differenceVerticesInto(
        const LocalVertexSet& input,
        Vertex vertex,
        std::vector<Vertex>& result) const {
        result.clear();
        const auto adjacent = row(vertex);
        for (std::size_t wordIndex = 0;
             wordIndex < wordsPerRow_;
             ++wordIndex) {
            std::uint64_t word = input.bits[wordIndex] & ~adjacent[wordIndex];
            while (word != 0U) {
                const int bit = std::countr_zero(word);
                const std::size_t local = wordIndex * 64U +
                    static_cast<std::size_t>(bit);
                if (local < localToGlobal_.size()) {
                    result.push_back(localToGlobal_[local]);
                }
                word &= word - 1U;
            }
        }
    }

    [[nodiscard]] int intersectionCount(
        const LocalVertexSet& input,
        Vertex vertex) const noexcept {
        const auto adjacent = row(vertex);
        std::size_t count = 0;
        for (std::size_t word = 0; word < wordsPerRow_; ++word) {
            count += static_cast<std::size_t>(std::popcount(input.bits[word] & adjacent[word]));
        }
        return static_cast<int>(count);
    }

    [[nodiscard]] bool contains(const LocalVertexSet& set, Vertex vertex) const noexcept {
        const int localValue = localIndex(vertex);
        if (localValue < 0) {
            return false;
        }
        const auto local = static_cast<std::size_t>(localValue);
        return (set.bits[local / 64U] & (std::uint64_t{1} << (local % 64U))) != 0;
    }

    void erase(LocalVertexSet& set, Vertex vertex) const {
        const int localValue = localIndex(vertex);
        if (localValue < 0) {
            return;
        }
        const auto local = static_cast<std::size_t>(localValue);
        set.bits[local / 64U] &= ~(std::uint64_t{1} << (local % 64U));
        const auto position = std::lower_bound(set.vertices.begin(), set.vertices.end(), vertex);
        if (position != set.vertices.end() && *position == vertex) {
            set.vertices.erase(position);
        }
    }

    void moveToExcludedBits(
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        Vertex vertex) const noexcept {
        const int localValue = localIndex(vertex);
        if (localValue < 0) {
            return;
        }
        const auto local = static_cast<std::size_t>(localValue);
        const std::uint64_t bit = std::uint64_t{1} << (local % 64U);
        candidates.bits[local / 64U] &= ~bit;
        excluded.bits[local / 64U] |= bit;
    }

    void insert(LocalVertexSet& set, Vertex vertex) const {
        const int localValue = localIndex(vertex);
        if (localValue < 0) {
            throw std::logic_error("vertex is outside the local search graph");
        }
        const auto local = static_cast<std::size_t>(localValue);
        set.bits[local / 64U] |= std::uint64_t{1} << (local % 64U);
        insertSortedUnique(set.vertices, vertex);
    }

    void synchronize(LocalVertexSet& set) const {
        set.bits.assign(wordsPerRow_, 0);
        for (Vertex vertex : set.vertices) {
            const int localValue = localIndex(vertex);
            if (localValue < 0) {
                throw std::logic_error("vertex set is outside the local search graph");
            }
            const auto local = static_cast<std::size_t>(localValue);
            set.bits[local / 64U] |= std::uint64_t{1} << (local % 64U);
        }
    }

private:
    [[nodiscard]] int localIndex(Vertex vertex) const noexcept {
        if (vertex < 0 || static_cast<std::size_t>(vertex) >= localIndexByVertex_.size()) {
            return -1;
        }
        return localIndexByVertex_[static_cast<std::size_t>(vertex)];
    }

    [[nodiscard]] std::span<const std::uint64_t> row(Vertex vertex) const noexcept {
        const auto local = static_cast<std::size_t>(localIndex(vertex));
        return std::span<const std::uint64_t>(
            adjacencyRows_.data() + local * wordsPerRow_, wordsPerRow_);
    }

    void rebuildVertices(LocalVertexSet& set) const {
        set.vertices.clear();
        for (std::size_t wordIndex = 0; wordIndex < set.bits.size(); ++wordIndex) {
            std::uint64_t word = set.bits[wordIndex];
            while (word != 0) {
                const int bit = std::countr_zero(word);
                const std::size_t local = wordIndex * 64U + static_cast<std::size_t>(bit);
                if (local < localToGlobal_.size()) {
                    set.vertices.push_back(localToGlobal_[local]);
                }
                word &= word - 1U;
            }
        }
    }

    std::vector<int> localIndexByVertex_;
    std::vector<Vertex> localToGlobal_;
    std::vector<std::uint64_t> adjacencyRows_;
    std::vector<int> parentLocalToCompact_;
    std::vector<std::uint64_t> parentMask_;
    std::vector<std::size_t> parentLocals_;
    std::vector<int> localDegrees_;
    bool trackDegrees_ = false;
    std::size_t wordsPerRow_ = 0;
};

}
