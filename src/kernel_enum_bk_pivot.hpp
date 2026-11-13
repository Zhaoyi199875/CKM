#pragma once

#include "maximal_clique_index.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace ckm {

struct KernelEnumBkPivotWideCounters {
    long long searches = 0;
    long long nodes = 0;
    long long leaves = 0;
};

template <class Emit>
class KernelEnumBkPivotWide {
public:
    template <class Edges>
    KernelEnumBkPivotWide(
        std::span<const Vertex> clique,
        const Edges& newEdges,
        std::size_t minimumSize,
        KernelEnumBkPivotWideCounters& counters,
        const Emit& emit)
        : kernel_(kernelOf(newEdges)),
          kernelSize_(kernel_.size()),
          words_((kernelSize_ + 63U) / 64U),
          universalSize_(clique.size() - kernelSize_),
          cliqueSize_(clique.size()),
          minimumSize_(minimumSize),
          counters_(counters),
          emit_(emit),
          rows_(kernelSize_ * words_, 0U),
          full_(words_, 0U),
          selected_(words_, 0U),
          arena_((kernelSize_ + 1U) * 3U * words_, 0U),
          positions_(kernelSize_),
          setWords_(kernelSize_),
          masks_((clique.size() + 63U) / 64U) {
        for (std::size_t i = 0; i < kernelSize_; ++i) {
            full_[i / 64U] |= std::uint64_t{1} << (i % 64U);
        }
        for (std::size_t i = 0; i < kernelSize_; ++i) {
            std::uint64_t* rowWords = row(i);
            std::ranges::copy(full_, rowWords);
            rowWords[i / 64U] &= ~(std::uint64_t{1} << (i % 64U));
            positions_[i] = static_cast<std::size_t>(
                std::lower_bound(clique.begin(), clique.end(), kernel_[i]) -
                clique.begin());
            setWords_[i] = vertexSetWord(kernel_[i]);
        }
        for (const auto& edge : newEdges) {
            const auto u = static_cast<std::size_t>(
                std::lower_bound(kernel_.begin(), kernel_.end(), edge.u) -
                kernel_.begin());
            const auto v = static_cast<std::size_t>(
                std::lower_bound(kernel_.begin(), kernel_.end(), edge.v) -
                kernel_.begin());
            row(u)[v / 64U] &= ~(std::uint64_t{1} << (v % 64U));
            row(v)[u / 64U] &= ~(std::uint64_t{1} << (u % 64U));
        }
    }

    void run() {
        ++counters_.searches;
        std::ranges::copy(full_, slot(0U, 0U));
        visit(0U);
    }

private:
    template <class Edges>
    [[nodiscard]] static std::vector<Vertex> kernelOf(const Edges& newEdges) {
        std::vector<Vertex> kernel;
        for (const auto& edge : newEdges) {
            kernel.push_back(edge.u);
            kernel.push_back(edge.v);
        }
        std::ranges::sort(kernel);
        kernel.erase(std::unique(kernel.begin(), kernel.end()), kernel.end());
        return kernel;
    }
    [[nodiscard]] std::uint64_t* row(std::size_t index) noexcept {
        return rows_.data() + index * words_;
    }
    [[nodiscard]] std::uint64_t* slot(std::size_t depth, std::size_t which) noexcept {
        return arena_.data() + (depth * 3U + which) * words_;
    }
    [[nodiscard]] int count(const std::uint64_t* set) const noexcept {
        int total = 0;
        for (std::size_t i = 0; i < words_; ++i) total += std::popcount(set[i]);
        return total;
    }
    [[nodiscard]] int countAnd(
        const std::uint64_t* lhs,
        const std::uint64_t* rhs) const noexcept {
        int total = 0;
        for (std::size_t i = 0; i < words_; ++i) {
            total += std::popcount(lhs[i] & rhs[i]);
        }
        return total;
    }
    [[nodiscard]] bool none(const std::uint64_t* set) const noexcept {
        for (std::size_t i = 0; i < words_; ++i) {
            if (set[i] != 0U) return false;
        }
        return true;
    }

    void visit(std::size_t depth) {
        ++counters_.nodes;
        std::uint64_t* p = slot(depth, 0U);
        std::uint64_t* x = slot(depth, 1U);
        const int pCount = count(p);
        if (universalSize_ + selectedCount_ + static_cast<std::size_t>(pCount) <
            minimumSize_) {
            return;
        }
        if (pCount == 0) {
            if (none(x)) {
                ++counters_.leaves;
                emitLeaf();
            }
            return;
        }
        const bool xEmpty = none(x);
        std::size_t pivot = 0;
        int best = -1;
        bool settled = false;
        for (std::size_t i = 0; i < words_ && !settled; ++i) {
            for (std::uint64_t rest = p[i] | x[i]; rest != 0U; rest &= rest - 1U) {
                const auto index =
                    i * 64U + static_cast<std::size_t>(std::countr_zero(rest));
                const int score = countAnd(p, row(index));
                if (score == pCount) {
                    return;
                }
                if (score > best) {
                    best = score;
                    pivot = index;
                    if (xEmpty && score == pCount - 1) {
                        settled = true;
                        break;
                    }
                }
            }
        }
        std::uint64_t* extensions = slot(depth, 2U);
        const std::uint64_t* pivotRow = row(pivot);
        for (std::size_t i = 0; i < words_; ++i) {
            extensions[i] = p[i] & ~pivotRow[i];
        }
        std::uint64_t* nextP = slot(depth + 1U, 0U);
        std::uint64_t* nextX = slot(depth + 1U, 1U);
        for (std::size_t i = 0; i < words_; ++i) {
            while (extensions[i] != 0U) {
                const auto bitIndex =
                    static_cast<std::size_t>(std::countr_zero(extensions[i]));
                const std::uint64_t bit = std::uint64_t{1} << bitIndex;
                extensions[i] &= extensions[i] - 1U;
                const std::uint64_t* chosen = row(i * 64U + bitIndex);
                for (std::size_t k = 0; k < words_; ++k) {
                    nextP[k] = p[k] & chosen[k];
                    nextX[k] = x[k] & chosen[k];
                }
                selected_[i] |= bit;
                ++selectedCount_;
                visit(depth + 1U);
                selected_[i] &= ~bit;
                --selectedCount_;
                p[i] &= ~bit;
                x[i] |= bit;
            }
        }
    }

    void emitLeaf() {
        std::ranges::fill(masks_, ~std::uint64_t{0});
        if (cliqueSize_ % 64U != 0U) {
            masks_.back() = (std::uint64_t{1} << (cliqueSize_ % 64U)) - 1U;
        }
        std::uint64_t omittedHash = 0U;
        for (std::size_t i = 0; i < words_; ++i) {
            for (std::uint64_t rest = full_[i] & ~selected_[i]; rest != 0U;
                 rest &= rest - 1U) {
                const auto index =
                    i * 64U + static_cast<std::size_t>(std::countr_zero(rest));
                const std::size_t position = positions_[index];
                masks_[position / 64U] &= ~(std::uint64_t{1} << (position % 64U));
                omittedHash ^= setWords_[index];
            }
        }
        emit_(std::span<const std::uint64_t>(masks_), universalSize_ + selectedCount_,
              omittedHash);
    }

    std::vector<Vertex> kernel_;
    std::size_t kernelSize_ = 0;
    std::size_t words_ = 0;
    std::size_t universalSize_ = 0;
    std::size_t cliqueSize_ = 0;
    std::size_t minimumSize_ = 0;
    std::size_t selectedCount_ = 0;
    KernelEnumBkPivotWideCounters& counters_;
    const Emit& emit_;
    std::vector<std::uint64_t> rows_;
    std::vector<std::uint64_t> full_;
    std::vector<std::uint64_t> selected_;
    std::vector<std::uint64_t> arena_;
    std::vector<std::size_t> positions_;
    std::vector<std::uint64_t> setWords_;
    std::vector<std::uint64_t> masks_;
};

}
