#pragma once
// KernelEnum's BK_pivot on G[K_C] (Section IV) for the new maximal cliques of
// more than 64 vertices.  `SubsumedSearch::kernelEnumBkPivotNarrow` is the same
// search for the other cliques, with every set in one word.
//
// For a new maximal clique C, F_C is the set of new edges inside C, the
// kernel K_C the endpoints of F_C, and the universal part U_C = C \ K_C.  Two
// kernel vertices are adjacent in G unless a new edge joins them.  The
// candidate subsumed cliques of C are U_C u I over the maximal cliques I of
// G[K_C]; every candidate keeps all of U_C, and its mask is over C.
#include "maximal_clique_index.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace ckm {

// The searches `KernelEnumBkPivotWide` ran, and their nodes and leaves.
struct KernelEnumBkPivotWideCounters {
    long long searches = 0;
    long long nodes = 0;
    long long leaves = 0;
};

// The kernel is numbered 0..|K_C|-1 and every set is a run of
// ceil(|K_C| / 64) words, so no kernel is too large.  P, X and the extension
// snapshot of each depth live in one arena (three sets per depth; a child is
// one deeper than its parent, so siblings reuse one slot), and the selected
// set is a single shared bit set kept in step with the recursion, so a node
// allocates nothing.
template <class Emit>
class KernelEnumBkPivotWide {
public:
    // C and its new edges F_C, whose endpoints are K_C.
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
          // A node at depth d has d selected vertices and a child only while
          // P is non-empty, so the depths run 0..|K_C|.
          arena_((kernelSize_ + 1U) * 3U * words_, 0U),
          positions_(kernelSize_),
          setWords_(kernelSize_),
          masks_((clique.size() + 63U) / 64U) {
        for (std::size_t i = 0; i < kernelSize_; ++i) {
            full_[i / 64U] |= std::uint64_t{1} << (i % 64U);
        }
        // A vertex's row: K_C less the vertex and its new-edge neighbours.
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

    // Emits every candidate U_C u I: a mask over C, its size, and the
    // set-hash words of the kernel vertices it leaves out.
    void run() {
        ++counters_.searches;
        std::ranges::copy(full_, slot(0U, 0U));
        visit(0U);
    }

private:
    // K_C: the endpoints of the new edges, sorted.
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
    // P (0), X (1) and the extension snapshot (2) of one depth.
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
            // Nothing left to add: the selected set is maximal exactly when
            // X is empty.
            if (none(x)) {
                ++counters_.leaves;
                emitLeaf();
            }
            return;
        }
        // Pivot: the vertex of P u X with the most neighbours in P, the first
        // in index order on ties.  Two scores end the scan early without
        // changing the pivot or anything emitted.  A vertex scoring |P| can
        // only be in X (no row holds its own vertex) and extends every set
        // below this call, so the call ends; with X empty, |P|-1 is the best
        // score left, so the first vertex reaching it is the pivot the full
        // scan would keep.
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
        // The extension set is a snapshot: P and X change as the loop goes.
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

    // The candidate U_C u I: the mask over C with the kernel vertices outside
    // I cleared, its size, and the set-hash words of the cleared vertices.
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

    std::vector<Vertex> kernel_;     // K_C, sorted
    std::size_t kernelSize_ = 0;     // |K_C|
    std::size_t words_ = 0;
    std::size_t universalSize_ = 0;  // |U_C| = |C| - |K_C|
    std::size_t cliqueSize_ = 0;     // |C|
    std::size_t minimumSize_ = 0;
    std::size_t selectedCount_ = 0;
    KernelEnumBkPivotWideCounters& counters_;
    const Emit& emit_;
    std::vector<std::uint64_t> rows_;      // adjacency in G[K_C], |K_C| rows
    std::vector<std::uint64_t> full_;
    std::vector<std::uint64_t> selected_;
    std::vector<std::uint64_t> arena_;
    std::vector<std::size_t> positions_;   // kernel index -> position in C
    std::vector<std::uint64_t> setWords_;  // kernel index -> set-hash word
    std::vector<std::uint64_t> masks_;     // one mask over C, reused
};

} // namespace ckm
