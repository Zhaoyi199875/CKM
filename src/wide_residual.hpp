#pragma once
// Pivot BK on the complement of the internal batch edges, over either the
// full clique or just the endpoints of those edges. Vertices outside the
// search domain are fixed in every emitted candidate.
//
// For a new maximal clique C with internal batch edges H[C], let T be the
// endpoints of H[C] and K = C \ T the mandatory core.  F = (T, H[C]) is the
// conflict graph; a candidate keeps an independent set of F, so the candidates
// are K union S over the maximal independent sets S of F, which are the
// maximal cliques of complement(F), enumerated here by Bron-Kerbosch.
//
//   * `runNarrow` is the single-word search, for |T| <= 64.
//   * `Bits<W>` is a fixed-width bitset over W 64-bit words, W in 2..4, so
//     65 <= |T| <= 256 is covered by `WideSearch`.
//   * `DynamicSearch` takes every wider kernel, with the word count set at run
//     time, so no kernel is ever too wide to search.
//   * The recursion carries P, X and the selected set by value in W words
//     (`DynamicSearch`: in an arena, three sets per depth); nothing is heap
//     allocated per node.
//
// |T| and |C| are separate everywhere.  The emitted mask is over C, which can
// pass 64 vertices while |T| stays small, and the other way round.
#include "maximal_clique_index.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace immce::wideresidual {

// Search counters: nodes and leaves over every search, and the
// search calls by the path they took.
struct Counters {
    long long nodes = 0;
    long long leaves = 0;
    long long kernelCalls = 0;
    long long kernelNarrow = 0;
    long long kernelWide = 0;
    long long kernelDynamic = 0;
    long long dynamicNodes = 0;
};

// The fixed-width search's limit; wider kernels take `DynamicSearch`.
inline constexpr std::size_t maximumWords = 4U;
inline constexpr std::size_t maximumWidth = 64U * maximumWords;

// A fixed-width bit set over W words.  Every operation is word-wise, so the
// tail word's unused bits stay zero as long as the set is only ever built from
// `prefix` and from operations on sets that already satisfy that.
template <std::size_t W>
struct Bits {
    std::array<std::uint64_t, W> word{};

    [[nodiscard]] static Bits prefix(std::size_t bits) noexcept {
        Bits value{};
        for (std::size_t i = 0; i < W; ++i) {
            const std::size_t low = i * 64U;
            if (low >= bits) break;
            const std::size_t width = bits - low >= 64U ? 64U : bits - low;
            value.word[i] = width == 64U
                ? ~std::uint64_t{0}
                : (std::uint64_t{1} << width) - 1U;
        }
        return value;
    }

    void set(std::size_t index) noexcept {
        word[index / 64U] |= std::uint64_t{1} << (index % 64U);
    }
    void reset(std::size_t index) noexcept {
        word[index / 64U] &= ~(std::uint64_t{1} << (index % 64U));
    }
    [[nodiscard]] bool test(std::size_t index) const noexcept {
        return ((word[index / 64U] >> (index % 64U)) & 1U) != 0U;
    }
    [[nodiscard]] bool none() const noexcept {
        for (std::uint64_t value : word) {
            if (value != 0U) return false;
        }
        return true;
    }
    [[nodiscard]] int count() const noexcept {
        int total = 0;
        for (std::uint64_t value : word) total += std::popcount(value);
        return total;
    }
    [[nodiscard]] int countAnd(const Bits& other) const noexcept {
        int total = 0;
        for (std::size_t i = 0; i < W; ++i) {
            total += std::popcount(word[i] & other.word[i]);
        }
        return total;
    }
    [[nodiscard]] Bits operator&(const Bits& other) const noexcept {
        Bits result{};
        for (std::size_t i = 0; i < W; ++i) result.word[i] = word[i] & other.word[i];
        return result;
    }
    [[nodiscard]] Bits operator|(const Bits& other) const noexcept {
        Bits result{};
        for (std::size_t i = 0; i < W; ++i) result.word[i] = word[i] | other.word[i];
        return result;
    }
    [[nodiscard]] Bits without(const Bits& other) const noexcept {
        Bits result{};
        for (std::size_t i = 0; i < W; ++i) result.word[i] = word[i] & ~other.word[i];
        return result;
    }
};

// One wide instance.  Held in the caller's frame, so the recursion allocates
// nothing: it carries three W-word sets by value and reads `rows` in place.
template <std::size_t W, class Emit>
struct WideSearch {
    static constexpr std::size_t capacity = 64U * W;

    std::size_t width = 0;      // |T|
    std::size_t coreSize = 0;   // |K| = |C| - |T|
    std::size_t minimumSize = 0;
    std::size_t cliqueSize = 0; // |C|, independent of |T|
    Bits<W> full{};
    std::array<Bits<W>, capacity> rows{};        // complement adjacency
    std::array<std::size_t, capacity> positions{};  // T index -> position in C
    std::array<std::uint64_t, capacity> words{};    // T index -> set-hash word
    std::vector<std::uint64_t>* masks = nullptr;    // one mask over C, reused
    Counters* counters = nullptr;
    const Emit* emit = nullptr;

    void visit(const Bits<W>& selected, Bits<W> p, Bits<W> x) {
        ++counters->nodes;
        // |P| serves the size bound and the pivot's early exits below.
        const int pCount = p.count();
        const std::size_t reachable =
            coreSize + static_cast<std::size_t>(selected.count() + pCount);
        if (reachable < minimumSize) {
            return;
        }
        if (pCount == 0) {
            // Nothing left to add: `selected` is maximal exactly when X is
            // empty; a non-empty X would only have left the pivot scan
            // without an extension.
            if (x.none()) {
                ++counters->leaves;
                emitLeaf(selected);
            }
            return;
        }
        // Pivot: the vertex of P u X with the most neighbours inside P, which
        // is what keeps the branching factor down; the first in index order
        // on ties.  The same two early exits as the single-word search, with
        // the same pivot and output: a vertex scoring |P| can only be in X
        // and extends everything below this call, so the call ends; with X
        // empty, the first vertex scoring |P|-1 cannot be beaten.
        const bool xEmpty = x.none();
        std::size_t pivot = 0;
        int best = -1;
        bool settled = false;
        const Bits<W> choices = p | x;
        for (std::size_t i = 0; i < W && !settled; ++i) {
            for (std::uint64_t rest = choices.word[i]; rest != 0U; rest &= rest - 1U) {
                const auto index =
                    i * 64U + static_cast<std::size_t>(std::countr_zero(rest));
                const int count = p.countAnd(rows[index]);
                if (count == pCount) {
                    return;
                }
                if (count > best) {
                    best = count;
                    pivot = index;
                    if (xEmpty && count == pCount - 1) {
                        settled = true;
                        break;
                    }
                }
            }
        }
        // The extension set is a snapshot: P and X are mutated as the loop
        // goes, exactly as in the single-word original.
        Bits<W> extensions = p.without(rows[pivot]);
        for (std::size_t i = 0; i < W; ++i) {
            while (extensions.word[i] != 0U) {
                const auto index =
                    i * 64U +
                    static_cast<std::size_t>(std::countr_zero(extensions.word[i]));
                extensions.word[i] &= extensions.word[i] - 1U;
                Bits<W> next = selected;
                next.set(index);
                visit(next, p & rows[index], x & rows[index]);
                p.reset(index);
                x.set(index);
            }
        }
    }

private:
    void emitLeaf(const Bits<W>& selected) {
        std::ranges::fill(*masks, ~std::uint64_t{0});
        if (cliqueSize % 64U != 0U) {
            masks->back() = (std::uint64_t{1} << (cliqueSize % 64U)) - 1U;
        }
        std::uint64_t omittedHash = 0U;
        const Bits<W> omitted = full.without(selected);
        for (std::size_t i = 0; i < W; ++i) {
            for (std::uint64_t rest = omitted.word[i]; rest != 0U; rest &= rest - 1U) {
                const auto index =
                    i * 64U + static_cast<std::size_t>(std::countr_zero(rest));
                const std::size_t position = positions[index];
                (*masks)[position / 64U] &= ~(std::uint64_t{1} << (position % 64U));
                omittedHash ^= words[index];
            }
        }
        (*emit)(std::span<const std::uint64_t>(*masks),
                coreSize + static_cast<std::size_t>(selected.count()), omittedHash);
    }
};

template <std::size_t W, class Edges, class Emit>
void runWide(
    std::span<const Vertex> clique,
    std::span<const Vertex> endpoints,
    const Edges& edges,
    std::size_t minimumSize,
    Counters& counters,
    const Emit& emit) {
    const std::size_t width = endpoints.size();
    auto state = std::make_unique<WideSearch<W, Emit>>();
    state->width = width;
    state->cliqueSize = clique.size();
    state->coreSize = clique.size() - width;
    state->minimumSize = minimumSize;
    state->full = Bits<W>::prefix(width);
    for (std::size_t i = 0; i < width; ++i) {
        state->rows[i] = state->full;
        state->rows[i].reset(i);
        state->positions[i] = static_cast<std::size_t>(
            std::lower_bound(clique.begin(), clique.end(), endpoints[i]) -
            clique.begin());
        state->words[i] = vertexSetWord(endpoints[i]);
    }
    for (const auto& edge : edges) {
        const auto u = static_cast<std::size_t>(
            std::lower_bound(endpoints.begin(), endpoints.end(), edge.u) -
            endpoints.begin());
        const auto v = static_cast<std::size_t>(
            std::lower_bound(endpoints.begin(), endpoints.end(), edge.v) -
            endpoints.begin());
        state->rows[u].reset(v);
        state->rows[v].reset(u);
    }
    std::vector<std::uint64_t> masks((clique.size() + 63U) / 64U);
    state->masks = &masks;
    state->counters = &counters;
    state->emit = &emit;
    state->visit(Bits<W>{}, state->full, Bits<W>{});
}

// `WideSearch` past four words: the same search -- the same pivot, early exits
// and size bound, so the same leaves in the same order -- with the word count
// fixed at run time instead of compile time.  Every set is a run of `words_`
// words.  P, X and the extension snapshot of each depth live in one arena
// (three sets per depth; a child is one deeper than its parent, so siblings
// reuse one slot), and the selected set is a single shared bit set kept in
// step with the recursion, so a node allocates nothing.
template <class Emit>
class DynamicSearch {
public:
    template <class Edges>
    DynamicSearch(
        std::span<const Vertex> clique,
        std::span<const Vertex> endpoints,
        const Edges& edges,
        std::size_t minimumSize,
        Counters& counters,
        const Emit& emit)
        : width_(endpoints.size()),
          words_((width_ + 63U) / 64U),
          coreSize_(clique.size() - width_),
          minimumSize_(minimumSize),
          cliqueSize_(clique.size()),
          counters_(counters),
          emit_(emit),
          rows_(width_ * words_, 0U),
          full_(words_, 0U),
          selected_(words_, 0U),
          // A node at depth d has d selected vertices and a child only while
          // P is non-empty, so the depths run 0..|T|.
          arena_((width_ + 1U) * 3U * words_, 0U),
          positions_(width_),
          setWords_(width_),
          masks_((clique.size() + 63U) / 64U) {
        for (std::size_t i = 0; i < width_; ++i) {
            full_[i / 64U] |= std::uint64_t{1} << (i % 64U);
        }
        for (std::size_t i = 0; i < width_; ++i) {
            std::uint64_t* rowWords = row(i);
            std::ranges::copy(full_, rowWords);
            rowWords[i / 64U] &= ~(std::uint64_t{1} << (i % 64U));
            positions_[i] = static_cast<std::size_t>(
                std::lower_bound(clique.begin(), clique.end(), endpoints[i]) -
                clique.begin());
            setWords_[i] = vertexSetWord(endpoints[i]);
        }
        for (const auto& edge : edges) {
            const auto u = static_cast<std::size_t>(
                std::lower_bound(endpoints.begin(), endpoints.end(), edge.u) -
                endpoints.begin());
            const auto v = static_cast<std::size_t>(
                std::lower_bound(endpoints.begin(), endpoints.end(), edge.v) -
                endpoints.begin());
            row(u)[v / 64U] &= ~(std::uint64_t{1} << (v % 64U));
            row(v)[u / 64U] &= ~(std::uint64_t{1} << (u % 64U));
        }
    }

    void run() {
        std::ranges::copy(full_, slot(0U, 0U));
        visit(0U);
    }

private:
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
        ++counters_.dynamicNodes;
        std::uint64_t* p = slot(depth, 0U);
        std::uint64_t* x = slot(depth, 1U);
        const int pCount = count(p);
        const std::size_t reachable =
            coreSize_ + selectedCount_ + static_cast<std::size_t>(pCount);
        if (reachable < minimumSize_) {
            return;
        }
        if (pCount == 0) {
            if (none(x)) {
                ++counters_.leaves;
                emitLeaf();
            }
            return;
        }
        // The pivot and its early exits, exactly as in `WideSearch::visit`.
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
        // The extension set is a snapshot: P and X are mutated as the loop
        // goes, exactly as in the single-word original.
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
        emit_(std::span<const std::uint64_t>(masks_), coreSize_ + selectedCount_,
              omittedHash);
    }

    std::size_t width_ = 0;        // |T|
    std::size_t words_ = 0;
    std::size_t coreSize_ = 0;     // |K| = |C| - |T|
    std::size_t minimumSize_ = 0;
    std::size_t cliqueSize_ = 0;   // |C|, independent of |T|
    std::size_t selectedCount_ = 0;
    Counters& counters_;
    const Emit& emit_;
    std::vector<std::uint64_t> rows_;      // complement adjacency, |T| rows
    std::vector<std::uint64_t> full_;
    std::vector<std::uint64_t> selected_;
    std::vector<std::uint64_t> arena_;
    std::vector<std::size_t> positions_;   // T index -> position in C
    std::vector<std::uint64_t> setWords_;  // T index -> set-hash word
    std::vector<std::uint64_t> masks_;     // one mask over C, reused
};

// The single-word search.  The pivot scan's early exits choose the same pivot,
// visit the same nodes and emit the same sets in the same order as the full
// scan; they only skip work whose result was already decided.
template <class Edges, class Emit>
void runNarrow(
    std::span<const Vertex> clique,
    std::span<const Vertex> endpoints,
    const Edges& edges,
    std::size_t minimumSize,
    Counters& counters,
    const Emit& emit) {
    const auto t = endpoints.size();
    const auto full = t == 64U ? ~std::uint64_t{0} : (std::uint64_t{1} << t) - 1U;
    std::array<std::uint64_t, 64U> rows{}, words{};
    std::array<std::size_t, 64U> positions{};
    for (std::size_t i = 0; i < t; ++i) {
        rows[i] = full & ~(std::uint64_t{1} << i);
        positions[i] = static_cast<std::size_t>(
            std::lower_bound(clique.begin(), clique.end(), endpoints[i]) -
            clique.begin());
        words[i] = vertexSetWord(endpoints[i]);
    }
    for (const auto& e : edges) {
        const auto u = static_cast<std::size_t>(
            std::lower_bound(endpoints.begin(), endpoints.end(), e.u) -
            endpoints.begin());
        const auto v = static_cast<std::size_t>(
            std::lower_bound(endpoints.begin(), endpoints.end(), e.v) -
            endpoints.begin());
        rows[u] &= ~(std::uint64_t{1} << v);
        rows[v] &= ~(std::uint64_t{1} << u);
    }
    const auto coreSize = clique.size() - t;
    std::vector<std::uint64_t> masks((clique.size() + 63U) / 64U);
    const auto visit = [&](auto&& self, std::uint64_t selected, std::uint64_t p,
                           std::uint64_t x) -> void {
        ++counters.nodes;
        // |P| serves the size bound and the pivot's early exits below.
        const int pCount = std::popcount(p);
        if (coreSize + static_cast<std::size_t>(std::popcount(selected) + pCount) <
            minimumSize) {
            return;
        }
        if (p == 0U) {
            // Nothing left to add: `selected` is maximal exactly when X is
            // empty; a non-empty X would only have left the pivot scan
            // without an extension.
            if (x != 0U) {
                return;
            }
            ++counters.leaves;
            std::ranges::fill(masks, ~std::uint64_t{0});
            if (clique.size() % 64U) {
                masks.back() = (std::uint64_t{1} << (clique.size() % 64U)) - 1U;
            }
            std::uint64_t omittedHash = 0;
            for (auto omitted = full & ~selected; omitted; omitted &= omitted - 1U) {
                const auto i = static_cast<std::size_t>(std::countr_zero(omitted));
                masks[positions[i] / 64U] &= ~(std::uint64_t{1} << (positions[i] % 64U));
                omittedHash ^= words[i];
            }
            emit(std::span<const std::uint64_t>(masks),
                 coreSize + static_cast<std::size_t>(std::popcount(selected)),
                 omittedHash);
            return;
        }
        // Pivot, first in bit order on ties.  A vertex scoring |P| can only be
        // in X (no row holds its own vertex) and extends every set below this
        // call, so the full scan would find no extension and the call ends;
        // with X empty, |P|-1 is the best score left, so the first vertex
        // reaching it is the pivot the full scan would keep.
        std::size_t pivot = 0;
        int best = -1;
        for (auto rest = p | x; rest; rest &= rest - 1U) {
            const auto i = static_cast<std::size_t>(std::countr_zero(rest));
            const int count = std::popcount(p & rows[i]);
            if (count == pCount) {
                return;
            }
            if (count > best) {
                best = count;
                pivot = i;
                if (x == 0U && count == pCount - 1) {
                    break;
                }
            }
        }
        auto extensions = p & ~rows[pivot];
        while (extensions) {
            const auto bit = extensions & (~extensions + 1U);
            extensions ^= bit;
            const auto i = static_cast<std::size_t>(std::countr_zero(bit));
            self(self, selected | bit, p & rows[i], x & rows[i]);
            p &= ~bit;
            x |= bit;
        }
    };
    visit(visit, 0U, full, 0U);
}

// `endpoints` is a sorted search domain containing both endpoints of every
// internal batch edge. An emitted mask is always indexed over the whole clique.
template <class Edges, class Emit>
void searchOnVertices(
    std::span<const Vertex> clique,
    std::span<const Vertex> endpoints,
    const Edges& edges,
    std::size_t minimumSize,
    Counters& counters,
    const Emit& emit) {
    ++counters.kernelCalls;
    const auto t = endpoints.size();
    if (t > maximumWidth) {
        ++counters.kernelDynamic;
        DynamicSearch<Emit> search(
            clique, endpoints, edges, minimumSize, counters, emit);
        search.run();
        return;
    }
    if (t <= 64U) {
        ++counters.kernelNarrow;
        runNarrow(clique, endpoints, edges, minimumSize, counters, emit);
        return;
    }
    ++counters.kernelWide;
    const std::size_t neededWords = (t + 63U) / 64U;
    if (neededWords == 2U) {
        runWide<2U>(clique, endpoints, edges, minimumSize, counters, emit);
    } else if (neededWords == 3U) {
        runWide<3U>(clique, endpoints, edges, minimumSize, counters, emit);
    } else {
        runWide<4U>(clique, endpoints, edges, minimumSize, counters, emit);
    }
}

template <class Edges, class Emit>
void searchKernel(
    std::span<const Vertex> clique,
    const Edges& edges,
    std::size_t minimumSize,
    Counters& counters,
    const Emit& emit) {
    static thread_local std::vector<Vertex> endpoints;
    endpoints.clear();
    for (const auto& edge : edges) {
        endpoints.push_back(edge.u);
        endpoints.push_back(edge.v);
    }
    std::ranges::sort(endpoints);
    endpoints.erase(std::unique(endpoints.begin(), endpoints.end()), endpoints.end());
    searchOnVertices(clique, endpoints, edges, minimumSize, counters, emit);
}

} // namespace immce::wideresidual
