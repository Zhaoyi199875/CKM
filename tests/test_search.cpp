#include "maintenance_engine.hpp"
#include "maximal_clique_index.hpp"
#include "new_clique_search.hpp"
#include "star_clique_search.hpp"
#include "subsumed_clique_search.hpp"
#include "vertex_clique_search.hpp"
#include "wide_residual.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace immce {
namespace {

using CliqueSet = std::set<Clique, CliqueLexicographicOrder>;

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void require(bool condition, const std::string& message) {
    if (!condition) {
        fail(message);
    }
}

[[nodiscard]] std::vector<Edge> completeCatalog(int vertexCount) {
    std::vector<Edge> result;
    for (Vertex u = 0; u < vertexCount; ++u) {
        for (Vertex v = u + 1; v < vertexCount; ++v) {
            result.push_back({u, v});
        }
    }
    return result;
}

[[nodiscard]] CliqueSet asSet(
    const std::vector<Clique>& cliques,
    const std::string& label) {
    CliqueSet result;
    for (Clique clique : cliques) {
        clique = canonicalClique(std::move(clique));
        if (!result.insert(std::move(clique)).second) {
            fail(label + " emitted a duplicate clique");
        }
    }
    return result;
}

[[nodiscard]] CliqueSet bruteForceMaximalCliques(
    int vertexCount,
    const std::vector<Edge>& catalog,
    const std::vector<unsigned char>& active) {
    std::vector<unsigned char> activeVertex(
        static_cast<std::size_t>(vertexCount), 0U);
    for (std::size_t edgeId = 0; edgeId < catalog.size(); ++edgeId) {
        if (active[edgeId] != 0U) {
            activeVertex[static_cast<std::size_t>(catalog[edgeId].u)] = 1U;
            activeVertex[static_cast<std::size_t>(catalog[edgeId].v)] = 1U;
        }
    }

    CliqueSet result;
    const std::uint64_t limit = std::uint64_t{1} << vertexCount;
    for (std::uint64_t mask = 1U; mask < limit; ++mask) {
        bool usesOnlyActiveVertices = true;
        for (Vertex vertex = 0; vertex < vertexCount; ++vertex) {
            if ((mask & (std::uint64_t{1} << vertex)) != 0U &&
                activeVertex[static_cast<std::size_t>(vertex)] == 0U) {
                usesOnlyActiveVertices = false;
                break;
            }
        }
        if (!usesOnlyActiveVertices) {
            continue;
        }

        bool clique = true;
        for (std::size_t edgeId = 0; edgeId < catalog.size(); ++edgeId) {
            const Edge& edge = catalog[edgeId];
            const bool both = (mask & (std::uint64_t{1} << edge.u)) != 0U &&
                (mask & (std::uint64_t{1} << edge.v)) != 0U;
            if (both && active[edgeId] == 0U) {
                clique = false;
                break;
            }
        }
        if (!clique) {
            continue;
        }

        bool maximal = true;
        for (Vertex extension = 0; extension < vertexCount && maximal; ++extension) {
            if (activeVertex[static_cast<std::size_t>(extension)] == 0U ||
                (mask & (std::uint64_t{1} << extension)) != 0U) {
                continue;
            }
            bool adjacentToAll = true;
            for (Vertex member = 0; member < vertexCount; ++member) {
                if ((mask & (std::uint64_t{1} << member)) == 0U) {
                    continue;
                }
                const Edge wanted = normalizedEdge(extension, member);
                const auto found = std::lower_bound(
                    catalog.begin(), catalog.end(), wanted);
                if (found == catalog.end() || *found != wanted ||
                    active[static_cast<std::size_t>(found - catalog.begin())] == 0U) {
                    adjacentToAll = false;
                    break;
                }
            }
            if (adjacentToAll) {
                maximal = false;
            }
        }
        if (!maximal) {
            continue;
        }

        Clique value;
        for (Vertex vertex = 0; vertex < vertexCount; ++vertex) {
            if ((mask & (std::uint64_t{1} << vertex)) != 0U) {
                value.push_back(vertex);
            }
        }
        result.insert(std::move(value));
    }
    return result;
}

void requireIndexEquals(
    const MaximalCliqueIndex& index,
    const CliqueSet& expected,
    const std::string& label) {
    require(index.size() == expected.size(), label + " index size mismatch");
    for (const Clique& clique : expected) {
        require(index.contains(clique), label + " index is missing a clique");
    }
}

// One row per (framework, pipeline) pair the binary can be asked for.  A new
// subsumed framework gets a row here and is then held to brute force over
// every state below, in both pipelines and behind both producers.
struct Variant {
    const char* name;
    SubsumedCliqueAlgorithm algorithm;
    bool streaming = false;
    // Probes are held back to the end of a batch only once the index holds
    // 4,096 live cliques, which no state here reaches.  `forceDeferral` drops
    // that floor to zero so the buffered path -- including the flush that
    // erases while later probes are still queued -- is held to brute force
    // like the inline one.
    bool forceDeferral = false;
    // Maintain the live clique-size histogram when the framework reads it;
    // false keeps it off, which switches the size gate off.
    bool sizeHistogram = true;
};

[[nodiscard]] std::span<const Variant> variants() {
    using Algorithm = SubsumedCliqueAlgorithm;
    constexpr auto optimized = Algorithm::Optimized;
    constexpr auto bitmap = Algorithm::OptimizedBitmap;
    constexpr auto rk = Algorithm::ResidualKernel;
    static constexpr std::array table{
        Variant{.name = "optimized", .algorithm = optimized},
        Variant{.name = "optimized/stream", .algorithm = optimized, .streaming = true},
        Variant{.name = "optimized/defer", .algorithm = optimized,
                .forceDeferral = true},
        Variant{.name = "optimized/defer/stream", .algorithm = optimized,
                .streaming = true, .forceDeferral = true},
        Variant{.name = "optimized-bitmap", .algorithm = bitmap},
        Variant{.name = "optimized-bitmap/stream", .algorithm = bitmap,
                .streaming = true},
        Variant{.name = "optimized-bitmap/defer", .algorithm = bitmap,
                .forceDeferral = true},
        Variant{.name = "optimized-bitmap/defer/stream", .algorithm = bitmap,
                .streaming = true, .forceDeferral = true},
        Variant{.name = "residual-kernel", .algorithm = rk},
        Variant{.name = "residual-kernel/stream", .algorithm = rk, .streaming = true},
        Variant{.name = "residual-kernel/defer", .algorithm = rk,
                .forceDeferral = true},
        Variant{.name = "residual-kernel/defer/stream", .algorithm = rk,
                .streaming = true, .forceDeferral = true},
        // The size gate switched off: the histogram is not kept, so every
        // test of it answers "possible".
        Variant{.name = "residual-kernel/no-size-histogram", .algorithm = rk,
                .sizeHistogram = false},
        Variant{.name = "residual-kernel/no-size-histogram/stream", .algorithm = rk,
                .streaming = true, .sizeHistogram = false},
        Variant{.name = "residual-kernel/no-size-histogram/defer", .algorithm = rk,
                .forceDeferral = true, .sizeHistogram = false},
        Variant{.name = "residual-kernel/no-size-histogram/defer/stream",
                .algorithm = rk, .streaming = true, .forceDeferral = true,
                .sizeHistogram = false},
        Variant{.name = "full-bk", .algorithm = Algorithm::FullBK,
                .sizeHistogram = false},
        Variant{.name = "full-bk/defer/stream", .algorithm = Algorithm::FullBK,
                .streaming = true, .forceDeferral = true, .sizeHistogram = false},
        Variant{.name = "kernel-bk", .algorithm = Algorithm::KernelBK,
                .sizeHistogram = false},
        Variant{.name = "kernel-bk/defer/stream", .algorithm = Algorithm::KernelBK,
                .streaming = true, .forceDeferral = true, .sizeHistogram = false},
    };
    return table;
}

// The vertex baseline's other configurations: its branchings besides the
// pivot, the ablation of 2026-10-21 (VertexBranching), and every branching
// from the cover vertices alone, the ablation of 2026-11-01 (VertexSeeds).
[[nodiscard]] std::vector<VertexSearchConfiguration> vertexAblationMethods() {
    std::vector<VertexSearchConfiguration> table;
    for (const VertexBranching branching :
         {VertexBranching::NewEdges, VertexBranching::RuleRoot,
          VertexBranching::RuleUnbound}) {
        table.push_back(
            VertexSearchConfiguration{branching, VertexSeeds::Endpoints});
    }
    for (const VertexBranching branching :
         {VertexBranching::Pivot, VertexBranching::NewEdges,
          VertexBranching::RuleRoot, VertexBranching::RuleUnbound}) {
        table.push_back(VertexSearchConfiguration{branching, VertexSeeds::Cover});
    }
    return table;
}

[[nodiscard]] std::string vertexAblationLabel(
    const VertexSearchConfiguration& method) {
    return std::string("vertex/") + vertexBranchingName(method.branching) +
        (method.seeds == VertexSeeds::Cover ? "/cover" : "");
}

// The star method with one of its two partitions fixed, with the root filter
// on and off: the branching ablation of 2026-10-31 (StarBranching). Each one
// has to find exactly what the method finds.
[[nodiscard]] std::vector<StarSearchConfiguration> starAblationMethods() {
    std::vector<StarSearchConfiguration> table;
    for (const auto branching : {StarBranching::Anchor, StarBranching::Pivot}) {
        for (const bool filter : {true, false}) {
            StarSearchConfiguration configuration;
            configuration.branching = branching;
            configuration.rootFilter = filter;
            table.push_back(configuration);
        }
    }
    return table;
}

[[nodiscard]] std::string starAblationLabel(const StarSearchConfiguration& method) {
    return std::string("star/") + starBranchingName(method.branching) +
        (method.rootFilter ? "" : "/no-root-filter");
}

// The two star methods. Re-indexing never opens on the handful of vertices
// the exhaustive states use -- a graph of one word cannot get narrower -- so
// the wide-domain tests below are where it is held to the reference.
[[nodiscard]] std::vector<StarSearchConfiguration> starMethods() {
    std::vector<StarSearchConfiguration> table;
    for (const auto scope : {StarBranchScope::Unbound, StarBranchScope::Root}) {
        StarSearchConfiguration configuration;
        configuration.branchScope = scope;
        table.push_back(configuration);
    }
    return table;
}

// The bookkeeping every star run has to balance: each batch edge placed in
// exactly one star, the centres drawn from the cover, and every star searched
// by exactly one of the three shapes.
void requireStarAccounts(const StarSearchStats& stats, const std::string& label) {
    require(stats.starEdgeOccurrences == stats.batchEdges &&
                stats.singletonStarEdges + stats.perPartnerStarEdges +
                        stats.mergedStarEdges ==
                    stats.batchEdges,
            label + " did not place every batch edge in exactly one star");
    require(stats.coverVertices >= stats.stars,
            label + " used more centres than the cover holds");
    require(stats.singletonStars + stats.branchRootEdgeStars +
                    stats.branchRootPivotStars ==
                stats.stars,
            label + " did not search every star exactly once");
}

[[nodiscard]] std::vector<Clique> enumerateStar(
    const IncrementalGraph& graph,
    std::span<const EdgeId> batch,
    const StarSearchConfiguration& configuration) {
    StarCliqueEnumerator enumerator(graph, configuration);
    StarSearchStats stats;
    std::vector<Clique> cliques = enumerator.enumerateBatch(batch, stats);
    requireStarAccounts(
        stats,
        std::string("star/") + starBranchScopeName(configuration.branchScope));
    return cliques;
}

// Maximal cliques by Bron-Kerbosch over bitset rows, so the structured cases
// below -- a clique of a hundred vertices, sixty-five changed endpoints --
// have an independent reference that the exponential enumerator cannot give.
// Isolated vertices are left out, matching `bruteForceMaximalCliques`.
class BitsetGraph {
public:
    BitsetGraph(int vertexCount, const std::vector<Edge>& edges)
        : vertexCount_(static_cast<std::size_t>(vertexCount)),
          words_((vertexCount_ + 63U) / 64U),
          rows_(vertexCount_ * words_, 0U) {
        for (const Edge& edge : edges) {
            set(static_cast<std::size_t>(edge.u),
                static_cast<std::size_t>(edge.v));
            set(static_cast<std::size_t>(edge.v),
                static_cast<std::size_t>(edge.u));
        }
    }

    [[nodiscard]] std::size_t vertexCount() const noexcept {
        return vertexCount_;
    }
    [[nodiscard]] std::size_t words() const noexcept { return words_; }
    [[nodiscard]] const std::uint64_t* row(std::size_t vertex) const noexcept {
        return rows_.data() + vertex * words_;
    }

private:
    void set(std::size_t from, std::size_t to) {
        rows_[from * words_ + (to / 64U)] |= std::uint64_t{1} << (to % 64U);
    }

    std::size_t vertexCount_ = 0U;
    std::size_t words_ = 0U;
    std::vector<std::uint64_t> rows_;
};

void bronKerbosch(
    const BitsetGraph& graph,
    Clique& current,
    std::vector<std::uint64_t> candidates,
    std::vector<std::uint64_t> excluded,
    CliqueSet& result) {
    const std::size_t words = graph.words();
    const auto empty = [words](const std::vector<std::uint64_t>& set) {
        for (std::size_t word = 0; word < words; ++word) {
            if (set[word] != 0U) {
                return false;
            }
        }
        return true;
    };
    if (empty(candidates) && empty(excluded)) {
        if (!current.empty()) {
            Clique value = current;
            std::ranges::sort(value);
            result.insert(std::move(value));
        }
        return;
    }
    // Pivot: the candidate or excluded vertex with the most neighbours among
    // the candidates.
    std::size_t pivot = graph.vertexCount();
    std::size_t best = 0U;
    for (std::size_t vertex = 0; vertex < graph.vertexCount(); ++vertex) {
        const std::uint64_t bit = std::uint64_t{1} << (vertex % 64U);
        const std::size_t word = vertex / 64U;
        if ((candidates[word] & bit) == 0U && (excluded[word] & bit) == 0U) {
            continue;
        }
        std::size_t count = 0U;
        for (std::size_t index = 0; index < words; ++index) {
            count += static_cast<std::size_t>(
                std::popcount(candidates[index] & graph.row(vertex)[index]));
        }
        if (pivot == graph.vertexCount() || count > best) {
            pivot = vertex;
            best = count;
        }
    }
    std::vector<std::uint64_t> branch(words, 0U);
    for (std::size_t word = 0; word < words; ++word) {
        branch[word] = candidates[word] & ~graph.row(pivot)[word];
    }
    for (std::size_t vertex = 0; vertex < graph.vertexCount(); ++vertex) {
        const std::uint64_t bit = std::uint64_t{1} << (vertex % 64U);
        const std::size_t word = vertex / 64U;
        if ((branch[word] & bit) == 0U) {
            continue;
        }
        std::vector<std::uint64_t> nextCandidates(words, 0U);
        std::vector<std::uint64_t> nextExcluded(words, 0U);
        for (std::size_t index = 0; index < words; ++index) {
            nextCandidates[index] =
                candidates[index] & graph.row(vertex)[index];
            nextExcluded[index] = excluded[index] & graph.row(vertex)[index];
        }
        current.push_back(static_cast<Vertex>(vertex));
        bronKerbosch(
            graph, current, std::move(nextCandidates),
            std::move(nextExcluded), result);
        current.pop_back();
        candidates[word] &= ~bit;
        excluded[word] |= bit;
    }
}

[[nodiscard]] CliqueSet referenceMaximalCliques(
    int vertexCount,
    const std::vector<Edge>& edges) {
    const BitsetGraph graph(vertexCount, edges);
    std::vector<std::uint64_t> candidates(graph.words(), 0U);
    for (const Edge& edge : edges) {
        candidates[static_cast<std::size_t>(edge.u) / 64U] |=
            std::uint64_t{1} << (static_cast<std::size_t>(edge.u) % 64U);
        candidates[static_cast<std::size_t>(edge.v) / 64U] |=
            std::uint64_t{1} << (static_cast<std::size_t>(edge.v) % 64U);
    }
    CliqueSet result;
    Clique current;
    bronKerbosch(
        graph, current, std::move(candidates),
        std::vector<std::uint64_t>(graph.words(), 0U), result);
    return result;
}

// Runs one batch the way the binary would, and leaves the index holding the
// updated maximal-clique set: subsumed cliques erased, new ones inserted.
[[nodiscard]] std::vector<Clique> runVariant(
    SubsumedCliqueEnumerator& enumerator,
    const Variant& variant,
    std::span<const EdgeId> batch,
    const std::vector<Clique>& newCliques,
    MaximalCliqueIndex& index,
    SubsumedCliqueSearchStats& stats) {
    if (!variant.streaming) {
        std::vector<Clique> deleted = enumerator.enumerateBatch(
            batch, newCliques, index, stats, variant.algorithm);
        const std::span<const std::uint64_t> hashes =
            enumerator.newCliqueHashes();
        require(hashes.size() == newCliques.size(),
                "the subsumed search did not hash every new clique");
        for (std::size_t position = 0; position < newCliques.size();
             ++position) {
            require(index.insertHashed(newCliques[position], hashes[position]),
                    "new clique was already present in old index");
        }
        return deleted;
    }

    CliqueBatchOutput output(true);
    enumerator.prepareBatch(batch, variant.algorithm, stats);
    for (std::size_t position = 0; position < newCliques.size(); ++position) {
        const Clique& clique = newCliques[position];
        enumerator.consumeNewClique(
            clique, index, variant.algorithm, stats, output);
        require(enumerator.newCliqueHashes().size() == position + 1U,
                "the subsumed search did not hash the new clique");
        require(index.insertHashed(clique, enumerator.newCliqueHashes().back()),
                "new clique was already present in old index");
    }
    enumerator.finishBatch(index, stats, output);
    return std::move(output).takeCliques();
}

void verifyState(int vertexCount, const std::vector<unsigned char>& state) {
    const std::vector<Edge> catalog = completeCatalog(vertexCount);
    std::vector<EdgeId> oldEdges;
    std::vector<EdgeId> batch;
    std::vector<unsigned char> before(catalog.size(), 0U);
    std::vector<unsigned char> after(catalog.size(), 0U);
    for (std::size_t edgeId = 0; edgeId < state.size(); ++edgeId) {
        if (state[edgeId] == 1U) {
            oldEdges.push_back(static_cast<EdgeId>(edgeId));
            before[edgeId] = 1U;
            after[edgeId] = 1U;
        } else if (state[edgeId] == 2U) {
            batch.push_back(static_cast<EdgeId>(edgeId));
            after[edgeId] = 1U;
        }
    }
    if (batch.empty()) {
        return;
    }

    const CliqueSet beforeCliques =
        bruteForceMaximalCliques(vertexCount, catalog, before);
    const CliqueSet afterCliques =
        bruteForceMaximalCliques(vertexCount, catalog, after);
    CliqueSet expectedNew = afterCliques;
    CliqueSet expectedDeleted = beforeCliques;
    for (const Clique& clique : beforeCliques) {
        expectedNew.erase(clique);
    }
    for (const Clique& clique : afterCliques) {
        expectedDeleted.erase(clique);
    }

    IncrementalGraph graph(vertexCount, catalog);
    graph.applyBatch(oldEdges);
    graph.applyBatch(batch);
    NewCliqueEnumerator newEnumerator(graph);
    NewCliqueSearchStats newStats;
    const std::vector<Clique> newCliques =
        newEnumerator.enumerateBatch(batch, newStats);
    require(asSet(newCliques, "new search") == expectedNew,
            "new clique search disagrees with brute force");

    for (const StarSearchConfiguration& method : starMethods()) {
        const std::string label =
            std::string("star/") + starBranchScopeName(method.branchScope);
        require(asSet(enumerateStar(graph, batch, method), label) == expectedNew,
                label + " disagrees with brute force");
        // The ablation: the same method on the centres' whole neighbourhoods.
        StarSearchConfiguration unfiltered = method;
        unfiltered.rootFilter = false;
        require(asSet(enumerateStar(graph, batch, unfiltered),
                      label + "/no-root-filter") == expectedNew,
                label + "/no-root-filter disagrees with brute force");
    }
    for (const StarSearchConfiguration& method : starAblationMethods()) {
        const std::string label = starAblationLabel(method);
        require(asSet(enumerateStar(graph, batch, method), label) == expectedNew,
                label + " disagrees with brute force");
    }
    VertexCliqueEnumerator vertexEnumerator(graph);
    VertexSearchStats vertexStats;
    const std::vector<Clique> vertexCliques =
        vertexEnumerator.enumerateBatch(batch, vertexStats);
    require(asSet(vertexCliques, "vertex") == expectedNew,
            "vertex disagrees with brute force");
    for (const VertexSearchConfiguration& method : vertexAblationMethods()) {
        const std::string label = vertexAblationLabel(method);
        VertexCliqueEnumerator ablated(graph, method);
        VertexSearchStats ablatedStats;
        require(asSet(ablated.enumerateBatch(batch, ablatedStats), label) ==
                    expectedNew,
                label + " disagrees with brute force");
    }

    // A batch can hand its new cliques to the subsumed search in any order:
    // no new maximal clique of a batch is a proper subset of another, so the
    // producer that walked the batch must not change the result. Running the
    // whole variant table against every producer's order is what holds the
    // two halves of this project to that.
    const std::vector<Clique> starCliques =
        enumerateStar(graph, batch, StarSearchConfiguration{});
    struct Producer {
        const char* name;
        const std::vector<Clique>* cliques;
    };
    const std::array<Producer, 3U> producers{
        Producer{"edge", &newCliques},
        Producer{"star", &starCliques},
        Producer{"vertex", &vertexCliques}};
    for (const auto& [producer, producedCliques] : producers) {
        for (const Variant& variant : variants()) {
            MaximalCliqueIndex index(
                usesSubsumedSizeGate(variant.algorithm) && variant.sizeHistogram);
            for (const Clique& clique : beforeCliques) {
                require(index.insert(clique),
                        "old clique index rejected a unique clique");
            }
            SubsumedCliqueEnumerator subsumedEnumerator(graph);
            if (variant.forceDeferral) {
                subsumedEnumerator.setDeferralMinimumIndexSize(0U);
            }
            SubsumedCliqueSearchStats subsumedStats;
            const std::string label =
                std::string(producer) + "/" + variant.name;
            const std::vector<Clique> deleted = runVariant(
                subsumedEnumerator, variant, batch, *producedCliques, index,
                subsumedStats);
            require(asSet(deleted, variant.name) == expectedDeleted,
                    label + " subsumed clique search disagrees with brute force");
            requireIndexEquals(index, afterCliques, label + " updated");
        }
    }
}

[[nodiscard]] CliqueSet setDifference(const CliqueSet& lhs, const CliqueSet& rhs) {
    CliqueSet result;
    std::ranges::set_difference(
        lhs, rhs, std::inserter(result, result.end()), CliqueLexicographicOrder{});
    return result;
}

void requireSameGraph(
    const IncrementalGraph& graph,
    const IncrementalGraph& reference,
    const std::string& label) {
    require(graph.activeEdgeCount() == reference.activeEdgeCount() &&
                std::ranges::equal(graph.activeVertices(),
                                   reference.activeVertices()),
            label + " left a different set of edges or active vertices");
    for (Vertex vertex = 0; vertex < graph.vertexCount(); ++vertex) {
        require(std::ranges::equal(
                    graph.neighbors(vertex), reference.neighbors(vertex),
                    [](const AdjacencyEntry& lhs, const AdjacencyEntry& rhs) {
                        return lhs.neighbor == rhs.neighbor &&
                            lhs.edgeId == rhs.edgeId;
                    }),
                label + " left a different adjacency");
    }
}

// The same states read as deletions: 1 marks an edge the batch keeps, 2 one
// it deletes.  The cliques the batch destroys are every producer's answer for
// the deleted edges on the graph that still holds them; then the edges go,
// and each subsumed algorithm, run backwards, must recover exactly the
// cliques the deletion creates and leave the index holding the maximal
// cliques after it.
void verifyDeletionState(int vertexCount, const std::vector<unsigned char>& state) {
    const std::vector<Edge> catalog = completeCatalog(vertexCount);
    std::vector<EdgeId> keptEdges;
    std::vector<EdgeId> deleted;
    std::vector<unsigned char> before(catalog.size(), 0U);
    std::vector<unsigned char> after(catalog.size(), 0U);
    for (std::size_t edgeId = 0; edgeId < state.size(); ++edgeId) {
        if (state[edgeId] == 1U) {
            keptEdges.push_back(static_cast<EdgeId>(edgeId));
            before[edgeId] = 1U;
            after[edgeId] = 1U;
        } else if (state[edgeId] == 2U) {
            deleted.push_back(static_cast<EdgeId>(edgeId));
            before[edgeId] = 1U;
        }
    }
    if (deleted.empty()) {
        return;
    }
    const CliqueSet beforeCliques =
        bruteForceMaximalCliques(vertexCount, catalog, before);
    const CliqueSet afterCliques =
        bruteForceMaximalCliques(vertexCount, catalog, after);
    const CliqueSet expectedDestroyed = setDifference(beforeCliques, afterCliques);
    const CliqueSet expectedCreated = setDifference(afterCliques, beforeCliques);

    IncrementalGraph graph(vertexCount, catalog);
    graph.applyBatch(keptEdges);
    graph.applyBatch(deleted);
    NewCliqueEnumerator edge(graph);
    NewCliqueSearchStats edgeStats;
    const std::vector<Clique> destroyed = edge.enumerateBatch(deleted, edgeStats);
    require(asSet(destroyed, "edge/delete") == expectedDestroyed,
            "edge disagrees with brute force on the destroyed cliques");
    for (const StarSearchConfiguration& method : starMethods()) {
        const std::string label =
            std::string("star/") + starBranchScopeName(method.branchScope);
        require(asSet(enumerateStar(graph, deleted, method), label) ==
                    expectedDestroyed,
                label + " disagrees with brute force on the destroyed cliques");
    }
    for (const StarSearchConfiguration& method : starAblationMethods()) {
        const std::string label = starAblationLabel(method);
        require(asSet(enumerateStar(graph, deleted, method), label) ==
                    expectedDestroyed,
                label + " disagrees with brute force on the destroyed cliques");
    }
    VertexCliqueEnumerator vertex(graph);
    VertexSearchStats vertexStats;
    require(asSet(vertex.enumerateBatch(deleted, vertexStats), "vertex/delete") ==
                expectedDestroyed,
            "vertex disagrees with brute force on the destroyed cliques");
    for (const VertexSearchConfiguration& method : vertexAblationMethods()) {
        const std::string label = vertexAblationLabel(method) + "/delete";
        VertexCliqueEnumerator ablated(graph, method);
        VertexSearchStats ablatedStats;
        require(asSet(ablated.enumerateBatch(deleted, ablatedStats), label) ==
                    expectedDestroyed,
                label + " disagrees with brute force on the destroyed cliques");
    }

    graph.removeBatch(deleted);
    IncrementalGraph reference(vertexCount, catalog);
    reference.applyBatch(keptEdges);
    requireSameGraph(graph, reference, "removeBatch");

    for (const auto algorithm :
         {SubsumedCliqueAlgorithm::Optimized,
          SubsumedCliqueAlgorithm::OptimizedBitmap,
          SubsumedCliqueAlgorithm::ResidualKernel}) {
        const std::string label =
            std::string("recover/") + subsumedCliqueAlgorithmName(algorithm);
        MaximalCliqueIndex index;
        for (const Clique& clique : beforeCliques) {
            require(index.insert(clique), "old clique index rejected a unique clique");
        }
        SubsumedCliqueEnumerator enumerator(graph);
        SubsumedCliqueSearchStats stats;
        const std::vector<Clique> created =
            enumerator.recoverBatch(deleted, destroyed, index, stats, algorithm);
        require(asSet(created, label) == expectedCreated,
                label + " disagrees with brute force on the created cliques");
        for (const Clique& clique : destroyed) {
            require(index.eraseCanonical(canonicalClique(clique)),
                    label + " index lost a destroyed clique");
        }
        requireIndexEquals(index, afterCliques, label + " after the deletion");
    }
}

// Whole streams there and back: random graphs cut into random batches,
// inserted with star and residual-kernel, then deleted again from the last
// batch on until the graph is empty.  Every deletion is held to brute force
// and the index to the maximal cliques after it.  The producer and the
// algorithm rotate with the trial.
void reverseReplayStreams() {
    std::mt19937 generator(20261010U);
    for (int trial = 0; trial < 300; ++trial) {
        const int vertexCount = std::uniform_int_distribution<int>(5, 8)(generator);
        const std::vector<Edge> catalog = completeCatalog(vertexCount);
        std::bernoulli_distribution keep(
            std::uniform_real_distribution<double>(0.3, 0.95)(generator));
        std::vector<EdgeId> stream;
        for (std::size_t id = 0; id < catalog.size(); ++id) {
            if (keep(generator)) {
                stream.push_back(static_cast<EdgeId>(id));
            }
        }
        std::ranges::shuffle(stream, generator);
        std::vector<std::vector<EdgeId>> batches;
        for (std::size_t begin = 0; begin < stream.size();) {
            const std::size_t size = std::min(
                stream.size() - begin,
                std::uniform_int_distribution<std::size_t>(1U, 6U)(generator));
            batches.emplace_back(
                stream.begin() + static_cast<std::ptrdiff_t>(begin),
                stream.begin() + static_cast<std::ptrdiff_t>(begin + size));
            begin += size;
        }

        IncrementalGraph graph(vertexCount, catalog);
        MaximalCliqueIndex index;
        SubsumedCliqueEnumerator subsumed(graph);
        StarCliqueEnumerator builder(graph);
        std::vector<unsigned char> active(catalog.size(), 0U);
        for (const std::vector<EdgeId>& batch : batches) {
            graph.applyBatch(batch);
            StarSearchStats starStats;
            std::vector<Clique> fresh = builder.enumerateBatch(batch, starStats);
            for (Clique& clique : fresh) {
                clique = canonicalClique(std::move(clique));
            }
            SubsumedCliqueSearchStats stats;
            static_cast<void>(subsumed.enumerateBatch(
                batch, fresh, index, stats, SubsumedCliqueAlgorithm::ResidualKernel));
            for (std::size_t position = 0; position < fresh.size(); ++position) {
                require(index.insertHashed(
                            fresh[position], subsumed.newCliqueHashes()[position]),
                        "the build inserted a clique twice");
            }
            for (const EdgeId id : batch) {
                active[static_cast<std::size_t>(id)] = 1U;
            }
        }

        const auto algorithm = std::array{
            SubsumedCliqueAlgorithm::Optimized,
            SubsumedCliqueAlgorithm::OptimizedBitmap,
            SubsumedCliqueAlgorithm::ResidualKernel}[
            static_cast<std::size_t>(trial % 3)];
        const int producer = (trial / 3) % 3;
        const std::string label = "reverse replay trial " + std::to_string(trial) +
            "/" + subsumedCliqueAlgorithmName(algorithm) + "/producer " +
            std::to_string(producer);
        NewCliqueEnumerator edge(graph);
        VertexCliqueEnumerator vertex(graph);
        StarCliqueEnumerator star(graph);
        for (auto batch = batches.rbegin(); batch != batches.rend(); ++batch) {
            const CliqueSet beforeCliques =
                bruteForceMaximalCliques(vertexCount, catalog, active);
            for (const EdgeId id : *batch) {
                active[static_cast<std::size_t>(id)] = 0U;
            }
            const CliqueSet afterCliques =
                bruteForceMaximalCliques(vertexCount, catalog, active);

            std::vector<Clique> destroyed;
            if (producer == 0) {
                StarSearchStats stats;
                destroyed = star.enumerateBatch(*batch, stats);
            } else if (producer == 1) {
                NewCliqueSearchStats stats;
                destroyed = edge.enumerateBatch(*batch, stats);
            } else {
                VertexSearchStats stats;
                destroyed = vertex.enumerateBatch(*batch, stats);
            }
            require(asSet(destroyed, label) == setDifference(beforeCliques, afterCliques),
                    label + ": destroyed cliques disagree with brute force");
            graph.removeBatch(*batch);
            SubsumedCliqueSearchStats stats;
            const std::vector<Clique> created =
                subsumed.recoverBatch(*batch, destroyed, index, stats, algorithm);
            require(asSet(created, label) == setDifference(afterCliques, beforeCliques),
                    label + ": created cliques disagree with brute force");
            for (const Clique& clique : destroyed) {
                require(index.eraseCanonical(canonicalClique(clique)),
                        label + ": the index lost a destroyed clique");
            }
            requireIndexEquals(index, afterCliques, label);
        }
        require(index.empty() && graph.activeEdgeCount() == 0 &&
                    graph.activeVertices().empty(),
                label + ": deleting every batch did not empty the graph");
    }
}

// The maintenance engine's configurations for batches that delete and insert
// at once: every producer (both star scopes), every subsumed algorithm (the
// residual-kernel size gate on and off), both strategies, and SOMEi (both
// update scopes).
struct EngineCase {
    std::string label;
    MaintenanceConfiguration configuration;
};

[[nodiscard]] const std::vector<EngineCase>& mixedEngineCases() {
    static const std::vector<EngineCase> cases = [] {
        struct Producer {
            const char* name;
            bool star;
            bool edge;
            bool vertex;
            StarBranchScope scope;
        };
        const std::array producers{
            Producer{"star-unbound", true, false, false, StarBranchScope::Unbound},
            Producer{"star-root", true, false, false, StarBranchScope::Root},
            Producer{"edge", false, true, false, StarBranchScope::Unbound},
            Producer{"vertex", false, false, true, StarBranchScope::Unbound}};
        struct Algorithm {
            const char* name;
            SubsumedCliqueAlgorithm algorithm;
            bool sizeHistogram;
        };
        const std::array algorithms{
            Algorithm{"opt", SubsumedCliqueAlgorithm::Optimized, false},
            Algorithm{"optbm", SubsumedCliqueAlgorithm::OptimizedBitmap, false},
            Algorithm{"rk", SubsumedCliqueAlgorithm::ResidualKernel, true},
            Algorithm{"rk-nogate", SubsumedCliqueAlgorithm::ResidualKernel, false}};
        std::vector<EngineCase> table;
        for (const MixedStrategy strategy :
             {MixedStrategy::SinglePass, MixedStrategy::TwoPhase}) {
            for (const Producer& producer : producers) {
                for (const Algorithm& algorithm : algorithms) {
                    MaintenanceConfiguration configuration;
                    configuration.runStar = producer.star;
                    configuration.runEdge = producer.edge;
                    configuration.runVertex = producer.vertex;
                    configuration.fullMaintenance = true;
                    configuration.builds = true;
                    configuration.starBranchScope = producer.scope;
                    configuration.subsumedAlgorithm = algorithm.algorithm;
                    configuration.sizeHistogram = algorithm.sizeHistogram;
                    configuration.mixedStrategy = strategy;
                    table.push_back(EngineCase{
                        std::string(producer.name) + "/" + algorithm.name + "/" +
                            mixedStrategyName(strategy),
                        configuration});
                }
            }
        }
        for (const bool pruning : {true, false}) {
            MaintenanceConfiguration configuration;
            configuration.runSomei = true;
            configuration.fullMaintenance = true;
            configuration.builds = true;
            configuration.somei.pruning = pruning;
            table.push_back(EngineCase{
                std::string("somei/") + (pruning ? "pruning" : "no pruning"),
                configuration});
        }
        return table;
    }();
    return cases;
}

// Marks the endpoints of `edges` in `appeared`.
void markEndpoints(
    const std::vector<Edge>& catalog,
    std::span<const EdgeId> edges,
    std::vector<unsigned char>& appeared) {
    for (const EdgeId id : edges) {
        const Edge& edge = catalog[static_cast<std::size_t>(id)];
        appeared[static_cast<std::size_t>(edge.u)] = 1U;
        appeared[static_cast<std::size_t>(edge.v)] = 1U;
    }
}

// The maximal cliques the maintenance engine keeps: brute force's, and {v}
// for every vertex that has appeared and has no active edge, which stays in
// the graph as an isolated vertex (see MaintenanceEngine).
[[nodiscard]] CliqueSet engineMaximalCliques(
    int vertexCount,
    const std::vector<Edge>& catalog,
    const std::vector<unsigned char>& active,
    const std::vector<unsigned char>& appeared) {
    CliqueSet cliques = bruteForceMaximalCliques(vertexCount, catalog, active);
    std::vector<unsigned char> hasEdge(appeared.size(), 0U);
    for (std::size_t id = 0; id < catalog.size(); ++id) {
        if (active[id] != 0U) {
            hasEdge[static_cast<std::size_t>(catalog[id].u)] = 1U;
            hasEdge[static_cast<std::size_t>(catalog[id].v)] = 1U;
        }
    }
    for (std::size_t vertex = 0; vertex < appeared.size(); ++vertex) {
        if (appeared[vertex] != 0U && hasEdge[vertex] == 0U) {
            cliques.insert(Clique{static_cast<Vertex>(vertex)});
        }
    }
    return cliques;
}

// One batch of an engine against brute force: the cliques it reports
// created and removed, the number it holds after, and for the subsumed
// search's index the cliques of two or more vertices it holds.
void requireEngineBatch(
    const MaintenanceEngine& engine,
    const BatchChanges& changes,
    const CliqueSet& beforeCliques,
    const CliqueSet& afterCliques,
    bool somei,
    const std::string& label) {
    require(asSet(changes.created, label) ==
                setDifference(afterCliques, beforeCliques),
            label + ": created cliques disagree with brute force");
    require(asSet(changes.removed, label) ==
                setDifference(beforeCliques, afterCliques),
            label + ": removed cliques disagree with brute force");
    require(engine.currentCliqueCount() == afterCliques.size(),
            label + ": the engine holds a different number of cliques");
    if (!somei) {
        CliqueSet indexed;
        for (const Clique& clique : afterCliques) {
            if (clique.size() >= 2U) {
                indexed.insert(clique);
            }
        }
        requireIndexEquals(engine.index(), indexed, label);
    }
}

// The states read a fourth way: 1 marks an edge the batch keeps, 2 one it
// deletes, 3 one it inserts.  Every engine configuration starts from the
// maximal cliques before the batch, runs it as one change and must report
// exactly the cliques it creates and removes, leave the graph after it and
// hold its maximal cliques.
void verifyMixedState(int vertexCount, const std::vector<unsigned char>& state) {
    const std::vector<Edge> catalog = completeCatalog(vertexCount);
    std::vector<EdgeId> beforeEdges;
    std::vector<EdgeId> afterEdges;
    std::vector<EdgeId> inserted;
    std::vector<EdgeId> deleted;
    std::vector<unsigned char> before(catalog.size(), 0U);
    std::vector<unsigned char> after(catalog.size(), 0U);
    for (std::size_t edgeId = 0; edgeId < state.size(); ++edgeId) {
        const auto id = static_cast<EdgeId>(edgeId);
        if (state[edgeId] == 1U || state[edgeId] == 2U) {
            beforeEdges.push_back(id);
            before[edgeId] = 1U;
        }
        if (state[edgeId] == 1U || state[edgeId] == 3U) {
            afterEdges.push_back(id);
            after[edgeId] = 1U;
        }
        if (state[edgeId] == 2U) {
            deleted.push_back(id);
        } else if (state[edgeId] == 3U) {
            inserted.push_back(id);
        }
    }
    if (inserted.empty() && deleted.empty()) {
        return;
    }
    // A vertex the batch inserts an edge at for the first time is an
    // isolated vertex of the graph before it.
    std::vector<unsigned char> appeared(static_cast<std::size_t>(vertexCount), 0U);
    markEndpoints(catalog, beforeEdges, appeared);
    markEndpoints(catalog, inserted, appeared);
    const CliqueSet beforeCliques =
        engineMaximalCliques(vertexCount, catalog, before, appeared);
    const CliqueSet afterCliques =
        engineMaximalCliques(vertexCount, catalog, after, appeared);
    IncrementalGraph reference(vertexCount, catalog);
    reference.applyBatch(afterEdges);
    for (const EngineCase& engineCase : mixedEngineCases()) {
        IncrementalGraph graph(vertexCount, catalog);
        MaintenanceEngine engine(graph, engineCase.configuration);
        if (!beforeEdges.empty()) {
            engine.build(beforeEdges);
        }
        BatchChanges changes;
        engine.updateBatch(0, inserted, deleted, changes);
        const std::string label = "mixed/" + engineCase.label;
        requireEngineBatch(engine, changes, beforeCliques, afterCliques,
                           engineCase.configuration.runSomei, label);
        requireSameGraph(graph, reference, label);
    }
}

void exhaustiveMixedFourVertices() {
    constexpr int vertexCount = 4;
    constexpr int edgeCount = 6;
    std::uint64_t assignments = 1U;
    for (int edge = 0; edge < edgeCount; ++edge) {
        assignments *= 4U;
    }
    for (std::uint64_t code = 0; code < assignments; ++code) {
        std::uint64_t value = code;
        std::vector<unsigned char> state(edgeCount, 0U);
        for (unsigned char& edgeState : state) {
            edgeState = static_cast<unsigned char>(value % 4U);
            value /= 4U;
        }
        verifyMixedState(vertexCount, state);
    }
}

void randomizedMixedFiveVertices() {
    constexpr int vertexCount = 5;
    constexpr int edgeCount = 10;
    std::mt19937 generator(20261011U);
    std::uniform_int_distribution<int> distribution(0, 3);
    for (int trial = 0; trial < 1000; ++trial) {
        std::vector<unsigned char> state(edgeCount, 0U);
        for (unsigned char& edgeState : state) {
            edgeState = static_cast<unsigned char>(distribution(generator));
        }
        verifyMixedState(vertexCount, state);
    }
}

// Long runs of mixed batches on random graphs of five to eight vertices:
// every batch deletes some present edges and inserts some absent ones --
// edges deleted earlier coming back too -- in random batch orders, and each
// is held to brute force.  The configuration rotates with the trial.
void randomMixedStreams() {
    std::mt19937 generator(20261012U);
    const std::vector<EngineCase>& cases = mixedEngineCases();
    for (int trial = 0; trial < static_cast<int>(cases.size()) * 10; ++trial) {
        const EngineCase& engineCase =
            cases[static_cast<std::size_t>(trial) % cases.size()];
        const int vertexCount = std::uniform_int_distribution<int>(5, 8)(generator);
        const std::vector<Edge> catalog = completeCatalog(vertexCount);
        std::bernoulli_distribution initially(
            std::uniform_real_distribution<double>(0.2, 0.9)(generator));
        std::vector<unsigned char> active(catalog.size(), 0U);
        std::vector<EdgeId> initial;
        for (std::size_t id = 0; id < catalog.size(); ++id) {
            if (initially(generator)) {
                initial.push_back(static_cast<EdgeId>(id));
                active[id] = 1U;
            }
        }
        IncrementalGraph graph(vertexCount, catalog);
        MaintenanceEngine engine(graph, engineCase.configuration);
        std::vector<unsigned char> appeared(
            static_cast<std::size_t>(vertexCount), 0U);
        if (!initial.empty()) {
            engine.build(initial);
            markEndpoints(catalog, initial, appeared);
        }
        const std::string label =
            "mixed stream " + std::to_string(trial) + "/" + engineCase.label;
        const int batchCount = std::uniform_int_distribution<int>(4, 10)(generator);
        std::bernoulli_distribution dropping(
            std::uniform_real_distribution<double>(0.1, 0.6)(generator));
        std::bernoulli_distribution adding(
            std::uniform_real_distribution<double>(0.1, 0.6)(generator));
        for (int batch = 0; batch < batchCount; ++batch) {
            const std::vector<unsigned char> activeBefore = active;
            std::vector<EdgeId> inserted;
            std::vector<EdgeId> deleted;
            for (std::size_t id = 0; id < catalog.size(); ++id) {
                if (active[id] != 0U ? dropping(generator) : adding(generator)) {
                    (active[id] != 0U ? deleted : inserted)
                        .push_back(static_cast<EdgeId>(id));
                }
            }
            std::ranges::shuffle(inserted, generator);
            std::ranges::shuffle(deleted, generator);
            for (const EdgeId id : inserted) {
                active[static_cast<std::size_t>(id)] = 1U;
            }
            for (const EdgeId id : deleted) {
                active[static_cast<std::size_t>(id)] = 0U;
            }
            // The batch's new endpoints are isolated vertices before it.
            markEndpoints(catalog, inserted, appeared);
            const CliqueSet beforeCliques =
                engineMaximalCliques(vertexCount, catalog, activeBefore, appeared);
            const CliqueSet afterCliques =
                engineMaximalCliques(vertexCount, catalog, active, appeared);
            BatchChanges changes;
            engine.updateBatch(batch, inserted, deleted, changes);
            requireEngineBatch(engine, changes, beforeCliques, afterCliques,
                               engineCase.configuration.runSomei,
                               label + "/batch " + std::to_string(batch));
        }
    }
}

// removeBatch checks the whole batch before it changes anything.
void removeBatchRejectsBadInput() {
    const std::vector<Edge> catalog = completeCatalog(4);
    IncrementalGraph graph(4, catalog);
    graph.applyBatch(std::vector<EdgeId>{0, 1, 2});
    const auto rejects = [&](std::vector<EdgeId> batch, const std::string& what) {
        bool threw = false;
        try {
            graph.removeBatch(batch);
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, "removeBatch accepted " + what);
        require(graph.activeEdgeCount() == 3,
                "a rejected removeBatch changed the graph");
    };
    rejects({3}, "an edge that is not in the graph");
    rejects({0, 0}, "a repeated edge");
    rejects({0, 99}, "an edge outside the catalog");
    graph.removeBatch(std::vector<EdgeId>{0, 1, 2});
    require(graph.activeEdgeCount() == 0 && graph.activeVertices().empty(),
            "removing every edge left something behind");
    graph.applyBatch(std::vector<EdgeId>{0});
    require(graph.activeEdgeCount() == 1 && graph.isActive(0),
            "a removed edge could not be applied again");
}

void exhaustiveFourVertices() {
    constexpr int vertexCount = 4;
    constexpr int edgeCount = 6;
    std::uint64_t assignments = 1U;
    for (int edge = 0; edge < edgeCount; ++edge) {
        assignments *= 3U;
    }
    for (std::uint64_t code = 0; code < assignments; ++code) {
        std::uint64_t value = code;
        std::vector<unsigned char> state(edgeCount, 0U);
        for (unsigned char& edgeState : state) {
            edgeState = static_cast<unsigned char>(value % 3U);
            value /= 3U;
        }
        verifyState(vertexCount, state);
        verifyDeletionState(vertexCount, state);
    }
}

void randomizedFiveVertices() {
    constexpr int vertexCount = 5;
    constexpr int edgeCount = 10;
    std::mt19937 generator(20260831U);
    std::uniform_int_distribution<int> distribution(0, 2);
    for (int trial = 0; trial < 1000; ++trial) {
        std::vector<unsigned char> state(edgeCount, 0U);
        for (unsigned char& edgeState : state) {
            edgeState = static_cast<unsigned char>(distribution(generator));
        }
        verifyState(vertexCount, state);
        verifyDeletionState(vertexCount, state);
    }
}

// The new-clique search alone against brute force -- the edge and vertex
// baselines and both star methods -- for states past the sizes the subsumed
// table runs on.
// `state` holds 0 for an edge that never appears, 1 for one present before
// the batch and 2 for one the batch inserts.
void verifyNewCliques(int vertexCount, const std::vector<unsigned char>& state) {
    const std::vector<Edge> catalog = completeCatalog(vertexCount);
    std::vector<EdgeId> oldEdges;
    std::vector<EdgeId> batch;
    std::vector<unsigned char> before(catalog.size(), 0U);
    std::vector<unsigned char> after(catalog.size(), 0U);
    for (std::size_t edgeId = 0; edgeId < state.size(); ++edgeId) {
        if (state[edgeId] == 1U) {
            oldEdges.push_back(static_cast<EdgeId>(edgeId));
            before[edgeId] = 1U;
            after[edgeId] = 1U;
        } else if (state[edgeId] == 2U) {
            batch.push_back(static_cast<EdgeId>(edgeId));
            after[edgeId] = 1U;
        }
    }
    if (batch.empty()) {
        return;
    }
    CliqueSet expected = bruteForceMaximalCliques(vertexCount, catalog, after);
    for (const Clique& clique :
         bruteForceMaximalCliques(vertexCount, catalog, before)) {
        expected.erase(clique);
    }

    IncrementalGraph graph(vertexCount, catalog);
    graph.applyBatch(oldEdges);
    graph.applyBatch(batch);
    const std::string where = " on a " + std::to_string(vertexCount) +
        "-vertex graph";
    NewCliqueEnumerator edge(graph);
    NewCliqueSearchStats edgeStats;
    require(asSet(edge.enumerateBatch(batch, edgeStats), "edge") == expected,
            "edge disagrees with brute force" + where);
    VertexCliqueEnumerator vertex(graph);
    VertexSearchStats vertexStats;
    require(asSet(vertex.enumerateBatch(batch, vertexStats), "vertex") ==
                expected,
            "vertex disagrees with brute force" + where);
    for (const VertexSearchConfiguration& method : vertexAblationMethods()) {
        const std::string label = vertexAblationLabel(method);
        VertexCliqueEnumerator ablated(graph, method);
        VertexSearchStats ablatedStats;
        require(asSet(ablated.enumerateBatch(batch, ablatedStats), label) ==
                    expected,
                label + " disagrees with brute force" + where);
    }
    for (const StarSearchConfiguration& method : starMethods()) {
        const std::string label =
            std::string("star/") + starBranchScopeName(method.branchScope);
        require(asSet(enumerateStar(graph, batch, method), label) == expected,
                label + " disagrees with brute force" + where);
    }
    for (const StarSearchConfiguration& method : starAblationMethods()) {
        const std::string label = starAblationLabel(method);
        require(asSet(enumerateStar(graph, batch, method), label) == expected,
                label + " disagrees with brute force" + where);
    }
}

// Six vertices are past exhaustive enumeration (3^15 states), so they are
// sampled. It is the first size at which a centre can own five edges.
void randomizedSixVertices() {
    constexpr int vertexCount = 6;
    constexpr int edgeCount = 15;
    std::mt19937 generator(20260917U);
    std::uniform_int_distribution<int> distribution(0, 2);
    for (int trial = 0; trial < 1500; ++trial) {
        std::vector<unsigned char> state(edgeCount, 0U);
        for (unsigned char& edgeState : state) {
            edgeState = static_cast<unsigned char>(distribution(generator));
        }
        verifyNewCliques(vertexCount, state);
    }
}

// Two hubs sharing a rim, which is where stars compete for the same cliques:
// the rim an old path, every rim vertex joined to one hub by a batch edge and
// to the other by an old one, the hubs joined by an old edge. Every split of
// the rim between the two hubs, for rims of two to six vertices.
void doubleHubStates() {
    for (int rim = 2; rim <= 6; ++rim) {
        const int vertexCount = rim + 2;
        const std::vector<Edge> catalog = completeCatalog(vertexCount);
        const Vertex hubA = 0;
        const Vertex hubB = 1;
        for (std::uint32_t rimMask = 0;
             rimMask < (std::uint32_t{1} << rim); ++rimMask) {
            std::vector<unsigned char> state(catalog.size(), 0U);
            const auto set = [&](Vertex u, Vertex v, unsigned char value) {
                const Edge wanted = normalizedEdge(u, v);
                const auto found =
                    std::lower_bound(catalog.begin(), catalog.end(), wanted);
                state[static_cast<std::size_t>(found - catalog.begin())] =
                    value;
            };
            for (int index = 0; index + 1 < rim; ++index) {
                set(2 + index, 3 + index, 1U);
            }
            set(hubA, hubB, 1U);
            for (int index = 0; index < rim; ++index) {
                const bool toA = (rimMask >> index & 1U) != 0U;
                set(toA ? hubA : hubB, 2 + index, 2U);
                set(toA ? hubB : hubA, 2 + index, 1U);
            }
            verifyNewCliques(vertexCount, state);
        }
    }
}

// The vertex baseline where ownership and the new-edge condition both act.
// The batch (0, 3), (1, 2) ranks the endpoints 0, 3, 1, 2, and old edges join
// 0 to 1, 2 and 4. Root 0 owns {0, 1, 2} through (1, 2), a batch edge not its
// own, and {0, 3}; its child on 4 reaches no batch edge and is cut, where the
// unpruned search would reach the old clique {0, 4}. Roots 3 and 2 have no
// later neighbour and are cut before their graphs are built; root 1 builds
// one and finds its only clique already owned by 0.
// The seed-cover ablation (VertexSeeds::Cover).  The batch (0, 1), (1, 2),
// (3, 4), (2, 4) puts 0, 1, 3 and 4 into the cover but not 2, which appears
// before 4: as no seed, 2 stays in P at root 4, which owns {2, 4}.
void vertexCoverSeedsLeaveOtherEndpointsInP() {
    const std::vector<Edge> catalog{{0, 1}, {1, 2}, {3, 4}, {2, 4}};
    IncrementalGraph graph(5, catalog);
    const std::vector<EdgeId> batch{0, 1, 2, 3};
    graph.applyBatch(batch);
    const CliqueSet expected{{0, 1}, {1, 2}, {3, 4}, {2, 4}};
    VertexCliqueEnumerator endpoints(graph);
    VertexSearchStats endpointStats;
    require(asSet(endpoints.enumerateBatch(batch, endpointStats), "vertex") ==
                expected,
            "vertex found the wrong new cliques");
    require(endpointStats.roots == 5, "vertex did not seed every endpoint");
    VertexCliqueEnumerator cover(
        graph,
        VertexSearchConfiguration{VertexBranching::Pivot, VertexSeeds::Cover});
    VertexSearchStats coverStats;
    require(asSet(cover.enumerateBatch(batch, coverStats), "vertex/cover") ==
                expected,
            "vertex/cover found the wrong new cliques");
    require(coverStats.roots == 4, "vertex/cover did not seed the cover alone");
}

void vertexRootsOwnOnlyTheirCliques() {
    const std::vector<Edge> catalog{{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}};
    IncrementalGraph graph(5, catalog);
    const std::vector<EdgeId> oldEdges{0, 1, 3};
    const std::vector<EdgeId> batch{2, 4};
    graph.applyBatch(oldEdges);
    graph.applyBatch(batch);
    VertexCliqueEnumerator enumerator(graph);
    VertexSearchStats stats;
    require(asSet(enumerator.enumerateBatch(batch, stats), "vertex") ==
                CliqueSet{{0, 1, 2}, {0, 3}},
            "vertex found the wrong new cliques");
    require(stats.roots == 4 && stats.prunedRoots == 2 &&
                stats.localGraph.builds == 2,
            "vertex did not cut the roots that own nothing");
    require(stats.newEdgePrunes == 1 && stats.leaves == 2,
            "vertex did not cut the branch that reaches no batch edge");
    // The ablations keep the ownership by seeds (here the cover holds every
    // endpoint): root 0's new edges include (1, 2), which it does not touch,
    // or {0, 1, 2} is lost.
    for (const VertexSearchConfiguration& method : vertexAblationMethods()) {
        const std::string label = vertexAblationLabel(method);
        VertexCliqueEnumerator ablated(graph, method);
        VertexSearchStats ablatedStats;
        require(asSet(ablated.enumerateBatch(batch, ablatedStats), label) ==
                    CliqueSet{{0, 1, 2}, {0, 3}},
                label + " found the wrong new cliques");
        require(ablatedStats.roots == 4 && ablatedStats.prunedRoots == 2,
                label + " did not cut the roots that own nothing");
    }
    VertexCliqueEnumerator newEdges(
        graph, VertexSearchConfiguration{VertexBranching::NewEdges});
    VertexSearchStats newEdgeStats;
    static_cast<void>(newEdges.enumerateBatch(batch, newEdgeStats));
    require(newEdgeStats.rootEdgeChoices == 2 &&
                newEdgeStats.edgeChildren == 3 &&
                newEdgeStats.skippedEdgeChildren == 0,
            "vertex/new-edges did not open one child per edge root 0 and root 1 own");
}

// What re-indexing may not change: the search itself, node for node. Pivots
// are chosen by global vertex order, which a re-indexed graph keeps, so even
// the pivot scans are the same.
void requireSameSearch(
    const StarSearchStats& plain,
    const StarSearchStats& moved,
    const std::string& label) {
    require(plain.recursiveCalls == moved.recursiveCalls &&
                plain.pivotCalls == moved.pivotCalls &&
                plain.pivotCandidatesEvaluated ==
                    moved.pivotCandidatesEvaluated &&
                plain.leaves == moved.leaves &&
                plain.partnerPrunes == moved.partnerPrunes &&
                plain.emittedCliques == moved.emittedCliques &&
                plain.branchRootEdgeStars == moved.branchRootEdgeStars &&
                plain.branchRootPivotStars == moved.branchRootPivotStars &&
                plain.branchUnboundBindings == moved.branchUnboundBindings,
            label + " changed the search");
}

// The re-indexing gate shut: no subtree is ever big enough to move.
[[nodiscard]] StarSearchConfiguration withGateClosed(
    StarSearchConfiguration configuration) {
    configuration.reindexMinCandidates = std::numeric_limits<std::size_t>::max();
    return configuration;
}

// Star domains several words wide, so that every re-indexing gate opens: the
// per-partner trees, the pivot partition's first-level subtrees and the levels
// below them all move onto graphs of their own. Three hubs reach most of the
// graph and each owns a few dozen batch edges, the edges between the hubs are
// new as well, and the batch is shuffled, so the stars also claim cliques from
// one another. Both methods must enumerate exactly the reference's new
// cliques, and re-indexing -- at the default gate or opened all the way --
// must leave the search exactly as it is with the gate shut.
void starWideDomainsMatchReference() {
    constexpr int vertexCount = 220;
    constexpr Vertex hubCount = 3;
    std::mt19937 generator(20261005U);
    std::bernoulli_distribution sparse(0.12);
    std::bernoulli_distribution hubLink(0.7);
    std::bernoulli_distribution hubLinkIsNew(0.3);
    std::bernoulli_distribution otherIsNew(0.02);

    std::vector<Edge> catalog;
    std::vector<Edge> oldEdges;
    std::vector<EdgeId> oldIds;
    std::vector<EdgeId> batch;
    for (Vertex u = 0; u < vertexCount; ++u) {
        for (Vertex v = u + 1; v < vertexCount; ++v) {
            bool present = false;
            bool isNew = false;
            if (v < hubCount) {
                present = true;
                isNew = true;
            } else if (u < hubCount) {
                present = hubLink(generator);
                isNew = present && hubLinkIsNew(generator);
            } else {
                present = sparse(generator);
                isNew = present && otherIsNew(generator);
            }
            if (!present) {
                continue;
            }
            const auto id = static_cast<EdgeId>(catalog.size());
            catalog.push_back({u, v});
            if (isNew) {
                batch.push_back(id);
            } else {
                oldIds.push_back(id);
                oldEdges.push_back({u, v});
            }
        }
    }
    std::ranges::shuffle(batch, generator);

    const CliqueSet before = referenceMaximalCliques(vertexCount, oldEdges);
    CliqueSet expectedNew = referenceMaximalCliques(vertexCount, catalog);
    for (const Clique& clique : before) {
        expectedNew.erase(clique);
    }

    IncrementalGraph graph(vertexCount, catalog);
    graph.applyBatch(oldIds);
    graph.applyBatch(batch);
    NewCliqueEnumerator edgeEnumerator(graph);
    NewCliqueSearchStats edgeStats;
    require(asSet(edgeEnumerator.enumerateBatch(batch, edgeStats), "wide edge") ==
                expectedNew,
            "wide domains: edge search disagrees with the reference");

    const auto run = [&](const std::string& label,
                         const StarSearchConfiguration& configuration) {
        StarCliqueEnumerator enumerator(graph, configuration);
        StarSearchStats stats;
        require(asSet(enumerator.enumerateBatch(batch, stats), label) ==
                    expectedNew,
                "wide domains: star/" + label + " disagrees with the reference");
        require(stats.maximumDomainVertices > 64U,
                "wide domains: star/" + label +
                    " never built a domain wider than one word");
        requireStarAccounts(stats, "wide domains: star/" + label);
        return stats;
    };

    long long pivotRoots = 0;
    long long edgeRoots = 0;
    for (const StarSearchConfiguration& method : starMethods()) {
        const std::string name = starBranchScopeName(method.branchScope);
        StarSearchConfiguration eager = method;
        eager.reindexMinCandidates = 1U;
        const StarSearchStats plain = run(name + "/gate-shut", withGateClosed(method));
        const StarSearchStats gated = run(name, method);
        const StarSearchStats opened = run(name + "/gate-open", eager);
        StarSearchConfiguration unfiltered = method;
        unfiltered.rootFilter = false;
        const StarSearchStats ablated = run(name + "/no-root-filter", unfiltered);
        require(ablated.domainVertexOccurrences > gated.domainVertexOccurrences,
                "wide domains: " + name +
                    " with the root filter off did not widen its domains");
        require(plain.reindexedTrees == 0 && plain.mergedReindexedSubtrees == 0 &&
                    plain.deepReindexedSubtrees == 0,
                "wide domains: " + name + " moved a subtree with the gate shut");
        require(gated.reindexedTrees + gated.mergedReindexedSubtrees > 0 &&
                    gated.deepReindexedSubtrees > 0,
                "wide domains: " + name + " moved no subtree at the default gate");
        require(opened.deepReindexedSubtrees > 0,
                "wide domains: " + name + " moved no subtree with the gate open");
        requireSameSearch(plain, gated, "wide domains: re-indexing of " + name);
        requireSameSearch(plain, opened, "wide domains: eager re-indexing of " + name);
        require(gated.nodeWords < plain.nodeWords &&
                    gated.pivotWords <= plain.pivotWords,
                "wide domains: re-indexing of " + name + " did not narrow the search");
        if (method.branchScope == StarBranchScope::Root) {
            require(gated.branchUnboundBindings == 0 &&
                        gated.branchUnitBindings == 0,
                    "wide domains: the root scope decided below a root");
        }
        pivotRoots += gated.branchRootPivotStars;
        edgeRoots += gated.branchRootEdgeStars;
    }
    require(pivotRoots > 0 && edgeRoots > 0,
            "wide domains: the rule never took one of the two partitions");

    // The ablation's fixed partitions take their own side at every centre and
    // never ask the rule below one.
    for (const StarSearchConfiguration& method : starAblationMethods()) {
        const std::string name = starAblationLabel(method).substr(5);
        const StarSearchStats stats = run(name, method);
        require(stats.branchUnboundBindings == 0 && stats.branchUnitBindings == 0,
                "wide domains: " + name + " asked the rule below a centre");
        if (method.branching == StarBranching::Anchor) {
            require(stats.branchRootPivotStars == 0 && stats.branchRootEdgeStars > 0,
                    "wide domains: " + name + " took the pivot at a centre");
        } else {
            require(stats.branchRootEdgeStars == 0 && stats.branchRootPivotStars > 0,
                    "wide domains: " + name + " took the new edges at a centre");
        }
    }
}

// One star built so that the stack of nested graphs grows while an unbound
// state above it still reads its own level's partner mask. Centre 0; partners
// 4..253, joined to the centre by the batch and to nothing else among
// themselves; old vertices a = 1, b = 2, c = 3, adjacent to the centre and to
// each other, with a adjacent to the first 180 partners, b to the first 100
// and c to the first 40. The centre's pivot is a (182 of its 253 candidates,
// against 250 partners), so the pivot partition is taken, and a's subtree --
// three words in the star's four -- moves to the first-level graph. There the
// pivot is b, whose subtree (two words) moves to the first nested level, and
// inside that the pivot is c, whose subtree (one word) moves to the second:
// created while the state on the first level still has its siblings to
// classify through its own mask. Output and search must match the run with
// the gate shut.
void starNestedGraphsKeepTheirMasks() {
    constexpr Vertex centre = 0;
    constexpr Vertex a = 1;
    constexpr Vertex b = 2;
    constexpr Vertex c = 3;
    constexpr Vertex firstPartner = 4;
    constexpr int partnerCount = 250;
    constexpr int vertexCount = firstPartner + partnerCount;

    std::vector<Edge> catalog;
    std::vector<Edge> oldEdges;
    std::vector<EdgeId> oldIds;
    std::vector<EdgeId> batch;
    const auto add = [&](Vertex u, Vertex v, bool isNew) {
        const auto id = static_cast<EdgeId>(catalog.size());
        catalog.push_back(normalizedEdge(u, v));
        if (isNew) {
            batch.push_back(id);
        } else {
            oldIds.push_back(id);
            oldEdges.push_back(normalizedEdge(u, v));
        }
    };
    for (const Vertex core : {a, b, c}) {
        add(centre, core, false);
    }
    add(a, b, false);
    add(a, c, false);
    add(b, c, false);
    for (int index = 0; index < partnerCount; ++index) {
        const Vertex partner = firstPartner + index;
        add(centre, partner, true);
        if (index < 180) {
            add(a, partner, false);
        }
        if (index < 100) {
            add(b, partner, false);
        }
        if (index < 40) {
            add(c, partner, false);
        }
    }

    const CliqueSet before = referenceMaximalCliques(vertexCount, oldEdges);
    CliqueSet expectedNew = referenceMaximalCliques(vertexCount, catalog);
    for (const Clique& clique : before) {
        expectedNew.erase(clique);
    }
    require(expectedNew.size() == static_cast<std::size_t>(partnerCount),
            "nested graphs: unexpected new cliques in the reference");
    IncrementalGraph graph(vertexCount, catalog);
    graph.applyBatch(oldIds);
    graph.applyBatch(batch);

    for (const StarSearchConfiguration& method : starMethods()) {
        const std::string name = starBranchScopeName(method.branchScope);
        StarSearchStats plain;
        StarSearchStats nested;
        for (const bool shut : {true, false}) {
            StarCliqueEnumerator enumerator(
                graph, shut ? withGateClosed(method) : method);
            StarSearchStats& stats = shut ? plain : nested;
            require(asSet(enumerator.enumerateBatch(batch, stats),
                          "nested masks") == expectedNew,
                    "nested graphs: " + name + " disagrees with the reference");
        }
        require(nested.branchRootPivotStars == 1 &&
                    nested.mergedReindexedSubtrees == 1 &&
                    nested.deepReindexedSubtrees == 2 &&
                    nested.reindexedTrees == 0,
                "nested graphs: " + name +
                    " did not nest the subtrees two levels deep");
        requireSameSearch(plain, nested, "nested graphs: nesting under " + name);
    }
}

void duplicateSubsumedCliqueIsEmittedOnce() {
    const std::vector<Edge> catalog = completeCatalog(4);
    const std::vector<EdgeId> oldEdges{0}; // (0,1)
    const std::vector<EdgeId> batch{1, 2, 3, 4};
    IncrementalGraph graph(4, catalog);
    graph.applyBatch(oldEdges);

    graph.applyBatch(batch);

    NewCliqueEnumerator newEnumerator(graph);
    NewCliqueSearchStats newStats;
    const std::vector<Clique> newCliques =
        newEnumerator.enumerateBatch(batch, newStats);
    require(asSet(newCliques, "overlap new") ==
                CliqueSet{Clique{0, 1, 2}, Clique{0, 1, 3}},
            "overlap setup produced unexpected new cliques");

    for (const Variant& variant : variants()) {
        MaximalCliqueIndex index(
            usesSubsumedSizeGate(variant.algorithm) && variant.sizeHistogram);
        require(index.insert(Clique{0, 1}), "failed to seed old clique");
        SubsumedCliqueEnumerator subsumedEnumerator(graph);
        if (variant.forceDeferral) {
            subsumedEnumerator.setDeferralMinimumIndexSize(0U);
        }
        SubsumedCliqueSearchStats stats;
        const std::vector<Clique> deleted = runVariant(
            subsumedEnumerator, variant, batch, newCliques, index, stats);
        require(asSet(deleted, "overlap deleted") == CliqueSet{Clique{0, 1}},
                "an old clique contained by two new cliques must be emitted once");
        require(stats.emittedCliques == 1,
                "subsumed statistics must count the old clique once");
    }
}

void wideCliqueUsesCorrectRepresentation(int vertexCount) {
    const std::vector<Edge> catalog{{0, 1}};
    const std::vector<EdgeId> batch{0};
    IncrementalGraph graph(vertexCount, catalog);
    graph.applyBatch(batch);

    Clique newClique(static_cast<std::size_t>(vertexCount));
    std::iota(newClique.begin(), newClique.end(), 0);
    Clique withoutZero(newClique.begin() + 1, newClique.end());
    Clique withoutOne;
    withoutOne.reserve(newClique.size() - 1U);
    withoutOne.push_back(0);
    withoutOne.insert(withoutOne.end(), newClique.begin() + 2, newClique.end());
    const CliqueSet expected{withoutZero, withoutOne};

    long long fullBKNodes = 0;
    long long kernelBKNodes = 0;
    for (SubsumedCliqueAlgorithm algorithm : {
             SubsumedCliqueAlgorithm::Optimized,
             SubsumedCliqueAlgorithm::ResidualKernel,
             SubsumedCliqueAlgorithm::OptimizedBitmap,
             SubsumedCliqueAlgorithm::FullBK,
             SubsumedCliqueAlgorithm::KernelBK}) {
        MaximalCliqueIndex index;
        require(index.insert(withoutZero), "failed to seed wide old clique");
        require(index.insert(withoutOne), "failed to seed second wide old clique");
        SubsumedCliqueEnumerator enumerator(graph);
        SubsumedCliqueSearchStats stats;
        const std::array<Clique, 1U> newCliques{newClique};
        const std::vector<Clique> deleted = enumerator.enumerateBatch(
            batch,
            newCliques,
            index,
            stats,
            algorithm);
        require(asSet(deleted, "wide deleted") == expected,
                "wide clique search produced incorrect candidates");
        if (algorithm == SubsumedCliqueAlgorithm::Optimized) {
            const bool expectsFallback = vertexCount > 256;
            require(stats.fallbackCliques == (expectsFallback ? 1 : 0),
                    "wide clique used the wrong candidate representation");
            require(stats.bitmaskCliques == (expectsFallback ? 0 : 1),
                    "wide clique bitmask statistics are incorrect");
        }
        if (isBKComparison(algorithm)) {
            require(stats.indexLookups == 2, "BK comparison must probe both candidates");
            if (vertexCount <= 64) {
                require(stats.fastPathCliques == 1 && stats.wide.kernelCalls == 0,
                        "BK comparison must share the one-edge shortcut");
            } else {
                require(stats.wide.kernelCalls == 1 && stats.wide.leaves == 2,
                        "BK comparison must enumerate two leaves on a wide clique");
                if (algorithm == SubsumedCliqueAlgorithm::FullBK) {
                    fullBKNodes = stats.wide.nodes;
                } else {
                    kernelBKNodes = stats.wide.nodes;
                }
            }
        }
    }
    if (vertexCount > 64) {
        require(kernelBKNodes > 0 && fullBKNodes - kernelBKNodes == vertexCount - 2,
                "FullBK must visit one extra node per universal vertex");
    }
}

// The index survives many table growths, and each one is allowed to throw
// away the dead records and renumber the survivors underneath the table it is
// re-laying.  No dataset reaches the >50%-dead threshold that triggers the
// renumbering, so the path is driven here instead: churn enough cliques
// through the index that both the plain rehash and the compacting one run,
// and check the contents against an independent set after every phase.
void indexSurvivesRehashAndCompaction() {
    constexpr int cliqueCount = 40000;
    MaximalCliqueIndex index(true);
    CliqueSet live;
    const auto cliqueAt = [](int seed) {
        return Clique{seed, seed + 1, seed + 3};
    };
    const auto checkAgrees = [&](const std::string& phase) {
        require(index.size() == live.size(),
                phase + ": index size disagrees with the reference set");
        for (const Clique& clique : live) {
            require(index.containsCanonical(clique),
                    phase + ": index lost a live clique");
        }
        std::size_t histogram = index.oversizeLiveCliqueCount();
        for (std::size_t k = 0; k < MaximalCliqueIndex::trackedSizeLimit; ++k) {
            histogram += index.liveCliqueCountOfSize(k);
        }
        require(histogram == index.size(),
                phase + ": size histogram disagrees with the live count");
    };

    for (int seed = 0; seed < cliqueCount; seed += 5) {
        const Clique clique = cliqueAt(seed);
        require(index.insert(clique), "index rejected a fresh clique");
        live.insert(clique);
    }
    checkAgrees("after growth");

    // Erase most of them, so the next growth has to compact.
    int erased = 0;
    for (int seed = 0; seed < cliqueCount; seed += 5) {
        if (erased % 5 == 4) {
            ++erased;
            continue;
        }
        const Clique clique = cliqueAt(seed);
        require(index.eraseCanonical(clique), "index lost a clique to erase");
        live.erase(clique);
        ++erased;
    }
    checkAgrees("after bulk erase");

    // Re-grow past the old table size, which forces the compacting rehash.
    for (int seed = 1; seed < cliqueCount; seed += 5) {
        const Clique clique = cliqueAt(seed);
        require(index.insert(clique), "index rejected a fresh clique");
        live.insert(clique);
    }
    checkAgrees("after regrowth");

    // Everything still erasable, including across the renumbering.
    for (const Clique& clique : CliqueSet(live)) {
        require(index.eraseCanonical(clique), "index lost a clique to erase");
        live.erase(clique);
    }
    checkAgrees("after drain");
    require(index.empty(), "drained index is not empty");
    require(index.minimumLiveCliqueSize() == static_cast<std::size_t>(-1),
            "drained index still reports a minimum live size");
}

// Residual kernel search on kernels the exhaustive small states never reach:
// |T| from 12 past 64 (the complement search over clique positions, then the
// multi-word search), |S| past 64 (the anchor lists replace the pair tests),
// |C| past 64 and 256 (multi-word masks and inline probes), a single kernel
// edge inside a wide clique, and batch endpoints in C whose batch edges leave
// it -- S larger than T.
//
// Vertices [0, n) form the new clique C: complete in G + H, with the batch
// edges inside C a random graph on a random t-subset.  Each outside vertex
// has old edges to all of C but two vertices: one it misses entirely (so C
// stays maximal) and one it meets by a batch edge (which leaves C).  Every
// maximal clique of the old graph inside C is subsumed by C; a random half
// of them is seeded into the index, and each framework has to delete exactly
// that half -- the other half are candidates that must miss, as are the
// kernel's independent sets that an outside vertex extends in the old graph.
void residualKernelRandomKernels() {
    std::mt19937_64 random(20260924U);
    struct Shape {
        int n;
        int t;
        double density;
        int extra;
        // Also run the controls, which split the kernel edge by edge: only
        // where their intermediate candidate sets stay small enough for a
        // unit test.
        bool control;
    };
    const std::array shapes{
        Shape{16, 12, 0.30, 0, true},
        Shape{20, 2, 1.00, 3, true},
        Shape{24, 16, 0.50, 3, true},
        Shape{40, 20, 0.25, 2, true},
        Shape{64, 30, 0.60, 4, false},
        Shape{64, 64, 0.97, 0, false},
        Shape{70, 20, 0.40, 2, false},
        Shape{80, 2, 1.00, 0, true},
        Shape{100, 70, 0.98, 3, false},
        Shape{130, 40, 0.90, 5, false},
        Shape{300, 12, 0.50, 2, true},
    };
    struct Arm {
        const char* name;
        SubsumedCliqueAlgorithm algorithm;
        bool sizeHistogram = true;
    };
    constexpr auto rk = SubsumedCliqueAlgorithm::ResidualKernel;
    const std::array arms{
        Arm{.name = "optimized", .algorithm = SubsumedCliqueAlgorithm::Optimized},
        Arm{.name = "optimized-bitmap",
            .algorithm = SubsumedCliqueAlgorithm::OptimizedBitmap},
        Arm{.name = "residual-kernel", .algorithm = rk},
        Arm{.name = "residual-kernel/no-size-histogram", .algorithm = rk,
            .sizeHistogram = false},
    };
    long long residualSearches = 0;
    long long wideKernels = 0;
    for (const Shape& shape : shapes) {
        int trials = 0;
        for (int attempt = 0; attempt < 40 && trials < 4; ++attempt) {
            const int n = shape.n;
            const int total = n + shape.extra;
            const auto catalog = completeCatalog(total);
            std::vector<Vertex> order(static_cast<std::size_t>(n));
            std::iota(order.begin(), order.end(), 0);
            std::shuffle(order.begin(), order.end(), random);
            std::vector<unsigned char> inKernel(static_cast<std::size_t>(n), 0U);
            for (int i = 0; i < shape.t; ++i) {
                inKernel[static_cast<std::size_t>(
                    order[static_cast<std::size_t>(i)])] = 1U;
            }
            std::bernoulli_distribution conflict(shape.density);
            std::bernoulli_distribution seed(0.5);
            std::uniform_int_distribution<Vertex> inside(0, n - 1);
            std::vector<Vertex> misses;
            std::vector<Vertex> partners;
            for (int x = n; x < total; ++x) {
                misses.push_back(inside(random));
                Vertex partner = inside(random);
                while (partner == misses.back()) {
                    partner = inside(random);
                }
                partners.push_back(partner);
            }
            std::vector<EdgeId> oldIds;
            std::vector<EdgeId> newIds;
            std::vector<Edge> oldEdges;
            bool internal = false;
            for (std::size_t id = 0; id < catalog.size(); ++id) {
                const Edge edge = catalog[id];
                const auto edgeId = static_cast<EdgeId>(id);
                if (edge.v < n) {
                    if (inKernel[static_cast<std::size_t>(edge.u)] != 0U &&
                        inKernel[static_cast<std::size_t>(edge.v)] != 0U &&
                        conflict(random)) {
                        newIds.push_back(edgeId);
                        internal = true;
                    } else {
                        oldIds.push_back(edgeId);
                        oldEdges.push_back(edge);
                    }
                } else if (edge.u < n) {
                    const auto outside = static_cast<std::size_t>(edge.v - n);
                    if (edge.u == misses[outside]) {
                        continue;
                    }
                    if (edge.u == partners[outside]) {
                        newIds.push_back(edgeId);
                    } else {
                        oldIds.push_back(edgeId);
                        oldEdges.push_back(edge);
                    }
                }
            }
            if (!internal) {
                continue;
            }
            ++trials;
            IncrementalGraph graph(total, catalog);
            graph.applyBatch(oldIds);
            CliqueSet seeded;
            for (const Clique& clique : referenceMaximalCliques(total, oldEdges)) {
                if (clique.back() < n && seed(random)) {
                    seeded.insert(clique);
                }
            }
            graph.applyBatch(newIds);
            Clique full(static_cast<std::size_t>(n));
            std::iota(full.begin(), full.end(), 0);
            for (const Arm& arm : arms) {
                if (arm.algorithm != rk && !shape.control) {
                    continue;
                }
                MaximalCliqueIndex index(
                    usesSubsumedSizeGate(arm.algorithm) && arm.sizeHistogram);
                for (const Clique& clique : seeded) {
                    require(index.insert(clique), "residual kernel seed");
                }
                SubsumedCliqueEnumerator enumerator(graph);
                SubsumedCliqueSearchStats stats;
                const std::vector<Clique> deleted = enumerator.enumerateBatch(
                    newIds, std::array<Clique, 1U>{full}, index, stats,
                    arm.algorithm);
                const std::string where = std::string("residual kernel n=") +
                    std::to_string(n) + " t=" + std::to_string(shape.t) +
                    " trial " + std::to_string(trials) + ": " + arm.name;
                require(asSet(deleted, arm.name) == seeded,
                        where + " deleted the wrong set");
                if (arm.algorithm != rk) {
                    continue;
                }
                // Only a split ever holds a candidate set.
                require(stats.maximumCandidates == 0U &&
                            stats.fallbackCliques == 0,
                        where + ": a kernel was split");
                if (arm.sizeHistogram) {
                    residualSearches += stats.residualKernelSearches;
                } else {
                    require(stats.sizeGateSkippedCliques == 0 &&
                                stats.sizeGateSkippedProbes == 0,
                            where + ": the size gate fired without a histogram");
                }
            }
            if (shape.t > 64) {
                ++wideKernels;
            }
        }
        require(trials == 4, "residual kernel shape produced no kernel");
    }
    require(residualSearches > 0,
            "the complement search over clique positions never ran");
    require(wideKernels > 0, "no kernel wider than one word was generated");
}

// The pivot scan's early exits in the three complement searches -- the one in
// `consumeKernelRows64`, `wideresidual::runNarrow` and the multi-word
// `WideSearch::visit` -- on a kernel built to reach all three.  On T = [0, t)
// the complement of the conflict graph is a K4 {0..3}, a triangle {4, 5, 6}
// and isolated vertices, so: the root pivot is vertex 0; inside the K4 and the
// triangle every pivot scores |P|-1 with X empty; vertex 5 meets 4 in X
// covering its P = {6}; vertex 6 has P empty and X = {4, 5}.  The random
// kernels above never reach the X-covers-P exit on two words.  One shape per
// search: |C| <= 64; |C| > 64 with |T| <= 64; |T| > 64.  The old maximal
// cliques inside C are K plus one of the complement's maximal cliques; all
// are seeded but every other isolated vertex's, which must miss, and seeding
// one isolated vertex's keeps the size bound from pruning any branch.
void residualKernelPivotEarlyExits() {
    struct Shape {
        int n;
        int t;
        // Also run the controls, which split: only on the narrow kernel.
        bool control;
    };
    const std::array shapes{
        Shape{12, 10, true},
        Shape{70, 10, true},
        Shape{70, 67, false},
    };
    struct Arm {
        const char* name;
        SubsumedCliqueAlgorithm algorithm;
    };
    const std::array arms{
        Arm{"residual-kernel", SubsumedCliqueAlgorithm::ResidualKernel},
        Arm{"optimized", SubsumedCliqueAlgorithm::Optimized},
        Arm{"optimized-bitmap", SubsumedCliqueAlgorithm::OptimizedBitmap},
    };
    // Pairs of T that are not batch edges.
    const auto compatible = [](Vertex a, Vertex b) {
        return (a < 4 && b < 4) || (a >= 4 && a < 7 && b >= 4 && b < 7);
    };
    for (const Shape& shape : shapes) {
        const int n = shape.n;
        const int t = shape.t;
        const auto catalog = completeCatalog(n);
        std::vector<EdgeId> oldIds;
        std::vector<EdgeId> newIds;
        std::vector<Edge> oldEdges;
        for (std::size_t id = 0; id < catalog.size(); ++id) {
            const Edge edge = catalog[id];
            if (edge.v < t && !compatible(edge.u, edge.v)) {
                newIds.push_back(static_cast<EdgeId>(id));
            } else {
                oldIds.push_back(static_cast<EdgeId>(id));
                oldEdges.push_back(edge);
            }
        }
        IncrementalGraph graph(n, catalog);
        graph.applyBatch(oldIds);
        const std::size_t isolatedSize = static_cast<std::size_t>(n - t + 1);
        CliqueSet seeded;
        std::size_t candidates = 0U;
        for (const Clique& clique : referenceMaximalCliques(n, oldEdges)) {
            ++candidates;
            const bool isolated = clique.size() == isolatedSize && clique.front() >= 7;
            if (!isolated || clique.front() % 2 == 1) {
                seeded.insert(clique);
            }
        }
        require(candidates == static_cast<std::size_t>(t - 5),
                "pivot early exits: unexpected old maximal cliques");
        graph.applyBatch(newIds);
        Clique full(static_cast<std::size_t>(n));
        std::iota(full.begin(), full.end(), 0);
        const std::string where =
            "pivot early exits n=" + std::to_string(n) + " t=" + std::to_string(t);
        for (const Arm& arm : arms) {
            if ((arm.algorithm == SubsumedCliqueAlgorithm::Optimized ||
                 arm.algorithm == SubsumedCliqueAlgorithm::OptimizedBitmap) &&
                !shape.control) {
                continue;
            }
            MaximalCliqueIndex index(usesSubsumedSizeGate(arm.algorithm));
            for (const Clique& clique : seeded) {
                require(index.insert(clique), where + ": seed rejected");
            }
            SubsumedCliqueEnumerator enumerator(graph);
            SubsumedCliqueSearchStats stats;
            const std::vector<Clique> deleted = enumerator.enumerateBatch(
                newIds, std::array<Clique, 1U>{full}, index, stats,
                arm.algorithm);
            require(asSet(deleted, arm.name) == seeded,
                    where + ": " + arm.name + " deleted the wrong set");
            if (arm.algorithm == SubsumedCliqueAlgorithm::ResidualKernel) {
                // The shape reached the search it was built for.
                if (n <= 64) {
                    require(stats.residualKernelSearches == 1,
                            where + ": the single-word search did not run");
                } else if (t <= 64) {
                    require(stats.wide.kernelNarrow == 1,
                            where + ": runNarrow did not run");
                } else {
                    require(stats.wide.kernelWide == 1,
                            where + ": the multi-word search did not run");
                }
            }
        }
    }
}

// A kernel past the fixed widths' 256 endpoints.  The complement search takes
// it with the word count set at run time (`DynamicSearch`), with or without
// the size gate, and never splits it; the controls split it edge by edge in
// batch order.
//
// C = [0, 300) is complete in G + H.  The batch is the edge (298, 299) and
// then a star from vertex 0 to the leaves 268 down to 1, so |T| = 271.  Its
// maximal independent sets are {0} or [1, 269), each with 298 or 299: four
// candidates.  Batch order splits (298, 299) first, so the star runs on two
// candidates: 72,093 candidate occurrences.  The old maximal cliques are C
// minus {0} or minus [1, 269), each also minus 298 or 299; all but
// C \ {0, 299} are seeded, and that one has to miss.
void residualKernelWideKernelPastFixedWidths() {
    constexpr int n = 300;
    constexpr int leaves = 268;
    const auto catalog = completeCatalog(n);
    std::vector<EdgeId> oldIds;
    std::vector<EdgeId> star(static_cast<std::size_t>(leaves));
    EdgeId lastPair = invalidEdgeId;
    for (std::size_t id = 0; id < catalog.size(); ++id) {
        const Edge edge = catalog[id];
        const auto edgeId = static_cast<EdgeId>(id);
        if (edge.u == 0 && edge.v <= leaves) {
            star[static_cast<std::size_t>(leaves - edge.v)] = edgeId;
        } else if (edge.u == n - 2 && edge.v == n - 1) {
            lastPair = edgeId;
        } else {
            oldIds.push_back(edgeId);
        }
    }
    std::vector<EdgeId> newIds{lastPair};
    newIds.insert(newIds.end(), star.begin(), star.end());
    IncrementalGraph graph(n, catalog);
    graph.applyBatch(oldIds);
    graph.applyBatch(newIds);

    Clique full(static_cast<std::size_t>(n));
    std::iota(full.begin(), full.end(), 0);
    const auto without = [&](std::initializer_list<std::pair<Vertex, Vertex>> ranges) {
        Clique clique;
        for (Vertex v = 0; v < n; ++v) {
            bool dropped = false;
            for (const auto& [low, high] : ranges) {
                dropped = dropped || (v >= low && v < high);
            }
            if (!dropped) {
                clique.push_back(v);
            }
        }
        return clique;
    };
    const Clique hubAnd298 = without({{0, 1}, {n - 2, n - 1}});
    const Clique hubAnd299 = without({{0, 1}, {n - 1, n}});
    const Clique leavesAnd298 = without({{1, leaves + 1}, {n - 2, n - 1}});
    const Clique leavesAnd299 = without({{1, leaves + 1}, {n - 1, n}});
    const CliqueSet seeded{hubAnd298, leavesAnd298, leavesAnd299};

    struct Arm {
        const char* name;
        SubsumedCliqueAlgorithm algorithm;
        bool sizeHistogram;
    };
    constexpr auto rk = SubsumedCliqueAlgorithm::ResidualKernel;
    const std::array arms{
        Arm{"optimized", SubsumedCliqueAlgorithm::Optimized, false},
        Arm{"optimized-bitmap", SubsumedCliqueAlgorithm::OptimizedBitmap, false},
        Arm{"residual-kernel", rk, true},
        Arm{"residual-kernel/no-size-histogram", rk, false},
    };
    for (const Arm& arm : arms) {
        MaximalCliqueIndex index(
            usesSubsumedSizeGate(arm.algorithm) && arm.sizeHistogram);
        for (const Clique& clique : seeded) {
            require(index.insert(clique), "wide kernel seed");
        }
        require(!index.containsCanonical(hubAnd299), "wide kernel: unseeded clique present");
        SubsumedCliqueEnumerator enumerator(graph);
        SubsumedCliqueSearchStats stats;
        const std::vector<Clique> deleted = enumerator.enumerateBatch(
            newIds, std::array<Clique, 1U>{full}, index, stats, arm.algorithm);
        const std::string where = std::string("wide kernel: ") + arm.name;
        require(asSet(deleted, arm.name) == seeded, where + " deleted the wrong set");
        if (arm.algorithm == SubsumedCliqueAlgorithm::OptimizedBitmap) {
            require(stats.candidateOccurrences == 72093,
                    where + ": the batch-order split did not take the expected work");
        }
        if (arm.algorithm == rk) {
            require(stats.wide.kernelDynamic == 1,
                    where + ": the run-time-width search did not take the kernel");
            require(stats.splitRounds == 0 && stats.fallbackCliques == 0 &&
                        stats.bitmaskCliques == 1,
                    where + ": the wide clique was split");
            require(stats.kernelCandidates == 4,
                    where + ": the search did not emit the four independent sets");
        }
    }
}

// A kernel more than twice the fixed widths, through the enumerator: |T| =
// 520, nine words, so every multi-word step of `DynamicSearch` crosses word
// boundaries.  C = [0, 600) is complete in G + H, and every pair of T = [0,
// 520) is a batch edge except inside three groups that were already joined:
// [0, 10), [100, 105) and {300, 301, 302, 515, ..., 519} (across words).  The
// kernel's maximal independent sets are those three groups and every other
// vertex of T alone, so the old maximal cliques are K = [520, 600) with each
// of them -- 500 cliques, all subsumed by C.  Every one is seeded but those of
// the vertices 200..209, which have to miss.
void residualKernelDynamicWidthKernel() {
    constexpr int n = 600;
    constexpr int t = 520;
    const auto groupOf = [](Vertex v) {
        if (v < 10) return 1;
        if (v >= 100 && v < 105) return 2;
        if ((v >= 300 && v < 303) || (v >= 515 && v < t)) return 3;
        return 0;
    };
    const auto catalog = completeCatalog(n);
    std::vector<EdgeId> oldIds;
    std::vector<EdgeId> newIds;
    for (std::size_t id = 0; id < catalog.size(); ++id) {
        const Edge edge = catalog[id];
        const bool inT = edge.v < t;
        const bool grouped =
            groupOf(edge.u) != 0 && groupOf(edge.u) == groupOf(edge.v);
        (inT && !grouped ? newIds : oldIds).push_back(static_cast<EdgeId>(id));
    }
    IncrementalGraph graph(n, catalog);
    graph.applyBatch(oldIds);
    graph.applyBatch(newIds);

    Clique core;
    for (Vertex v = t; v < n; ++v) {
        core.push_back(v);
    }
    const auto withCore = [&](const std::vector<Vertex>& kept) {
        Clique clique = kept;
        clique.insert(clique.end(), core.begin(), core.end());
        std::ranges::sort(clique);
        return clique;
    };
    CliqueSet seeded;
    CliqueSet missing;
    for (int group = 1; group <= 3; ++group) {
        std::vector<Vertex> members;
        for (Vertex v = 0; v < t; ++v) {
            if (groupOf(v) == group) {
                members.push_back(v);
            }
        }
        seeded.insert(withCore(members));
    }
    for (Vertex v = 0; v < t; ++v) {
        if (groupOf(v) == 0) {
            (v >= 200 && v < 210 ? missing : seeded).insert(withCore({v}));
        }
    }
    require(seeded.size() + missing.size() == 500U,
            "dynamic-width kernel: wrong number of old cliques");
    Clique full(static_cast<std::size_t>(n));
    std::iota(full.begin(), full.end(), 0);
    for (const bool histogram : {true, false}) {
        const std::string where = std::string("dynamic-width kernel, size histogram ") +
            (histogram ? "on" : "off");
        MaximalCliqueIndex index(histogram);
        for (const Clique& clique : seeded) {
            require(index.insert(clique), where + ": seed rejected");
        }
        SubsumedCliqueEnumerator enumerator(graph);
        SubsumedCliqueSearchStats stats;
        const std::vector<Clique> deleted = enumerator.enumerateBatch(
            newIds, std::array<Clique, 1U>{full}, index, stats,
            SubsumedCliqueAlgorithm::ResidualKernel);
        require(asSet(deleted, where) == seeded, where + ": deleted the wrong set");
        require(stats.wide.kernelDynamic == 1,
                where + ": the run-time-width search did not take the kernel");
        require(stats.kernelCandidates == 500 && stats.splitRounds == 0,
                where + ": did not emit the 500 independent sets without a split");
    }
}

// `DynamicSearch` has to be `WideSearch` with the word count set at run time:
// on random kernels the fixed widths cover (65 to 256 endpoints, with and
// without a size bound), both emit the same leaves -- masks, sizes and hashes
// -- in the same order.
void dynamicSearchMatchesWideSearch() {
    struct Leaf {
        std::vector<std::uint64_t> mask;
        std::size_t count = 0;
        std::uint64_t hash = 0;
        bool operator==(const Leaf&) const = default;
    };
    std::mt19937 generator(20261007U);
    for (int trial = 0; trial < 120; ++trial) {
        const int t = std::uniform_int_distribution<int>(65, 256)(generator);
        const int coreCount = std::uniform_int_distribution<int>(0, 24)(generator);
        // Dense kernels: their maximal independent sets stay few.
        const double density =
            std::uniform_real_distribution<double>(0.85, 0.995)(generator);
        Clique clique(static_cast<std::size_t>(t + coreCount));
        std::iota(clique.begin(), clique.end(), 0);
        std::vector<Edge> edges;
        std::bernoulli_distribution keep(density);
        std::vector<int> touched(static_cast<std::size_t>(t), 0);
        for (Vertex u = 0; u < t; ++u) {
            for (Vertex v = u + 1; v < t; ++v) {
                if (keep(generator)) {
                    edges.push_back({u, v});
                    ++touched[static_cast<std::size_t>(u)];
                    ++touched[static_cast<std::size_t>(v)];
                }
            }
        }
        // T is the set of endpoints, so every vertex of [0, t) needs an edge.
        for (Vertex u = 0; u < t; ++u) {
            if (touched[static_cast<std::size_t>(u)] == 0) {
                const Vertex partner = u == 0 ? 1 : 0;
                edges.push_back(normalizedEdge(u, partner));
                ++touched[static_cast<std::size_t>(u)];
                ++touched[static_cast<std::size_t>(partner)];
            }
        }
        std::vector<Vertex> endpoints(static_cast<std::size_t>(t));
        std::iota(endpoints.begin(), endpoints.end(), 0);
        const std::size_t minimumSize = trial % 2 == 0
            ? 0U
            : static_cast<std::size_t>(coreCount) + 2U;
        std::vector<Leaf> fixed;
        std::vector<Leaf> dynamic;
        const auto into = [](std::vector<Leaf>& leaves) {
            return [&leaves](std::span<const std::uint64_t> mask, std::size_t count,
                             std::uint64_t hash) {
                leaves.push_back(Leaf{{mask.begin(), mask.end()}, count, hash});
            };
        };
        const auto fixedEmit = into(fixed);
        const std::size_t words = (static_cast<std::size_t>(t) + 63U) / 64U;
        wideresidual::Counters counters;
        if (words == 2U) {
            wideresidual::runWide<2U>(
                clique, endpoints, edges, minimumSize, counters, fixedEmit);
        } else if (words == 3U) {
            wideresidual::runWide<3U>(
                clique, endpoints, edges, minimumSize, counters, fixedEmit);
        } else {
            wideresidual::runWide<4U>(
                clique, endpoints, edges, minimumSize, counters, fixedEmit);
        }
        const auto dynamicEmit = into(dynamic);
        wideresidual::DynamicSearch<decltype(dynamicEmit)> search(
            clique, endpoints, edges, minimumSize, counters, dynamicEmit);
        search.run();
        require(!fixed.empty() || minimumSize != 0U,
                "dynamic search check: a kernel without independent sets");
        require(fixed == dynamic,
                "the run-time-width search differs from the fixed-width one at |T| = " +
                    std::to_string(t));
    }
}

} // namespace
} // namespace immce

int main() {
    try {
        immce::indexSurvivesRehashAndCompaction();
        immce::duplicateSubsumedCliqueIsEmittedOnce();
        immce::wideCliqueUsesCorrectRepresentation(64);
        immce::wideCliqueUsesCorrectRepresentation(65);
        immce::wideCliqueUsesCorrectRepresentation(257);
        immce::residualKernelRandomKernels();
        immce::residualKernelPivotEarlyExits();
        immce::residualKernelWideKernelPastFixedWidths();
        immce::residualKernelDynamicWidthKernel();
        immce::dynamicSearchMatchesWideSearch();
        immce::starWideDomainsMatchReference();
        immce::starNestedGraphsKeepTheirMasks();
        immce::exhaustiveFourVertices();
        immce::randomizedFiveVertices();
        immce::randomizedSixVertices();
        immce::doubleHubStates();
        immce::vertexRootsOwnOnlyTheirCliques();
        immce::vertexCoverSeedsLeaveOtherEndpointsInP();
        immce::removeBatchRejectsBadInput();
        immce::reverseReplayStreams();
        immce::exhaustiveMixedFourVertices();
        immce::randomizedMixedFiveVertices();
        immce::randomMixedStreams();
        std::cout << "all star and subsumed maximal-clique tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
