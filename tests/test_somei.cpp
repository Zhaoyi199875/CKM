#include "incremental_graph.hpp"
#include "somei.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// The SOMEi baseline against three references: brute force for the maximal
// cliques and for the v-header maximal cliques (Definition 2), and
// CM-Constructor (Yu & Liu 2017, Algorithms 1-2) run from scratch on the
// current graph, which fixes not only the bodies but their creation
// timestamps -- the state an incremental update must reproduce exactly.

namespace immce {
namespace {

using CliqueSet = std::set<Clique>;
// Per header: (body, creation vertex).
using HeaderState = std::map<Vertex, std::set<std::pair<Clique, Vertex>>>;

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

[[nodiscard]] EdgeId catalogId(int vertexCount, Vertex u, Vertex v) {
    const Edge wanted = normalizedEdge(u, v);
    EdgeId id = 0;
    for (Vertex a = 0; a < vertexCount; ++a) {
        for (Vertex b = a + 1; b < vertexCount; ++b) {
            if (a == wanted.u && b == wanted.v) {
                return id;
            }
            ++id;
        }
    }
    fail("edge outside the catalog");
}

// The graph as adjacency bit masks; vertexCount <= 16.
struct MaskGraph {
    int vertexCount = 0;
    std::vector<std::uint32_t> adjacency;

    MaskGraph(int count, const std::vector<Edge>& catalog,
              const std::vector<unsigned char>& active)
        : vertexCount(count), adjacency(static_cast<std::size_t>(count), 0U) {
        for (std::size_t id = 0; id < catalog.size(); ++id) {
            if (active[id] != 0U) {
                const Edge& edge = catalog[id];
                adjacency[static_cast<std::size_t>(edge.u)] |= 1U << edge.v;
                adjacency[static_cast<std::size_t>(edge.v)] |= 1U << edge.u;
            }
        }
    }

    [[nodiscard]] bool isClique(std::uint32_t mask) const {
        for (Vertex v = 0; v < vertexCount; ++v) {
            if ((mask >> v & 1U) != 0U &&
                (adjacency[static_cast<std::size_t>(v)] & mask) !=
                    (mask & ~(1U << v))) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] std::uint32_t commonNeighbours(std::uint32_t mask) const {
        std::uint32_t common = (1U << vertexCount) - 1U;
        for (Vertex v = 0; v < vertexCount; ++v) {
            if ((mask >> v & 1U) != 0U) {
                common &= adjacency[static_cast<std::size_t>(v)];
            }
        }
        return common & ~mask;
    }
};

[[nodiscard]] Clique members(std::uint32_t mask) {
    Clique result;
    for (Vertex v = 0; v < 32; ++v) {
        if ((mask >> v & 1U) != 0U) {
            result.push_back(v);
        }
    }
    return result;
}

// Maximal cliques with at least one edge, as the change stream counts them.
[[nodiscard]] CliqueSet bruteMaximalCliques(const MaskGraph& graph) {
    CliqueSet result;
    for (std::uint32_t mask = 1U; mask < (1U << graph.vertexCount); ++mask) {
        if (std::popcount(mask) >= 2 && graph.isClique(mask) &&
            graph.commonNeighbours(mask) == 0U) {
            result.insert(members(mask));
        }
    }
    return result;
}

// Definition 2: the v-header maximal cliques are the cliques with smallest
// vertex v and no common neighbour above v; as bodies (v left out).
[[nodiscard]] std::map<Vertex, CliqueSet> bruteHeaderBodies(
    const MaskGraph& graph) {
    std::map<Vertex, CliqueSet> result;
    for (std::uint32_t mask = 1U; mask < (1U << graph.vertexCount); ++mask) {
        if (std::popcount(mask) < 2 || !graph.isClique(mask)) {
            continue;
        }
        const Vertex header = std::countr_zero(mask);
        const std::uint32_t above = ~((2U << header) - 1U);
        if ((graph.commonNeighbours(mask) & above) != 0U) {
            continue;
        }
        Clique body = members(mask & ~(1U << header));
        result[header].insert(std::move(body));
    }
    return result;
}

// CM-Constructor, Algorithms 1 and 2 of Yu & Liu (2017), written out from
// that paper with its else-if chain: every vertex v_i in order, every small
// neighbour v_f of it as a header, R = L(N^s_{v_i}, v_f) n N^b_{v_f}.
[[nodiscard]] HeaderState cmConstructor(const MaskGraph& graph) {
    struct ListBody {
        Clique vertices;
        Vertex createdAt = 0;
    };
    std::vector<std::vector<ListBody>> bodies(
        static_cast<std::size_t>(graph.vertexCount));
    const auto small = [&](Vertex v) {
        return members(graph.adjacency[static_cast<std::size_t>(v)] &
                       ((1U << v) - 1U));
    };
    const auto big = [&](Vertex v) {
        return members(graph.adjacency[static_cast<std::size_t>(v)] &
                       ~((2U << v) - 1U));
    };
    for (Vertex vi = 0; vi < graph.vertexCount; ++vi) {
        const Clique smallNeighbours = small(vi);
        for (std::size_t position = 0; position < smallNeighbours.size();
             ++position) {
            const Vertex vf = smallNeighbours[position];
            const Clique bigOfHeader = big(vf);
            Clique renovation;
            std::set_intersection(
                smallNeighbours.begin() +
                    static_cast<std::ptrdiff_t>(position) + 1,
                smallNeighbours.end(), bigOfHeader.begin(), bigOfHeader.end(),
                std::back_inserter(renovation));
            std::vector<ListBody>& list = bodies[static_cast<std::size_t>(vf)];
            if (list.empty() || renovation.empty()) {
                list.push_back({{vi}, vi});
                continue;
            }
            bool perfectMatch = false;
            std::set<Clique> partTemp;
            std::vector<Clique> fullTemp;
            const std::size_t existing = list.size();
            for (std::size_t index = 0; index < existing; ++index) {
                Clique& body = list[index].vertices;
                Clique common;
                std::set_intersection(
                    body.begin(), body.end(), renovation.begin(),
                    renovation.end(), std::back_inserter(common));
                if (body == renovation) {
                    body.push_back(vi);
                    perfectMatch = true;
                    break;
                } else if (common == renovation) {
                    Clique created = renovation;
                    created.push_back(vi);
                    list.push_back({std::move(created), vi});
                    perfectMatch = true;
                    break;
                } else if (common == body) {
                    body.push_back(vi);
                    fullTemp.push_back(std::move(common));
                } else if (!common.empty()) {
                    partTemp.insert(std::move(common));
                }
            }
            if (perfectMatch) {
                continue;
            }
            for (const Clique& part : partTemp) {
                bool dropped = false;
                for (const Clique& other : partTemp) {
                    if (other != part && std::includes(
                            other.begin(), other.end(), part.begin(),
                            part.end())) {
                        dropped = true;
                    }
                }
                for (const Clique& other : fullTemp) {
                    if (std::includes(other.begin(), other.end(), part.begin(),
                                      part.end())) {
                        dropped = true;
                    }
                }
                if (!dropped) {
                    Clique created = part;
                    created.push_back(vi);
                    list.push_back({std::move(created), vi});
                }
            }
        }
    }
    HeaderState result;
    for (Vertex header = 0; header < graph.vertexCount; ++header) {
        for (const ListBody& body : bodies[static_cast<std::size_t>(header)]) {
            if (!result[header].insert({body.vertices, body.createdAt}).second) {
                fail("CM-Constructor produced a duplicate body");
            }
        }
    }
    return result;
}

[[nodiscard]] HeaderState someiState(const SomeiIndex& index, int vertexCount) {
    HeaderState result;
    for (Vertex header = 0; header < vertexCount; ++header) {
        for (const SomeiIndex::BodyView& body : index.bodies(header)) {
            Clique vertices(body.vertices.begin(), body.vertices.end());
            require(std::ranges::is_sorted(vertices) &&
                        (vertices.empty() || vertices.front() > header),
                    "SOMEi body is not an ascending list above its header");
            if (!result[header].insert({std::move(vertices), body.createdAt})
                     .second) {
                fail("SOMEi holds a duplicate body");
            }
        }
    }
    return result;
}

[[nodiscard]] CliqueSet asSet(const std::vector<Clique>& cliques,
                              const std::string& label) {
    CliqueSet result;
    for (const Clique& clique : cliques) {
        require(std::ranges::is_sorted(clique), label + " clique is not sorted");
        if (!result.insert(clique).second) {
            fail(label + " listed a clique twice");
        }
    }
    return result;
}

[[nodiscard]] CliqueSet difference(const CliqueSet& lhs, const CliqueSet& rhs) {
    CliqueSet result;
    std::ranges::set_difference(lhs, rhs, std::inserter(result, result.end()));
    return result;
}

// One engine on one edge stream: after every batch its change list, its
// ICTree and its candidate cliques are held to the references.
class Run {
public:
    Run(int vertexCount, SomeiConfiguration configuration, std::string label)
        : vertexCount_(vertexCount),
          catalog_(completeCatalog(vertexCount)),
          active_(catalog_.size(), 0U),
          graph_(vertexCount, catalog_),
          index_(graph_, configuration),
          label_(std::move(label)) {}

    void apply(const std::vector<EdgeId>& batch) {
        update(batch, {});
    }

    // A batch of the untimed build (MaintenanceEngine::build): an insertion
    // with SomeiIndex::setBuilding on, held to the same references.
    void build(const std::vector<EdgeId>& batch) {
        index_.setBuilding(true);
        update(batch, {});
        index_.setBuilding(false);
    }

    // A batch of deletions: the update scope is read while the edges are
    // still in the graph, the headers are redone after they are gone.
    void remove(const std::vector<EdgeId>& batch) {
        update({}, batch);
    }

    // A batch that inserts and deletes at once: the update scope is read off
    // the graph holding both (the union of the graphs before and after), the
    // headers are redone after the deletions.
    void change(
        const std::vector<EdgeId>& inserted,
        const std::vector<EdgeId>& deleted) {
        update(inserted, deleted);
    }

    [[nodiscard]] const SomeiIndex& index() const noexcept { return index_; }
    [[nodiscard]] EdgeId id(Vertex u, Vertex v) const {
        return catalogId(vertexCount_, u, v);
    }

private:
    void update(
        const std::vector<EdgeId>& inserted,
        const std::vector<EdgeId>& deleted) {
        const CliqueSet before =
            bruteMaximalCliques(MaskGraph(vertexCount_, catalog_, active_));
        for (const EdgeId id : inserted) {
            active_[static_cast<std::size_t>(id)] = 1U;
        }
        for (const EdgeId id : deleted) {
            active_[static_cast<std::size_t>(id)] = 0U;
        }
        const MaskGraph after(vertexCount_, catalog_, active_);
        const CliqueSet cliques = bruteMaximalCliques(after);

        std::vector<Clique> added;
        std::vector<Clique> removed;
        if (deleted.empty()) {
            graph_.applyBatch(inserted);
            index_.applyBatch(inserted, stats_, added, removed);
        } else {
            std::vector<EdgeId> changed(inserted);
            changed.insert(changed.end(), deleted.begin(), deleted.end());
            if (!inserted.empty()) {
                graph_.applyBatch(inserted);
            }
            index_.beginBatch(changed, stats_);
            graph_.removeBatch(deleted);
            index_.endBatch(stats_, added, removed);
        }

        if (index_.configuration().trackChanges) {
            require(asSet(added, label_) == difference(cliques, before),
                    label_ + ": new maximal cliques disagree with brute force");
            require(asSet(removed, label_) == difference(before, cliques),
                    label_ + ": deleted maximal cliques disagree with brute force");
        } else {
            require(added.empty() && removed.empty(),
                    label_ + ": change tracking off still listed changes");
        }
        require(asSet(index_.maximalCliques(), label_ + " ICTree") == cliques,
                label_ + ": ICTree branches disagree with brute force");
        require(index_.maximalCliqueCount() == cliques.size(),
                label_ + ": ICTree leaf count is off");

        const HeaderState state = someiState(index_, vertexCount_);
        require(state == cmConstructor(after),
                label_ + ": candidate cliques differ from CM-Constructor "
                         "run from scratch");
        const std::map<Vertex, CliqueSet> headerBodies = bruteHeaderBodies(after);
        std::size_t bodyCount = 0U;
        for (Vertex header = 0; header < vertexCount_; ++header) {
            CliqueSet bodies;
            const auto found = state.find(header);
            if (found != state.end()) {
                for (const auto& [body, createdAt] : found->second) {
                    bodies.insert(body);
                    require(body.front() <= createdAt &&
                                std::ranges::binary_search(body, createdAt),
                            label_ + ": a body does not hold its creation vertex");
                }
            }
            bodyCount += bodies.size();
            const auto expected = headerBodies.find(header);
            require(bodies == (expected == headerBodies.end() ? CliqueSet{}
                                                              : expected->second),
                    label_ + ": v-header maximal cliques disagree with "
                             "Definition 2");
        }
        require(index_.candidateCliqueCount() == bodyCount,
                label_ + ": candidate clique count is off");
    }

    int vertexCount_;
    std::vector<Edge> catalog_;
    std::vector<unsigned char> active_;
    IncrementalGraph graph_;
    SomeiIndex index_;
    SomeiStats stats_;
    std::string label_;
};

[[nodiscard]] std::vector<SomeiConfiguration> configurations() {
    std::vector<SomeiConfiguration> result;
    for (const bool pruning : {true, false}) {
        for (const bool trackChanges : {true, false}) {
            SomeiConfiguration configuration;
            configuration.pruning = pruning;
            configuration.trackChanges = trackChanges;
            configuration.profilePhases = !trackChanges;
            result.push_back(configuration);
        }
    }
    return result;
}

[[nodiscard]] std::string describe(const SomeiConfiguration& configuration) {
    return std::string("pruning=") + (configuration.pruning ? "on" : "off") +
        " changes=" + (configuration.trackChanges ? "on" : "off");
}

// Figs. 2-3 and Section V-B of the paper: the SOMEi of G, then the insertion
// e+(v3, v5) of Fig. 2(b).  Vertex v_k is k - 1 here.
void paperExample() {
    // G inserted as a measured batch, and as the untimed build whose steps
    // enumerate G[R]: the same SOMEi either way.
    for (const bool enumeratedBuild : {false, true})
    for (const bool pruning : {true, false}) {
        SomeiConfiguration configuration;
        configuration.pruning = pruning;
        Run run(7, configuration, std::string("paper example/") +
                                     (pruning ? "pruning" : "no pruning") +
                                     (enumeratedBuild ? "/enumerated build" : ""));
        const std::vector<std::pair<int, int>> paperEdges{
            {1, 2}, {1, 5}, {2, 3}, {2, 4}, {2, 5}, {2, 6}, {2, 7},
            {3, 4}, {3, 6}, {4, 5}, {4, 6}, {5, 6}, {6, 7}};
        std::vector<EdgeId> graphEdges;
        for (const auto& [u, v] : paperEdges) {
            graphEdges.push_back(run.id(u - 1, v - 1));
        }
        if (enumeratedBuild) {
            run.build(graphEdges);
        } else {
            run.apply(graphEdges);
        }

        // Fig. 3(a): eight v-header maximal cliques with their timestamps.
        const HeaderState figure3a{
            {0, {{{1, 4}, 1}}},
            {1, {{{2, 3, 5}, 2}, {{3, 4, 5}, 4}, {{5, 6}, 6}}},
            {2, {{{3, 5}, 3}}},
            {3, {{{4, 5}, 4}}},
            {4, {{{5}, 5}}},
            {5, {{{6}, 6}}}};
        require(someiState(run.index(), 7) == figure3a,
                "the SOMEi of G differs from Fig. 3(a)");
        require(run.index().maximalCliqueCount() == 4U,
                "G has four maximal cliques");

        run.apply({run.id(2, 4)});
        // Lemma 4's example: headers v1, v2, v3 roll back to T4 and redo
        // from v5; Lemma 8 prunes v1 (v3 is not in N^b_{v1}).
        std::vector<std::pair<Vertex, Vertex>> expectedScope{{1, 4}, {2, 4}};
        if (!pruning) {
            expectedScope.insert(expectedScope.begin(), {0, 4});
        }
        const auto scope = run.index().lastUpdateScope();
        require(std::vector<std::pair<Vertex, Vertex>>(scope.begin(), scope.end()) ==
                    expectedScope,
                "the update scope of e+(v3, v5) differs from the paper's");
        // Fig. 3(c).
        const HeaderState figure3c{
            {0, {{{1, 4}, 1}}},
            {1, {{{2, 3, 4, 5}, 2}, {{5, 6}, 6}}},
            {2, {{{3, 4, 5}, 3}}},
            {3, {{{4, 5}, 4}}},
            {4, {{{5}, 5}}},
            {5, {{{6}, 6}}}};
        require(someiState(run.index(), 7) == figure3c,
                "the SOMEi of G' differs from Fig. 3(c)");
        // Fig. 2(b): {v1,v2,v5}, {v2,v3,v4,v5,v6}, {v2,v6,v7}.
        require(asSet(run.index().maximalCliques(), "G'") ==
                    CliqueSet{{0, 1, 4}, {1, 2, 3, 4, 5}, {1, 5, 6}},
                "MC(G') differs from Fig. 2(b)");
    }
}

// Every assignment of absent / old / new to the edges of a small complete
// graph: the old edges as one batch, the new ones as the next.
void exhaustiveStates(int vertexCount) {
    const std::vector<Edge> catalog = completeCatalog(vertexCount);
    std::uint64_t assignments = 1U;
    for (std::size_t edge = 0; edge < catalog.size(); ++edge) {
        assignments *= 3U;
    }
    for (const SomeiConfiguration& configuration : configurations()) {
        for (std::uint64_t code = 0; code < assignments; ++code) {
            std::vector<EdgeId> oldEdges;
            std::vector<EdgeId> newEdges;
            std::uint64_t value = code;
            for (std::size_t edge = 0; edge < catalog.size(); ++edge) {
                if (value % 3U == 1U) {
                    oldEdges.push_back(static_cast<EdgeId>(edge));
                } else if (value % 3U == 2U) {
                    newEdges.push_back(static_cast<EdgeId>(edge));
                }
                value /= 3U;
            }
            Run run(vertexCount, configuration,
                    std::to_string(vertexCount) + " vertices/state " +
                        std::to_string(code) + "/" + describe(configuration));
            if (!oldEdges.empty()) {
                run.apply(oldEdges);
            }
            if (!newEdges.empty()) {
                run.apply(newEdges);
                // Read backwards: the new edges deleted again, from the
                // graph that holds both.
                run.remove(newEdges);
            }
        }
    }
}

// Every assignment of absent / kept / deleted / inserted to the edges of a
// small complete graph: the kept and the deleted edges as one batch, then one
// batch that deletes the ones and inserts the others.
void exhaustiveMixedStates(int vertexCount) {
    const std::vector<Edge> catalog = completeCatalog(vertexCount);
    std::uint64_t assignments = 1U;
    for (std::size_t edge = 0; edge < catalog.size(); ++edge) {
        assignments *= 4U;
    }
    for (const SomeiConfiguration& configuration : configurations()) {
        for (std::uint64_t code = 0; code < assignments; ++code) {
            std::vector<EdgeId> beforeEdges;
            std::vector<EdgeId> inserted;
            std::vector<EdgeId> deleted;
            std::uint64_t value = code;
            for (std::size_t edge = 0; edge < catalog.size(); ++edge) {
                const auto id = static_cast<EdgeId>(edge);
                if (value % 4U == 1U || value % 4U == 2U) {
                    beforeEdges.push_back(id);
                }
                if (value % 4U == 2U) {
                    deleted.push_back(id);
                } else if (value % 4U == 3U) {
                    inserted.push_back(id);
                }
                value /= 4U;
            }
            if (inserted.empty() || deleted.empty()) {
                continue;
            }
            Run run(vertexCount, configuration,
                    std::to_string(vertexCount) + " vertices/mixed state " +
                        std::to_string(code) + "/" + describe(configuration));
            run.apply(beforeEdges);
            run.change(inserted, deleted);
        }
    }
}

// Random graphs built in one batch, then a run of batches that each delete
// some present edges and insert some absent ones, deleted edges coming back
// too, and finally everything deleted.
void randomMixedStreams() {
    std::mt19937 generator(20261013U);
    for (int trial = 0; trial < 300; ++trial) {
        const int vertexCount = std::uniform_int_distribution<int>(5, 9)(generator);
        const std::size_t catalogSize =
            static_cast<std::size_t>(vertexCount * (vertexCount - 1) / 2);
        std::bernoulli_distribution initially(
            std::uniform_real_distribution<double>(0.2, 0.9)(generator));
        std::bernoulli_distribution dropping(
            std::uniform_real_distribution<double>(0.1, 0.6)(generator));
        std::bernoulli_distribution adding(
            std::uniform_real_distribution<double>(0.1, 0.6)(generator));
        std::vector<unsigned char> present(catalogSize, 0U);
        std::vector<EdgeId> initial;
        for (std::size_t id = 0; id < catalogSize; ++id) {
            if (initially(generator)) {
                initial.push_back(static_cast<EdgeId>(id));
                present[id] = 1U;
            }
        }
        std::vector<std::pair<std::vector<EdgeId>, std::vector<EdgeId>>> changes;
        const int batchCount = std::uniform_int_distribution<int>(3, 8)(generator);
        for (int batch = 0; batch < batchCount; ++batch) {
            std::vector<EdgeId> inserted;
            std::vector<EdgeId> deleted;
            for (std::size_t id = 0; id < catalogSize; ++id) {
                if (present[id] != 0U ? dropping(generator) : adding(generator)) {
                    (present[id] != 0U ? deleted : inserted)
                        .push_back(static_cast<EdgeId>(id));
                }
            }
            std::ranges::shuffle(inserted, generator);
            std::ranges::shuffle(deleted, generator);
            for (const EdgeId id : inserted) {
                present[static_cast<std::size_t>(id)] = 1U;
            }
            for (const EdgeId id : deleted) {
                present[static_cast<std::size_t>(id)] = 0U;
            }
            changes.emplace_back(std::move(inserted), std::move(deleted));
        }
        std::vector<EdgeId> remaining;
        for (std::size_t id = 0; id < catalogSize; ++id) {
            if (present[id] != 0U) {
                remaining.push_back(static_cast<EdgeId>(id));
            }
        }
        for (const SomeiConfiguration& configuration : configurations()) {
            Run run(vertexCount, configuration,
                    "mixed trial " + std::to_string(trial) + "/" +
                        describe(configuration));
            if (!initial.empty()) {
                run.apply(initial);
            }
            for (const auto& [inserted, deleted] : changes) {
                if (!inserted.empty() || !deleted.empty()) {
                    run.change(inserted, deleted);
                }
            }
            if (!remaining.empty()) {
                run.remove(remaining);
            }
            require(run.index().maximalCliqueCount() == 0U &&
                        run.index().candidateCliqueCount() == 0U,
                    "deleting every edge left SOMEi holding cliques");
        }
    }
}

// Random graphs fed as random edge orders cut into random batches.
void randomStreams() {
    std::mt19937 generator(20261002U);
    for (int trial = 0; trial < 400; ++trial) {
        const int vertexCount = std::uniform_int_distribution<int>(5, 10)(generator);
        const double density =
            std::uniform_real_distribution<double>(0.25, 0.95)(generator);
        const std::size_t largestBatch =
            std::uniform_int_distribution<std::size_t>(1U, 12U)(generator);
        std::vector<EdgeId> stream;
        std::bernoulli_distribution keep(density);
        const std::size_t catalogSize =
            static_cast<std::size_t>(vertexCount * (vertexCount - 1) / 2);
        for (std::size_t id = 0; id < catalogSize; ++id) {
            if (keep(generator)) {
                stream.push_back(static_cast<EdgeId>(id));
            }
        }
        std::ranges::shuffle(stream, generator);
        std::vector<std::vector<EdgeId>> batches;
        for (std::size_t begin = 0; begin < stream.size();) {
            const std::size_t size = std::min(
                stream.size() - begin,
                std::uniform_int_distribution<std::size_t>(1U, largestBatch)(
                    generator));
            batches.emplace_back(
                stream.begin() + static_cast<std::ptrdiff_t>(begin),
                stream.begin() + static_cast<std::ptrdiff_t>(begin + size));
            begin += size;
        }
        for (const SomeiConfiguration& configuration : configurations()) {
            Run run(vertexCount, configuration,
                    "random trial " + std::to_string(trial) + "/" +
                        describe(configuration));
            for (const std::vector<EdgeId>& batch : batches) {
                run.apply(batch);
            }
            // And back to the empty graph, the last batch first.
            for (auto batch = batches.rbegin(); batch != batches.rend(); ++batch) {
                run.remove(*batch);
            }
            require(run.index().maximalCliqueCount() == 0U &&
                        run.index().candidateCliqueCount() == 0U,
                    "deleting every batch left SOMEi holding cliques");
        }
    }
}

// Two SOMEi in the very same state: every header's bodies in the same order
// with the same creation vertices, and the same ICTree node by node.
void requireSameState(
    const SomeiIndex& expected,
    const SomeiIndex& actual,
    int vertexCount,
    const std::string& label) {
    for (Vertex header = 0; header < vertexCount; ++header) {
        const std::vector<SomeiIndex::BodyView> lhs = expected.bodies(header);
        const std::vector<SomeiIndex::BodyView> rhs = actual.bodies(header);
        require(lhs.size() == rhs.size(), label + ": a header has other bodies");
        for (std::size_t index = 0; index < lhs.size(); ++index) {
            require(std::ranges::equal(lhs[index].vertices, rhs[index].vertices) &&
                        lhs[index].createdAt == rhs[index].createdAt,
                    label + ": bodies differ or come in another order");
        }
    }
    require(expected.maximalCliques() == actual.maximalCliques() &&
                expected.treeNodeCount() == actual.treeNodeCount() &&
                expected.maximalCliqueCount() == actual.maximalCliqueCount() &&
                expected.candidateCliqueCount() == actual.candidateCliqueCount(),
            label + ": the ICTrees differ");
}

// The untimed build by enumeration (SomeiConfiguration::enumerateBuild)
// against the paper's construction, on random graphs built in one batch or
// a few: after every build batch the same state, both held to the
// references, and the same again after the paper's update batches (mixed
// insertions and deletions) that follow.  A small enumerateBuildLimit sends
// the steps with a larger R back to the paper's lines.
void enumeratedBuilds() {
    std::mt19937 generator(20261106U);
    for (int trial = 0; trial < 400; ++trial) {
        const int vertexCount = std::uniform_int_distribution<int>(5, 12)(generator);
        const std::size_t catalogSize =
            static_cast<std::size_t>(vertexCount * (vertexCount - 1) / 2);
        std::bernoulli_distribution keep(
            std::uniform_real_distribution<double>(0.25, 0.95)(generator));
        std::vector<EdgeId> stream;
        std::vector<unsigned char> present(catalogSize, 0U);
        for (std::size_t id = 0; id < catalogSize; ++id) {
            if (keep(generator)) {
                stream.push_back(static_cast<EdgeId>(id));
                present[id] = 1U;
            }
        }
        if (stream.empty()) {
            continue;
        }
        std::ranges::shuffle(stream, generator);
        std::vector<std::vector<EdgeId>> builds;
        const int buildCount = std::uniform_int_distribution<int>(1, 3)(generator);
        for (int part = 0; part < buildCount; ++part) {
            const std::size_t begin = stream.size() * static_cast<std::size_t>(part) /
                static_cast<std::size_t>(buildCount);
            const std::size_t end = stream.size() * static_cast<std::size_t>(part + 1) /
                static_cast<std::size_t>(buildCount);
            if (end > begin) {
                builds.emplace_back(
                    stream.begin() + static_cast<std::ptrdiff_t>(begin),
                    stream.begin() + static_cast<std::ptrdiff_t>(end));
            }
        }
        std::bernoulli_distribution flip(0.3);
        std::vector<std::pair<std::vector<EdgeId>, std::vector<EdgeId>>> changes;
        const int batchCount = std::uniform_int_distribution<int>(1, 4)(generator);
        for (int batch = 0; batch < batchCount; ++batch) {
            std::vector<EdgeId> inserted;
            std::vector<EdgeId> deleted;
            for (std::size_t id = 0; id < catalogSize; ++id) {
                if (flip(generator)) {
                    (present[id] != 0U ? deleted : inserted)
                        .push_back(static_cast<EdgeId>(id));
                    present[id] ^= 1U;
                }
            }
            changes.emplace_back(std::move(inserted), std::move(deleted));
        }
        const std::size_t limit = std::uniform_int_distribution<int>(0, 2)(generator) == 0
            ? std::uniform_int_distribution<std::size_t>(1U, 3U)(generator)
            : std::size_t{4096};
        for (const bool pruning : {true, false}) {
            SomeiConfiguration paper;
            paper.pruning = pruning;
            paper.enumerateBuild = false;
            SomeiConfiguration enumerated = paper;
            enumerated.enumerateBuild = true;
            enumerated.enumerateBuildLimit = limit;
            const std::string label = "enumerated build trial " +
                std::to_string(trial) + "/" + describe(paper) + "/limit " +
                std::to_string(limit);
            Run expected(vertexCount, paper, label + "/paper");
            Run actual(vertexCount, enumerated, label + "/enumerated");
            for (const std::vector<EdgeId>& batch : builds) {
                expected.build(batch);
                actual.build(batch);
                requireSameState(expected.index(), actual.index(), vertexCount,
                                 label + " after a build batch");
            }
            for (const auto& [inserted, deleted] : changes) {
                if (inserted.empty() && deleted.empty()) {
                    continue;
                }
                expected.change(inserted, deleted);
                actual.change(inserted, deleted);
                requireSameState(expected.index(), actual.index(), vertexCount,
                                 label + " after an update batch");
            }
        }
    }
}

} // namespace
} // namespace immce

int main() {
    try {
        immce::paperExample();
        immce::exhaustiveStates(4);
        immce::exhaustiveStates(5);
        immce::randomStreams();
        immce::exhaustiveMixedStates(4);
        immce::randomMixedStreams();
        immce::enumeratedBuilds();
        std::cout << "all SOMEi tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "test failure: " << error.what() << '\n';
        return 1;
    }
}
