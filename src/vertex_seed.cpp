#include "vertex_seed.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ckm {
namespace {

// What a VertexSeed search carries down: whether the clique holds a batch
// edge yet. Once it does, every maximal clique below qualifies and nothing is
// checked again.
struct VertexSearchState {
    bool holdsBatchEdge = false;
};

} // namespace

VertexSeed::VertexSeed(const IncrementalGraph& graph)
    : graph_(graph),
      batchIndex_(graph),
      localGraph_(graph.vertexCount()) {
    const auto vertexCount = static_cast<std::size_t>(
        graph.vertexCount() < 0 ? 0 : graph.vertexCount());
    endpointRankByVertex_.assign(vertexCount, -1);
    epochByVertex_.assign(vertexCount, 0U);
}

std::vector<Clique> VertexSeed::enumerateBatch(
    std::span<const EdgeId> batch,
    VertexSeedStats& stats) {
    CliqueBatchOutput output(true);
    enumerateBatchInto(batch, stats, output);
    return std::move(output).takeCliques();
}

long long VertexSeed::countBatch(
    std::span<const EdgeId> batch,
    VertexSeedStats& stats) {
    CliqueBatchOutput output(false);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

long long VertexSeed::visitBatch(
    std::span<const EdgeId> batch,
    VertexSeedStats& stats,
    CliqueVisitor visitor) {
    if (!visitor) {
        throw std::invalid_argument("vertex-clique visitor must not be empty");
    }
    CliqueBatchOutput output(false, visitor);
    enumerateBatchInto(batch, stats, output);
    return output.count();
}

std::uint32_t VertexSeed::beginEpoch() {
    if (epoch_ == std::numeric_limits<std::uint32_t>::max()) {
        std::ranges::fill(epochByVertex_, 0U);
        epoch_ = 0U;
    }
    return ++epoch_;
}

void VertexSeed::enumerateBatchInto(
    std::span<const EdgeId> batch,
    VertexSeedStats& stats,
    CliqueBatchOutput& output) {
    if (batch.empty()) {
        return;
    }
    auto batchScope = batchIndex_.begin(batch);
    ++stats.batches;
    stats.batchEdges += static_cast<long long>(batch.size());
    indexBatch(batch);
    for (std::size_t rank = 0; rank < endpoints_.size(); ++rank) {
        enumerateRoot(endpoints_[rank], static_cast<int>(rank), stats, output);
    }
    for (Vertex endpoint : endpoints_) {
        endpointRankByVertex_[static_cast<std::size_t>(endpoint)] = -1;
    }
}

void VertexSeed::indexBatch(std::span<const EdgeId> batch) {
    endpoints_.clear();
    batchDegrees_.clear();
    for (const EdgeId edgeId : batch) {
        const Edge& edge = graph_.edge(edgeId);
        for (const Vertex endpoint : {edge.u, edge.v}) {
            int& rank = endpointRankByVertex_[static_cast<std::size_t>(endpoint)];
            if (rank < 0) {
                rank = static_cast<int>(endpoints_.size());
                endpoints_.push_back(endpoint);
                batchDegrees_.push_back(0U);
            }
            ++batchDegrees_[static_cast<std::size_t>(rank)];
        }
    }
    incidenceOffsets_.assign(endpoints_.size() + 1U, 0U);
    for (std::size_t index = 0; index < batchDegrees_.size(); ++index) {
        incidenceOffsets_[index + 1U] =
            incidenceOffsets_[index] + batchDegrees_[index];
    }
    incidenceNeighbors_.assign(incidenceOffsets_.back(), 0);
    incidenceFill_.assign(
        incidenceOffsets_.begin(), incidenceOffsets_.end() - 1);
    for (const EdgeId edgeId : batch) {
        const Edge& edge = graph_.edge(edgeId);
        const auto uRank = static_cast<std::size_t>(
            endpointRankByVertex_[static_cast<std::size_t>(edge.u)]);
        const auto vRank = static_cast<std::size_t>(
            endpointRankByVertex_[static_cast<std::size_t>(edge.v)]);
        const std::size_t uSlot = incidenceFill_[uRank]++;
        const std::size_t vSlot = incidenceFill_[vRank]++;
        incidenceNeighbors_[uSlot] = edge.v;
        incidenceNeighbors_[vSlot] = edge.u;
    }
}

std::span<const Vertex> VertexSeed::batchNeighbors(
    Vertex vertex) const noexcept {
    const int rank = endpointRankByVertex_[static_cast<std::size_t>(vertex)];
    if (rank < 0) {
        return {};
    }
    const auto index = static_cast<std::size_t>(rank);
    return std::span<const Vertex>(
        incidenceNeighbors_.data() + incidenceOffsets_[index],
        incidenceOffsets_[index + 1U] - incidenceOffsets_[index]);
}

void VertexSeed::enumerateRoot(
    Vertex root,
    int rootRank,
    VertexSeedStats& stats,
    CliqueBatchOutput& output) {
    ++stats.roots;
    // P is the whole neighbourhood but the earlier endpoints, which own every
    // clique holding them and so start in X.
    const auto startsInX = [&](int rank) {
        return rank >= 0 && rank < rootRank;
    };
    std::vector<Vertex> candidates;
    std::vector<Vertex> excluded;
    std::vector<Vertex> candidateEndpoints;
    const std::uint32_t epoch = beginEpoch();
    for (const AdjacencyEntry& entry : graph_.neighbors(root)) {
        const Vertex neighbor = entry.neighbor;
        const int rank =
            endpointRankByVertex_[static_cast<std::size_t>(neighbor)];
        if (startsInX(rank)) {
            excluded.push_back(neighbor);
            continue;
        }
        candidates.push_back(neighbor);
        epochByVertex_[static_cast<std::size_t>(neighbor)] = epoch;
        if (rank >= 0) {
            candidateEndpoints.push_back(neighbor);
        }
    }

    // The live test at the root, before its graph exists: one of the root's
    // own batch edges into P -- to a later endpoint -- or a batch edge between
    // two vertices of P.
    const auto rootReachesBatchEdge = [&] {
        for (const Vertex partner : batchNeighbors(root)) {
            ++stats.batchEdgeChecks;
            if (!startsInX(
                    endpointRankByVertex_[static_cast<std::size_t>(partner)])) {
                return true;
            }
        }
        for (const Vertex endpoint : candidateEndpoints) {
            for (const Vertex neighbor : batchNeighbors(endpoint)) {
                ++stats.batchEdgeChecks;
                if (epochByVertex_[static_cast<std::size_t>(neighbor)] ==
                    epoch) {
                    return true;
                }
            }
        }
        return false;
    };
    if (!rootReachesBatchEdge()) {
        ++stats.prunedRoots;
        return;
    }

    Clique current{root};
    localGraph_.rebuild(
        graph_, current, candidates, excluded, stats.localGraph);
    localGraph_.setBits(candidateEndpoints, endpointMask_);

    // A vertex entering the clique completes a batch edge exactly when it is
    // the batch neighbour of a member.
    const auto updateState = [&, this](
        const LocalSearchGraph&,
        VertexSearchState parent,
        std::span<const Vertex> clique,
        Vertex vertex) {
        if (!parent.holdsBatchEdge) {
            for (const Vertex neighbor : batchNeighbors(vertex)) {
                ++stats.batchEdgeChecks;
                if (std::ranges::find(clique, neighbor) != clique.end()) {
                    return VertexSearchState{true};
                }
            }
        }
        return parent;
    };
    const auto pruneNode = [&, this](
        const LocalSearchGraph& graph,
        std::span<const Vertex> clique,
        const LocalVertexSet& candidateSet,
        const LocalVertexSet&,
        const VertexSearchState& state) {
        // The root passed the same test before its graph was built.
        if (state.holdsBatchEdge || clique.size() == 1U) {
            return false;
        }
        if (batchEdgeReachable(graph, clique, candidateSet, stats)) {
            return false;
        }
        ++stats.liveTestPrunes;
        return true;
    };
    const auto onLeaf = [&](
        std::span<const Vertex> clique,
        const VertexSearchState& state) {
        if (!state.holdsBatchEdge) {
            throw std::logic_error("a VertexSeed leaf holds no batch edge");
        }
        ++stats.leaves;
        output.emit(clique);
        ++stats.emittedCliques;
    };
    bronKerboschLocal(
        localGraph_,
        current,
        localGraph_.makeSet(std::move(candidates)),
        localGraph_.makeSet(std::move(excluded)),
        VertexSearchState{},
        updateState,
        NoChildRefinement{},
        pruneNode,
        onLeaf,
        stats);
}

bool VertexSeed::batchEdgeReachable(
    const LocalSearchGraph& graph,
    std::span<const Vertex> clique,
    const LocalVertexSet& candidates,
    VertexSeedStats& stats) const {
    // A batch edge from the clique into P.
    for (const Vertex member : clique) {
        for (const Vertex neighbor : batchNeighbors(member)) {
            ++stats.batchEdgeChecks;
            if (graph.contains(candidates, neighbor)) {
                return true;
            }
        }
    }
    // A batch edge inside P, which only the endpoints in the root's P can
    // carry.
    const std::span<const Vertex> locals = graph.localVertices();
    for (std::size_t word = 0; word < candidates.bits.size(); ++word) {
        std::uint64_t bits = candidates.bits[word] & endpointMask_.bits[word];
        while (bits != 0U) {
            const Vertex endpoint = locals[word * 64U +
                static_cast<std::size_t>(std::countr_zero(bits))];
            bits &= bits - 1U;
            for (const Vertex neighbor : batchNeighbors(endpoint)) {
                ++stats.batchEdgeChecks;
                if (graph.contains(candidates, neighbor)) {
                    return true;
                }
            }
        }
    }
    return false;
}

} // namespace ckm
