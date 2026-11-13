#pragma once

#include "incremental_graph.hpp"
#include "maximal_clique_index.hpp"
#include "kernel_enum_bk_pivot.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace ckm {

struct SubsumedSearchStats {
    long long batches = 0;
    long long newCliquesProcessed = 0;
    long long internalBatchEdges = 0;
    long long indexLookups = 0;
    long long emittedCliques = 0;
    long long batchEdgeVisits = 0;
    long long sizeGateSkippedCliques = 0;
    long long sizeGateSkippedProbes = 0;
    long long oneEdgeKernels = 0;
    long long deferredProbes = 0;
    long long inlineProbes = 0;
    long long sharedEndpointKernels = 0;
    long long edgelessKernels = 0;
    long long kernelCandidates = 0;
    std::size_t maximumDeferredProbes = 0U;
    long long kernelSearches = 0;
    long long kernelSearchNodes = 0;
    KernelEnumBkPivotWideCounters bkPivotWide;
    struct Recovery {
        long long candidates = 0;
        long long belowTwo = 0;
        long long duplicates = 0;
        long long checks = 0;
        long long rejected = 0;
        long long accepted = 0;
        long long neighborReads = 0;
        long long adjacencyTests = 0;
        long long groups = 0;
        long long groupedChecks = 0;
        long long rowReads = 0;
        long long maskedSearches = 0;
        long long universalExtenders = 0;
        long long extensionMasks = 0;
        long long maskCuts = 0;
        long long extendedWhole = 0;
    };
    Recovery recovery;
};

class VertexScratch {
public:
    void ensure(std::size_t size);
    void clear() noexcept;
    [[nodiscard]] std::uint32_t& operator[](std::size_t index) noexcept {
        return data_.get()[index];
    }

private:
    struct Free {
        void operator()(std::uint32_t* pointer) const noexcept {
            std::free(pointer);
        }
    };
    std::unique_ptr<std::uint32_t, Free> data_;
    std::size_t size_ = 0U;
};

class SubsumedSearch {
public:
    struct InternalBatchEdge {
        Vertex u = 0;
        Vertex v = 0;
    };

    explicit SubsumedSearch(const IncrementalGraph& graph);

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        std::span<const Clique> newCliques,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats);

    void prepareBatch(
        std::span<const EdgeId> batch,
        SubsumedSearchStats& stats);
    void consumeNewClique(
        const Clique& canonicalNewClique,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void finishBatch(
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);

    [[nodiscard]] std::vector<Clique> recoverBatch(
        std::span<const EdgeId> deletedBatch,
        std::span<const Clique> destroyed,
        MaximalCliqueIndex& cliques,
        SubsumedSearchStats& stats);

    [[nodiscard]] std::span<const std::uint64_t> newCliqueHashes()
        const noexcept {
        return newCliqueHashes_;
    }

private:
    void validateBatchEdges() const;

    void consumeNewCliqueKernelEnum(
        const Clique& clique,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void findNewEdgesInside(
        const Clique& clique,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void enumerateCandidates(
        const Clique& clique,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void enumerateCandidatesNarrow(
        const Clique& clique,
        const std::array<std::uint64_t, 64U>& newEdgeRows,
        std::uint64_t kernel,
        std::size_t edgeCount,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void kernelEnumBkPivotNarrow(
        const Clique& clique,
        const std::array<std::uint64_t, 64U>& newEdgeRows,
        std::uint64_t kernel,
        std::uint64_t universal,
        std::uint64_t kernelHash,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void enumerateCandidatesWide(
        const Clique& clique,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);

    void prepareBitmapExtraction();
    void ensureBitmapExtraction() {
        if (!bitmapExtractionPrepared_) {
            prepareBitmapExtraction();
        }
    }
    void releaseBitmapExtraction();
    static constexpr std::size_t noBatchEdgeSlot = ~std::size_t{0};
    [[nodiscard]] std::size_t findBatchEdge(Vertex u, Vertex v) const noexcept;
    [[nodiscard]] bool batchEdgeSetContains(Vertex u, Vertex v) const noexcept;
    void probeOneEdgeKernel(
        const Clique& clique,
        std::size_t uIndex,
        std::size_t vIndex,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void submitProbe(
        std::span<const Vertex> clique,
        std::span<const std::uint64_t> maskWords,
        std::size_t memberCount,
        std::uint64_t setHash,
        std::uint32_t& arenaOffset,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    static constexpr std::uint32_t noArenaOffset = 0xFFFFFFFFU;
    [[nodiscard]] std::size_t minimumCandidateSize(
        const MaximalCliqueIndex& cliques) const noexcept {
        return recovering_ ? 2U : cliques.minimumLiveCliqueSize();
    }
    [[nodiscard]] bool candidateSizePossible(
        const MaximalCliqueIndex& cliques,
        std::size_t size) const noexcept {
        return recovering_ ? size >= 2U : cliques.hasLiveCliqueOfSize(size);
    }

    void acceptRecovered(
        std::span<const Vertex> clique,
        std::span<const std::uint64_t> maskWords,
        std::size_t memberCount,
        std::uint64_t setHash,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void acceptRecoveredClique(
        std::span<const Vertex> candidate,
        std::uint64_t setHash,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    [[nodiscard]] bool maximalInGraph(
        std::span<const Vertex> candidate,
        SubsumedSearchStats& stats);
    std::size_t readMemberLists(std::span<const Vertex> candidate);
    [[nodiscard]] bool maximalOnMemberLists(
        std::span<const Vertex> candidate,
        std::span<const Vertex> skipped,
        SubsumedSearchStats& stats);
    [[nodiscard]] bool acceptedInBatch(
        std::span<const Vertex> candidate,
        std::uint64_t setHash) const noexcept;
    void recordAccepted(std::span<const Vertex> candidate, std::uint64_t setHash);
    void insertAcceptedInBatch(MaximalCliqueIndex& cliques);
    template <typename OnExtender>
    void scanUniversalExtenders(
        const Clique& clique,
        std::span<const std::uint32_t> universal,
        std::span<const std::uint32_t> kernel,
        SubsumedSearchStats& stats,
        OnExtender&& onExtender);
    [[nodiscard]] bool collectExtensionMasks(
        const Clique& clique,
        std::span<const std::uint32_t> universal,
        std::span<const std::uint32_t> kernel,
        SubsumedSearchStats& stats);

    struct GatheredCandidate {
        std::uint64_t setHash = 0U;
        std::size_t offset = 0U;
        std::uint32_t size = 0U;
        Vertex centre = 0;
        bool onRows = false;
    };
    void checkGatheredCandidates(
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    void sortKeysByHighHalf();
    [[nodiscard]] std::span<const Vertex> gatheredVertices(
        const GatheredCandidate& candidate) const noexcept {
        return std::span<const Vertex>(
            gatheredVertices_.data() + candidate.offset, candidate.size);
    }
    [[nodiscard]] const std::uint64_t* groupRow(
        Vertex centre,
        Vertex vertex,
        std::size_t words,
        SubsumedSearchStats& stats);

    const IncrementalGraph& graph_;
    bool recovering_ = false;
    std::vector<Vertex> recoveredScratch_;
    std::vector<std::pair<std::size_t, std::uint32_t>> memberOrder_;
    std::vector<std::span<const AdjacencyEntry>> memberLists_;
    std::vector<std::size_t> memberCursors_;
    std::span<const Vertex> consumedClique_;
    struct AcceptedRecord {
        std::uint64_t setHash = 0U;
        std::size_t offset = 0U;
        std::uint32_t size = 0U;
    };
    std::vector<AcceptedRecord> acceptedRecords_;
    std::vector<Vertex> acceptedVertices_;
    std::vector<std::uint64_t> acceptedSlots_;
    std::uint32_t acceptedEpoch_ = 0U;
    bool knownMaximal_ = false;
    bool filterGathered_ = false;
    std::vector<std::pair<std::size_t, std::uint32_t>> universalOrder_;
    std::vector<std::span<const AdjacencyEntry>> universalLists_;
    std::vector<std::size_t> universalCursors_;
    std::vector<std::span<const AdjacencyEntry>> kernelLists_;
    std::vector<std::size_t> kernelCursors_;
    std::vector<std::uint64_t> extensionMasks_;
    std::vector<std::uint32_t> universalPositions_;
    std::vector<std::uint32_t> kernelPositions_;
    std::vector<std::uint32_t> newEdgeDegree_;
    std::vector<std::uint64_t> candidateMask_;
    std::vector<std::pair<std::uint64_t, std::uint64_t>> heldLeaves_;
    std::vector<Vertex> gatheredVertices_;
    std::vector<GatheredCandidate> gathered_;
    std::vector<std::uint32_t> gatheredOrder_;
    std::vector<std::uint64_t> gatheredKeys_;
    std::vector<std::uint64_t> radixScratch_;
    std::vector<std::uint32_t> radixCounts_;
    std::vector<std::uint64_t> keptPositions_;
    std::vector<std::uint32_t> centreOrder_;
    std::vector<std::pair<std::size_t, std::size_t>> rowGroups_;
    std::vector<std::uint32_t> oneByOneOrder_;
    VertexScratch holdCount_;
    VertexScratch holdEpoch_;
    std::uint32_t holdEpochValue_ = 0U;
    std::vector<std::uint64_t> centreBits_;
    std::vector<std::uint32_t> centreRank_;
    std::vector<std::uint32_t> rowHits_;
    bool centreNumbered_ = false;
    VertexScratch rowOfVertex_;
    VertexScratch rowEpoch_;
    std::uint32_t rowEpochValue_ = 0U;
    std::vector<std::uint64_t> rowArena_;
    std::vector<std::uint64_t> rowScratch_;
    std::vector<std::uint32_t> liveWords_;
    std::vector<EdgeId> batchEdges_;
    bool batchPrepared_ = false;

    std::vector<std::uint64_t> endpointBits_;
    std::vector<Vertex> endpointList_;
    std::vector<std::uint64_t> batchEdgeSet_;
    std::vector<std::uint32_t> batchEdgeSlots_;
    unsigned batchEdgeSetShift_ = 64U;
    bool bitmapExtractionPrepared_ = false;
    std::vector<std::uint32_t> batchEndpointPositions_;
    std::vector<InternalBatchEdge> newEdges_;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> newEdgePositions_;
    struct DeferredProbe {
        std::uint64_t setHash = 0U;
        std::uint64_t mask = 0U;
        std::uint32_t offset = 0U;
        std::uint16_t size = 0U;
        std::uint16_t memberCount = 0U;
    };
    void flushDeferredProbes(
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    struct ResolvedProbe {
        std::size_t slot = 0U;
        std::uint32_t probeIndex = 0U;
        std::uint32_t record = 0U;
    };
    std::vector<DeferredProbe> probeBuffer_;
    std::vector<ResolvedProbe> resolvedProbes_;
    std::vector<Vertex> probeArena_;
    std::vector<std::uint64_t> newCliqueHashes_;
};

}
