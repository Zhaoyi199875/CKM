#pragma once

#include "batch_index.hpp"
#include "local_bk.hpp"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <span>
#include <vector>

namespace ckm {

struct CoverSeedConfiguration {
    bool reindex = true;
};

struct SeedIncidence {
    Vertex neighbor = 0;
    int ownerRank = 0;
    int order = 0;
};

struct CoverSeedStats {
    long long batches = 0;
    long long batchEdges = 0;
    long long seeds = 0;
    long long coverVertices = 0;
    long long coverBothEndpoints = 0;
    long long assignedEdges = 0;
    long long maximumSeedEdges = 0;
    long long singleEdgeSeeds = 0;
    long long singleEdgeSeedEdges = 0;
    long long anchorRootEdges = 0;
    long long pivotRootEdges = 0;
    long long reindexedTrees = 0;
    long long rootChildReindexedSubtrees = 0;
    long long deepReindexedSubtrees = 0;
    long long nodeWords = 0;
    long long pivotWords = 0;
    long long anchorRootSeeds = 0;
    long long pivotRootSeeds = 0;
    long long rootAnchorBranches = 0;
    long long rootPivotBranches = 0;
    long long rootChosenBranches = 0;
    long long rootDegreeDecisions = 0;
    long long rootPivotScans = 0;
    long long unanchoredAnchorChoices = 0;
    long long singlePartnerAnchors = 0;
    long long domainVertexOccurrences = 0;
    std::size_t maximumDomainVertices = 0U;
    double assignMilliseconds = 0.0;
    long long recursiveCalls = 0;
    long long pivotCalls = 0;
    long long pivotCandidatesEvaluated = 0;
    long long leaves = 0;
    long long ownershipMoves = 0;
    long long ownershipProbes = 0;
    long long unanchoredNodes = 0;
    long long liveTestPrunes = 0;
    long long emittedCliques = 0;
    LocalGraphStats localGraph;
};

class CoverSeed {
public:
    CoverSeed(
        const IncrementalGraph& graph,
        CoverSeedConfiguration configuration = {});

    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        CoverSeedStats& stats);

private:
    enum class BranchPartition {
        AnchorOrPivot,
        Pivot,
    };

    void enumerateBatchInto(
        std::span<const EdgeId> batch,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    void assignEdges(std::span<const EdgeId> batch, CoverSeedStats& stats);
    void buildIncidenceIndex(std::span<const EdgeId> batch);
    void buildDomain(
        Vertex seed,
        std::size_t seedIndex,
        CoverSeedStats& stats);
    void searchSingleEdgeSeed(
        Vertex seed,
        std::size_t seedIndex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    void searchSeed(
        Vertex seed,
        std::size_t seedIndex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    void prepareRoot(
        Vertex seed,
        int seedRank,
        CoverSeedStats& stats);
    void searchPartnerRoots(
        Vertex seed,
        std::size_t seedIndex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    void branchUnanchored(
        const LocalSearchGraph& graph,
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool atRoot,
        BranchPartition partition,
        Vertex knownPivot,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    void branchAnchored(
        const LocalSearchGraph& graph,
        Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        bool mayReindex,
        int seedRank,
        CoverSeedStats& stats,
        CliqueBatchOutput& output);
    template <typename Recurse>
    void descendInto(
        const LocalSearchGraph& graph,
        const Clique& current,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        CoverSeedStats& stats,
        const Recurse& recurse);
    [[nodiscard]] bool reindexPays(
        const LocalSearchGraph& graph,
        const Clique& current,
        const LocalVertexSet& candidates,
        const LocalVertexSet& excluded) const noexcept;
    [[nodiscard]] const LocalVertexSet& partnerMaskFor(
        const LocalSearchGraph& graph) const noexcept;
    void excludeClaimedNeighbors(
        const LocalSearchGraph& localGraph,
        Vertex vertex,
        int thresholdRank,
        int thresholdOrder,
        LocalVertexSet& candidates,
        LocalVertexSet& excluded,
        CoverSeedStats& stats) const;
    [[nodiscard]] Vertex partnerAt(std::size_t position, Vertex seed) const;
    [[nodiscard]] bool claimedEarlier(
        EdgeId edgeId,
        int thresholdRank,
        int thresholdOrder) const noexcept;
    [[nodiscard]] std::uint32_t beginEpoch();
    [[nodiscard]] std::uint32_t beginUnionEpoch();

    const IncrementalGraph& graph_;
    CoverSeedConfiguration configuration_;
    BatchIndex batchIndex_;
    LocalSearchGraph localGraph_;
    LocalSearchGraph treeGraph_;
    std::vector<Vertex> treeVertices_;
    std::vector<std::unique_ptr<LocalSearchGraph>> deepGraphs_;
    std::deque<LocalVertexSet> deepMasks_;
    std::deque<std::vector<Vertex>> deepVertices_;
    std::size_t deepLevel_ = 0U;
    bool maskActive_ = false;
    struct UnanchoredFrame {
        LocalVertexSet childCandidates;
        LocalVertexSet childExcluded;
        std::vector<Vertex> extensions;
    };
    std::deque<UnanchoredFrame> unanchoredFrames_;
    std::size_t unanchoredDepth_ = 0U;

    std::span<const EdgeId> batch_;
    std::vector<Vertex> seeds_;
    std::vector<std::size_t> seedOffsets_;
    std::vector<std::size_t> seedFill_;
    std::vector<std::size_t> ownedCounts_;
    std::vector<int> seedEdgePositions_;
    std::vector<int> ownerRankByPosition_;
    std::vector<int> orderByPosition_;
    std::vector<int> seedIndexByVertex_;
    std::vector<unsigned char> inCover_;
    std::vector<Vertex> endpoints_;
    std::vector<int> endpointRankByVertex_;
    std::vector<std::size_t> batchDegrees_;
    std::vector<std::size_t> incidenceOffsets_;
    std::vector<SeedIncidence> incidenceEntries_;
    std::vector<std::size_t> incidenceFill_;

    std::vector<Vertex> domain_;
    std::vector<Vertex> partners_;
    std::vector<Vertex> candidateVertices_;
    std::vector<Vertex> excludedVertices_;
    std::vector<std::uint32_t> epochByVertex_;
    std::uint32_t epoch_ = 0U;
    std::vector<std::uint32_t> unionEpochByVertex_;
    std::uint32_t unionEpoch_ = 0U;

    LocalVertexSet partnerMask_;
    LocalVertexSet allVertices_;
    LocalVertexSet rootCandidates_;
    LocalVertexSet rootExcluded_;
    LocalVertexSet finishedPartners_;
    LocalVertexSet partnerCandidates_;
    LocalVertexSet partnerExcluded_;
    Clique partnerClique_;
    LocalVertexSet compactCandidates_;
    LocalVertexSet compactExcluded_;
    LocalVertexSet treePartnerMask_;
};

}
