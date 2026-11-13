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
    // Cliques of at most 64 vertices with a single new edge, whose
    // candidates are C - u and C - v.
    long long oneEdgeKernels = 0;
    long long deferredProbes = 0;
    long long inlineProbes = 0;
    // Cliques whose new edges share one endpoint, and those whose new edges
    // join every pair of K_C (G[K_C] has no edge).
    long long sharedEndpointKernels = 0;
    long long edgelessKernels = 0;
    // The candidate subsumed cliques, but for those of `oneEdgeKernels`.
    long long kernelCandidates = 0;
    std::size_t maximumDeferredProbes = 0U;
    // The runs of `kernelEnumBkPivotNarrow` (cliques of at most 64 vertices) and
    // their nodes; `bkPivotWide` counts those of `KernelEnumBkPivotWide`.
    long long kernelSearches = 0;
    long long kernelSearchNodes = 0;
    KernelEnumBkPivotWideCounters bkPivotWide;
    // Deletion batches (`recoverBatch`): the candidates the destroyed cliques
    // gave, those below two vertices, those already accepted in this batch,
    // the maximality checks run on the rest, the candidates the checks
    // rejected and accepted, the neighbour entries and adjacency tests the
    // checks read, the centre groups checked on bit rows (`groups`), the
    // candidates checked that way (`groupedChecks`), the neighbour entries
    // the rows read (`rowReads`), the searches run under extension masks
    // (`maskedSearches`), the vertices extending a universal part
    // (`universalExtenders`), the masks kept (`extensionMasks`), the search
    // nodes a mask cut (`maskCuts`), and the destroyed cliques one vertex
    // extends whole (`extendedWhole`).
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

// Zeros over the vertices, from calloc: an array this large is mapped
// fresh, and its pages come into being as they are first touched, so a
// batch that reads a few vertices does not pay for clearing all of them.
class VertexScratch {
public:
    void ensure(std::size_t size);
    // Back to zeros, for an epoch counter that wrapped.
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

// KernelEnum (Section IV): the subsumed cliques of a batch, the old maximal
// cliques that a new maximal clique C strictly contains.  Each is a maximal
// clique of G[C], G being the graph before the batch, so KernelEnum
// enumerates those as candidates and keeps the ones the old maximal-clique
// index still holds.  With F_C the new edges inside C, the kernel K_C their
// endpoints and the universal part U_C = C \ K_C, the maximal cliques of
// G[C] are U_C u I over the maximal cliques I of G[K_C], which BK_pivot
// enumerates: two kernel vertices are adjacent in G unless a new edge joins
// them.  F_C is found among the batch endpoints in C, which a bitmap marked
// once per batch picks out by one bit test per vertex.  Three kernels need
// no search: one new edge (the candidates C - u and C - v), new edges
// joining every pair of K_C (G[K_C] has no edge: U_C with one kernel vertex
// each), and new edges sharing one endpoint c (U_C u {c} and
// U_C u (K_C \ {c})).
//
// Deletion is the same search run backwards (Section II-D).  The maximal
// cliques a deletion batch destroys are those holding a deleted edge -- the
// new-clique search on the graph before the deletion finds them -- and every
// maximal clique the deletion creates is a candidate of one of them.  What
// differs is the test: those candidates are not in any index yet, so
// `recoverBatch` keeps the ones that are maximal in the graph after the
// deletion.
class SubsumedSearch {
public:
    // The endpoints of one internal batch edge.  Carrying them instead of an
    // edge id keeps the whole per-clique path off the global edge catalog,
    // which on the larger graphs is hundreds of megabytes of cold memory.
    struct InternalBatchEdge {
        Vertex u = 0;
        Vertex v = 0;
    };

    explicit SubsumedSearch(const IncrementalGraph& graph);

    // The whole batch of new maximal cliques at once.  The subsumed cliques
    // come back canonical and each exactly once, in the order the search
    // found them, so a caller that compares or prints them sorts them itself.
    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        std::span<const Clique> newCliques,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats);

    // Clique by clique: `prepareBatch` once per batch, then
    // `consumeNewClique` for each new maximal clique as the Lambda_new search
    // emits it, then `finishBatch`.  Subsumed cliques are reported through
    // `output`.
    void prepareBatch(
        std::span<const EdgeId> batch,
        SubsumedSearchStats& stats);
    void consumeNewClique(
        const Clique& canonicalNewClique,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    // Closes the batch: issues every probe the batch held back.  Holding them
    // to the end of the batch is sound because no new maximal clique of a
    // batch is a proper subset of another, so the insertions that run between
    // them cannot meet a candidate.
    void finishBatch(
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);

    // Deletion: the maximal cliques the deletion of `deletedBatch` creates.
    // The graph must already be without the batch.  `destroyed` are the
    // maximal cliques the batch destroyed (the new-clique search on the graph
    // before the deletion).  Each of their candidates with at least two
    // vertices that is maximal in the graph now is inserted into `cliques`
    // and returned, once however many destroyed cliques produce it: no
    // candidate can equal a clique the index held before the batch (it is a
    // proper subset of one), so the index itself tells a repeat.
    [[nodiscard]] std::vector<Clique> recoverBatch(
        std::span<const EdgeId> deletedBatch,
        std::span<const Clique> destroyed,
        MaximalCliqueIndex& cliques,
        SubsumedSearchStats& stats);

    // Set hashes of the new cliques this batch has consumed, in order, so the
    // caller can insert them without hashing a second time.
    [[nodiscard]] std::span<const std::uint64_t> newCliqueHashes()
        const noexcept {
        return newCliqueHashes_;
    }

private:
    void validateBatchEdges() const;

    // The batch endpoints S in C, the new edges F_C among them, then the
    // candidates from the kernel K_C.
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
    // The candidates of a clique of at most 64 vertices from the rows of F_C
    // over clique positions (`kernel` marks K_C): a closed form, or
    // `kernelEnumBkPivotNarrow`.
    void enumerateCandidatesNarrow(
        const Clique& clique,
        const std::array<std::uint64_t, 64U>& newEdgeRows,
        std::uint64_t kernel,
        std::size_t edgeCount,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    // KernelEnum's BK_pivot on G[K_C] for a clique of at most 64 vertices,
    // every set one word on the stack (`KernelEnumBkPivotWide` is the same
    // search for the wider cliques).  In a recovery it also runs under the
    // extension masks.
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
    // The candidates of a clique of more than 64 vertices: the closed form of
    // new edges joining every pair of K_C, or `KernelEnumBkPivotWide`.
    void enumerateCandidatesWide(
        const Clique& clique,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);

    // The batch structures for the bitmap extraction, built on first use.
    void prepareBitmapExtraction();
    void ensureBitmapExtraction() {
        if (!bitmapExtractionPrepared_) {
            prepareBitmapExtraction();
        }
    }
    void releaseBitmapExtraction();
    // The slot of batch edge {u, v} in `batchEdgeSet_`, or `noBatchEdgeSlot`.
    static constexpr std::size_t noBatchEdgeSlot = ~std::size_t{0};
    [[nodiscard]] std::size_t findBatchEdge(Vertex u, Vertex v) const noexcept;
    [[nodiscard]] bool batchEdgeSetContains(Vertex u, Vertex v) const noexcept;
    // A clique with a single new edge: C minus either endpoint, the
    // endpoints at clique positions `uIndex < vIndex`.
    void probeOneEdgeKernel(
        const Clique& clique,
        std::size_t uIndex,
        std::size_t vIndex,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedSearchStats& stats,
        CliqueBatchOutput& output);
    // Buffers one candidate probe for the end of the batch, or issues it on
    // the spot when its mask takes more than one word.  `arenaOffset` is
    // where the producing clique's vertices were copied into `probeArena_`,
    // or `noArenaOffset` if they have not been copied yet.
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

    // Recovering: one candidate, as its clique under a mask or materialized.
    // What it accepts enters the index at the end of the batch.
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
    // Whether a clique of the graph as it is now has no vertex extending it.
    [[nodiscard]] bool maximalInGraph(
        std::span<const Vertex> candidate,
        SubsumedSearchStats& stats);
    // The candidate's adjacency lists into `memberLists_`, returning the
    // least degree; and the check on them.
    std::size_t readMemberLists(std::span<const Vertex> candidate);
    // `skipped` holds the candidate and is skipped in the scan: the
    // candidate itself, or the whole destroyed clique it came from.
    [[nodiscard]] bool maximalOnMemberLists(
        std::span<const Vertex> candidate,
        std::span<const Vertex> skipped,
        SubsumedSearchStats& stats);
    // Recovering: whether the batch already accepted an equal candidate,
    // and the record of one it accepts; every candidate the batch accepts
    // is in this table until the end of the batch, when they all enter the
    // index.
    [[nodiscard]] bool acceptedInBatch(
        std::span<const Vertex> candidate,
        std::uint64_t setHash) const noexcept;
    void recordAccepted(std::span<const Vertex> candidate, std::uint64_t setHash);
    void insertAcceptedInBatch(MaximalCliqueIndex& cliques);
    // KernelEnum's recovery: each vertex outside a destroyed clique D that
    // extends its universal part U_D (the clique positions `universal`), as
    // the kernel vertices it misses (bits over `kernel`, at most 64 of them),
    // handed to `onExtender`, which returns true to end the scan.
    template <typename OnExtender>
    void scanUniversalExtenders(
        const Clique& clique,
        std::span<const std::uint32_t> universal,
        std::span<const std::uint32_t> kernel,
        SubsumedSearchStats& stats,
        OnExtender&& onExtender);
    // The extension masks of a destroyed clique: into `extensionMasks_`, the
    // minimal ones only.  False when some vertex misses none: it extends
    // every candidate.
    [[nodiscard]] bool collectExtensionMasks(
        const Clique& clique,
        std::span<const std::uint32_t> universal,
        std::span<const std::uint32_t> kernel,
        SubsumedSearchStats& stats);

    // A recovery checks a candidate on the spot unless every member has a
    // high degree; those wait for the end of the batch, where each distinct
    // one is checked once.  Each then gets the member the most of them hold
    // -- as a rule a vertex that many destroyed cliques share -- as its
    // centre, and every vertex extending a candidate is a neighbour of its
    // centre.  A large group of one centre is checked on bit rows over the
    // centre's neighbourhood, each member's row read once for the whole
    // group, when that reads less than checking its candidates one by one;
    // the rest are checked one by one in the order they were made.  A batch
    // that gathers more than `gatherLimit` vertices has them checked then and
    // there, and gathers on: coauthor's window of batch 200000 gathers 5.2
    // billion in one batch.
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
    // The row of `vertex` over the current centre's neighbourhood, read on
    // first use within the group.
    [[nodiscard]] const std::uint64_t* groupRow(
        Vertex centre,
        Vertex vertex,
        std::size_t words,
        SubsumedSearchStats& stats);

    const IncrementalGraph& graph_;
    // Set for the length of `recoverBatch`.  A recovery checks a candidate
    // on the spot, by a galloping scan that stops at the first extender, or
    // at the end of the batch with the other deferred ones, grouped on bit
    // rows or one by one.  The duplicates are found in the batch's own table,
    // and the index is filled at the end of the batch.  A scan on the spot
    // skips the whole destroyed clique: the candidates are its maximal
    // subsets free of deleted edges, so a vertex of it left out of one cannot
    // extend that one.
    bool recovering_ = false;
    std::vector<Vertex> recoveredScratch_;
    // A maximality check's members by increasing degree: (degree, position
    // in the candidate), with their adjacency lists and the galloping
    // cursors into them.
    std::vector<std::pair<std::size_t, std::uint32_t>> memberOrder_;
    std::vector<std::span<const AdjacencyEntry>> memberLists_;
    std::vector<std::size_t> memberCursors_;
    // The destroyed clique a recovery is consuming (the checks on the spot
    // skip it).
    std::span<const Vertex> consumedClique_;
    // The candidates a recovering batch accepted, in an open-addressing
    // table of record numbers tagged with the batch.
    struct AcceptedRecord {
        std::uint64_t setHash = 0U;
        std::size_t offset = 0U;
        std::uint32_t size = 0U;
    };
    std::vector<AcceptedRecord> acceptedRecords_;
    std::vector<Vertex> acceptedVertices_;
    std::vector<std::uint64_t> acceptedSlots_;
    std::uint32_t acceptedEpoch_ = 0U;
    // Set while a search runs under extension masks: what it emits is
    // maximal, so the candidate is accepted unless a duplicate.
    bool knownMaximal_ = false;
    // Whether the batch accepted a candidate that a deferred one can equal
    // (under extension masks, or in an earlier round of the deferred ones),
    // so that those are looked up in the batch's table.
    bool filterGathered_ = false;
    // The lists the extension scan reads: those of U_D by increasing length,
    // and those of K_D, each with a galloping cursor.
    std::vector<std::pair<std::size_t, std::uint32_t>> universalOrder_;
    std::vector<std::span<const AdjacencyEntry>> universalLists_;
    std::vector<std::size_t> universalCursors_;
    std::vector<std::span<const AdjacencyEntry>> kernelLists_;
    std::vector<std::size_t> kernelCursors_;
    std::vector<std::uint64_t> extensionMasks_;
    // A clique's universal part and kernel as clique positions, and a wide
    // clique's new-edge degrees and candidate mask.
    std::vector<std::uint32_t> universalPositions_;
    std::vector<std::uint32_t> kernelPositions_;
    std::vector<std::uint32_t> newEdgeDegree_;
    std::vector<std::uint64_t> candidateMask_;
    // A recovering search's candidates, held (the kernel vertices kept and
    // their hash) while checking them one by one still looks cheaper than
    // the masks.
    std::vector<std::pair<std::uint64_t, std::uint64_t>> heldLeaves_;
    std::vector<Vertex> gatheredVertices_;
    std::vector<GatheredCandidate> gathered_;
    std::vector<std::uint32_t> gatheredOrder_;
    // The sort keys of the deferred candidates, radix-sorted by their high
    // half (with its scratch and counts): the hash's low half above the
    // position for the duplicates, then the centre above it for the groups;
    // and the positions kept, as a bitmap.
    std::vector<std::uint64_t> gatheredKeys_;
    std::vector<std::uint64_t> radixScratch_;
    std::vector<std::uint32_t> radixCounts_;
    std::vector<std::uint64_t> keptPositions_;
    // The deferred candidates by centre, and the groups among them that go
    // on rows, as ranges of that order.
    std::vector<std::uint32_t> centreOrder_;
    std::vector<std::pair<std::size_t, std::size_t>> rowGroups_;
    // The deferred candidates left to check one by one, in generation order.
    std::vector<std::uint32_t> oneByOneOrder_;
    // Per-vertex scratch, each valid while its epoch is current: how many
    // deferred candidates hold the vertex; its row in `rowArena_`.
    VertexScratch holdCount_;
    VertexScratch holdEpoch_;
    std::uint32_t holdEpochValue_ = 0U;
    // The centre's neighbourhood as a bit per vertex, with the position in
    // it of each word's first vertex, for the rows read by lookup: a
    // position is that plus the bits below in the word.  Set when a group
    // first needs it (`centreNumbered_`), cleared word by word after.
    std::vector<std::uint64_t> centreBits_;
    std::vector<std::uint32_t> centreRank_;
    // A row's positions read by lookup, before they are marked.
    std::vector<std::uint32_t> rowHits_;
    bool centreNumbered_ = false;
    VertexScratch rowOfVertex_;
    VertexScratch rowEpoch_;
    std::uint32_t rowEpochValue_ = 0U;
    std::vector<std::uint64_t> rowArena_;
    std::vector<std::uint64_t> rowScratch_;
    // The positions of the words a candidate's AND still has set.
    std::vector<std::uint32_t> liveWords_;
    std::vector<EdgeId> batchEdges_;
    bool batchPrepared_ = false;

    // The bitmap extraction's batch structures.  The bitmap has one bit per
    // vertex (a few hundred kilobytes on the largest graph, so it stays in
    // L2), and the edge set holds the batch edges as (min << 32 | max) keys in
    // an open-addressed table at most a quarter full.
    std::vector<std::uint64_t> endpointBits_;
    std::vector<Vertex> endpointList_;
    std::vector<std::uint64_t> batchEdgeSet_;
    std::vector<std::uint32_t> batchEdgeSlots_;
    unsigned batchEdgeSetShift_ = 64U;
    bool bitmapExtractionPrepared_ = false;
    // A clique's batch endpoints S as clique positions, its new edges F_C,
    // and their endpoints' clique positions.
    std::vector<std::uint32_t> batchEndpointPositions_;
    std::vector<InternalBatchEdge> newEdges_;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> newEdgePositions_;
    // A candidate is a proper subset of the new clique that produced it, and
    // no other new maximal clique of the same batch can equal it -- both are
    // maximal in G+H, so neither contains the other.  Every probe of a batch
    // is therefore held back to the end of it, which turns two dependent
    // random accesses per clique into one buffer that can be prefetched
    // arbitrarily deep.
    //
    // One candidate whose probe has not been issued yet.  The clique's
    // vertices live once in `probeArena_`; the candidate is that clique under
    // `mask`, which is why a held-back probe costs 24 bytes and not a set.
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
    // A deferred probe that reached a tag match, on its way to verification.
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

} // namespace ckm
