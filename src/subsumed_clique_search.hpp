#pragma once

#include "incremental_graph.hpp"
#include "maximal_clique_index.hpp"
#include "wide_residual.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace immce {

// The subsumed clique searches this binary can run.
//
// `Optimized` is the control: the paper's split-and-probe search, with the
// batch-edge anchor index and the fixed-width bitmask candidate set, and
// nothing else.  The other two are measured against it under the same
// Lambda_new producer.
//
// Adding a framework means adding a value here, a name in
// `subsumedCliqueAlgorithmName`, and a branch in `consumeNewCliqueResolved`.
// Both pipelines dispatch through that one function, so the implementation is
// written once.
enum class SubsumedCliqueAlgorithm {
    Optimized,
    // Residual kernel search.  The candidates for a new clique C are K union
    // a maximal independent set of its conflict kernel F (T = V(F), K = C \ T),
    // with the kernel itself found output-sensitively: batch endpoints are
    // marked in a bitmap once per batch, so S = C n U costs one bit test per
    // vertex of C, and the kernel is the batch edges inside S.  |S| = 2 is the
    // one-edge kernel and needs no lookup at all.  One-edge, star and complete
    // kernels have closed forms; every other kernel goes to the complement
    // Bron-Kerbosch, whatever its density, edge count or width, so a kernel is
    // never split.  The index's live clique-size histogram gates cliques,
    // candidates and the search's size bound.
    ResidualKernel,
    // `optimized` with residual-kernel's kernel extraction and nothing else:
    // the internal batch edges come from the endpoint bitmap and the batch
    // edge set (two-endpoint lemma, pair tests) instead of the anchor lists,
    // and the one-edge kernel is handed its endpoint positions instead of
    // searching for them.  No size gate, no star or complete closed forms, no
    // complement search: the candidates, their order and the probes are the
    // control's.  So `optimized` vs this isolates the extraction, and this vs
    // `residual-kernel` is what residual-kernel adds beyond it -- the fair
    // baseline.
    OptimizedBitmap,
    // Insertion controls with the same BK implementation, one-edge shortcut,
    // and index probes. Only the BK search domain differs.
    FullBK,
    KernelBK,
};

[[nodiscard]] constexpr bool isBKComparison(
    SubsumedCliqueAlgorithm algorithm) noexcept {
    return algorithm == SubsumedCliqueAlgorithm::FullBK ||
        algorithm == SubsumedCliqueAlgorithm::KernelBK;
}

[[nodiscard]] constexpr bool usesSubsumedSizeGate(
    SubsumedCliqueAlgorithm algorithm) noexcept {
    // The comparison arms start with the size histogram off by default.
    return algorithm == SubsumedCliqueAlgorithm::ResidualKernel;
}

// Whether the framework finds kernels through the endpoint bitmap (and so
// builds the anchor lists only on demand).
[[nodiscard]] constexpr bool usesBitmapKernelExtraction(
    SubsumedCliqueAlgorithm algorithm) noexcept {
    return algorithm == SubsumedCliqueAlgorithm::ResidualKernel ||
        algorithm == SubsumedCliqueAlgorithm::OptimizedBitmap ||
        isBKComparison(algorithm);
}

[[nodiscard]] const char* subsumedCliqueAlgorithmName(
    SubsumedCliqueAlgorithm algorithm) noexcept;

struct SubsumedCliqueSearchStats {
    long long batches = 0;
    long long newCliquesProcessed = 0;
    long long internalBatchEdges = 0;
    long long splitRounds = 0;
    long long candidateOccurrences = 0;
    long long duplicateCandidatesRemoved = 0;
    long long indexLookups = 0;
    long long emittedCliques = 0;
    long long batchEdgeVisits = 0;
    long long bitmaskCliques = 0;
    long long fallbackCliques = 0;
    // The size gate (residual-kernel only): cliques retired whole, and
    // candidates dropped before their probe.
    long long sizeGateSkippedCliques = 0;
    long long sizeGateSkippedProbes = 0;
    // Cliques whose kernel is a single batch edge, probed without a split.
    long long fastPathCliques = 0;
    long long deferredProbes = 0;
    long long inlineProbes = 0;
    long long structuralStarCliques = 0;
    long long structuralCompleteCliques = 0;
    // Candidates residual-kernel's star and complete closed forms and its
    // complement search emitted.
    long long kernelCandidates = 0;
    std::size_t maximumDeferredProbes = 0U;
    // The largest candidate set a split held: the controls only.
    std::size_t maximumCandidates = 0U;
    // Complement searches over kernels of cliques of at most 64 vertices, and
    // their nodes; the wider ones are counted in `wide`.
    long long residualKernelSearches = 0;
    long long residualKernelNodes = 0;
    wideresidual::Counters wide;
    // Deletion batches (`recoverBatch`): the candidates the destroyed cliques
    // gave, those below two vertices, those already accepted in this batch,
    // the maximality checks run on the rest, the candidates the checks
    // rejected and accepted, the neighbour entries and adjacency tests the
    // checks read, the centre groups checked on bit rows (`groups`), the
    // candidates checked that way (`groupedChecks`), and the neighbour
    // entries the rows read (`rowReads`).  Residual-kernel only: the
    // complement searches run under extension masks (`maskedSearches`), the
    // vertices extending their cores (`coreExtenders`), the masks kept
    // (`extensionMasks`), the search nodes a mask cut (`maskCuts`), and the
    // destroyed cliques one vertex extends whole (`extendedWhole`).
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
        long long coreExtenders = 0;
        long long extensionMasks = 0;
        long long maskCuts = 0;
        long long extendedWhole = 0;
    };
    Recovery recovery;
    // Per-batch phase split of the batched pipeline, read four times per
    // batch and only in a build with -DIMMCE_PHASES.
    double phasePrepareMilliseconds = 0.0;
    double phaseCliquesMilliseconds = 0.0;
    double phaseFlushMilliseconds = 0.0;
    double phaseFinishMilliseconds = 0.0;
};

// Finds the old maximal cliques a batch destroys.
//
// An old maximal clique M stops being maximal only inside a new maximal
// clique C that strictly contains it, and every vertex of C \ M is separated
// from M by an internal batch edge of C.  So the candidates for one new
// clique are obtained by splitting C on its internal batch edges -- each such
// edge keeps at most one of its endpoints -- and every candidate that the old
// maximal-clique index still holds is subsumed.
//
// Deletion is the same search run backwards (Das et al., Observation 2).  The
// maximal cliques a deletion batch destroys are those holding a deleted edge
// -- the new-clique search on the graph before the deletion finds them -- and
// every maximal clique the deletion creates is a candidate of one of them,
// produced by exactly the same split.  What differs is the test: those
// candidates are not in any index yet, so `recoverBatch` keeps the ones that
// are maximal in the graph after the deletion.
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

class SubsumedCliqueEnumerator {
public:
    // The endpoints of one internal batch edge.  Carrying them instead of an
    // edge id keeps the whole per-clique path off the global edge catalog,
    // which on the larger graphs is hundreds of megabytes of cold memory.
    struct InternalBatchEdge {
        Vertex u = 0;
        Vertex v = 0;
    };

    explicit SubsumedCliqueEnumerator(const IncrementalGraph& graph);

    // Batched pipeline: the whole batch of new maximal cliques at once.  The
    // subsumed cliques come back canonical and each exactly once, in the order
    // the search found them -- which differs between frameworks, so a caller
    // that compares or prints them sorts them itself.
    [[nodiscard]] std::vector<Clique> enumerateBatch(
        std::span<const EdgeId> batch,
        std::span<const Clique> newCliques,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        SubsumedCliqueAlgorithm algorithm);

    // Streaming pipeline: `prepareBatch` once per batch, then
    // `consumeNewClique` for each new maximal clique as the Lambda_new search
    // emits it, then `finishBatch`.  Subsumed cliques are reported through
    // `output`.
    void prepareBatch(
        std::span<const EdgeId> batch,
        SubsumedCliqueAlgorithm algorithm,
        SubsumedCliqueSearchStats& stats);
    void consumeNewClique(
        const Clique& canonicalNewClique,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueAlgorithm algorithm,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // Closes the batch: issues every probe the batch held back.  Holding them
    // to the end of the batch is sound because no new maximal clique of a
    // batch is a proper subset of another, so the insertions that run between
    // them cannot meet a candidate.
    void finishBatch(
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);

    // Deletion: the maximal cliques the deletion of `deletedBatch` creates.
    // The graph must already be without the batch.  `destroyed` are the
    // maximal cliques the batch destroyed (the new-clique search on the graph
    // before the deletion).  Each of their candidates with at least two
    // vertices that is maximal in the graph now is inserted into `cliques`
    // and returned, once however many destroyed cliques produce it: no
    // candidate can equal a clique the index held before the batch (it is a
    // proper subset of one), so the index itself tells a repeat.  The size
    // gate does not apply -- the targets are not in the index -- so the only
    // size bound any algorithm uses is the floor of two vertices.
    [[nodiscard]] std::vector<Clique> recoverBatch(
        std::span<const EdgeId> deletedBatch,
        std::span<const Clique> destroyed,
        MaximalCliqueIndex& cliques,
        SubsumedCliqueSearchStats& stats,
        SubsumedCliqueAlgorithm algorithm);

    // A candidate is a proper subset of the new clique that produced it, and
    // no other new maximal clique of the same batch can equal it -- both are
    // maximal in G+H, so neither contains the other.  Every probe of a batch
    // is therefore held back to the end of it, which turns two dependent
    // random accesses per clique into one buffer that can be prefetched
    // arbitrarily deep -- once the index is too large to stay cached, which no
    // unit-test index ever is.  Lowering the floor is what lets the tests
    // reach the buffered path at all; production leaves it at the measured
    // default.
    void setDeferralMinimumIndexSize(std::size_t minimumSize) noexcept {
        deferralMinimumIndexSize_ = minimumSize;
    }
    // Set hashes of the new cliques this batch has consumed, in order, so the
    // caller can insert them without hashing a second time.
    [[nodiscard]] std::span<const std::uint64_t> newCliqueHashes()
        const noexcept {
        return newCliqueHashes_;
    }

private:
    void validateBatchEdges() const;
    void consumeNewCliqueResolved(
        const Clique& clique,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueAlgorithm algorithm,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);

    // optimized: the anchor-index extraction, then the split.
    void consumeNewCliqueOptimized(
        const Clique& clique,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // Anchors every batch edge at its lower-degree endpoint, so a clique only
    // has to walk the incidence lists of its own vertices to find the batch
    // edges inside it.  Built once per batch.
    //
    // The opposite endpoint is stored next to each incidence entry, so the
    // per-clique loop never reads the global edge catalog -- on these graphs
    // that array is megabytes, well past L2, while the |V|-sized epoch arrays
    // below are a few hundred kilobytes and stay cached.  (A compact
    // open-addressed table over the batch's own vertices was measured here
    // and is a net loss at these vertex counts: it replaces a cached array
    // index with a hash and a probe.)
    void prepareBatchIndex();
    void ensureBatchAnchors() {
        if (!batchAnchorsPrepared_) {
            prepareBatchIndex();
            batchAnchorsPrepared_ = true;
        }
    }
    // Returns a view of `internalEdgeBuffer_`, valid until the next call on
    // this enumerator.
    [[nodiscard]] std::span<const InternalBatchEdge> internalBatchEdges(
        std::span<const Vertex> clique,
        SubsumedCliqueSearchStats* stats);
    [[nodiscard]] std::uint32_t beginBatchIndexEpoch();

    // optimized-bitmap: the bitmap extraction, then the control's split.
    void consumeNewCliqueOptimizedBitmap(
        const Clique& clique,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // The kernel of S = `kernelPositions_` in the control's batch order, into
    // `kernelEdgeBuffer_` (with positions in `kernelEdgePositions_`).
    void collectKernelBatchOrder(
        const Clique& clique,
        SubsumedCliqueSearchStats& stats);

    // residual-kernel: the size gate and the bitmap extraction, then the
    // kernel's closed form or complement search.
    void consumeNewCliqueResidualKernel(
        const Clique& clique,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    void consumeResidualKernelPositions(
        const Clique& clique,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    void consumeResidualKernelEdges(
        const Clique& clique,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // The kernel of a clique of at most 64 vertices, given as conflict rows
    // over clique positions: its star or complete closed form, or the
    // complement Bron-Kerbosch.
    void consumeKernelRows64(
        const Clique& clique,
        const std::array<std::uint64_t, 64U>& conflicts,
        std::uint64_t touched,
        std::size_t edgeCount,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);

    // The bitmap extraction's batch structures, shared by optimized-bitmap
    // and residual-kernel, built on first use so a batch whose every clique
    // the size gate retires pays for none of them.
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
    // The edge's position in the batch, or -1 (optimized-bitmap only: the
    // control splits in batch order).
    [[nodiscard]] std::int64_t batchEdgeSetIndex(
        Vertex u, Vertex v) const noexcept;

    // The one-edge kernel, every framework's: C minus either endpoint, the
    // endpoints at clique positions `uIndex < vIndex`.  `sizeGate` drops both
    // candidates when the index holds no clique of their size.
    void probeOneEdgeKernel(
        const Clique& clique,
        std::size_t uIndex,
        std::size_t vIndex,
        std::uint64_t setAccumulator,
        bool sizeGate,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // The controls' split of C on its internal batch edges, in their order.
    void consumeNewCliqueSplit(
        const Clique& clique,
        std::span<const InternalBatchEdge> cliqueBatchEdges,
        std::uint64_t setAccumulator,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // Issues one candidate probe, or buffers it when this batch is deferring.
    // `arenaOffset` is where the producing clique's vertices were copied into
    // `probeArena_`, or `noArenaOffset` if they have not been copied yet.
    void submitProbe(
        std::span<const Vertex> clique,
        std::span<const std::uint64_t> maskWords,
        std::size_t memberCount,
        std::uint64_t setHash,
        bool defer,
        std::uint32_t& arenaOffset,
        MaximalCliqueIndex& oldCliques,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    [[nodiscard]] bool deferringProbes(
        const MaximalCliqueIndex& oldCliques) const noexcept;
    static constexpr std::uint32_t noArenaOffset = 0xFFFFFFFFU;

    // The candidate-size floor and test.  Inserting, they are the size gate:
    // the old cliques' live sizes (with the histogram off, no bound at all).
    // Recovering, the targets are not in the index, so only the fixed floor
    // of two vertices holds.
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
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    void acceptRecoveredClique(
        std::span<const Vertex> candidate,
        std::uint64_t setHash,
        SubsumedCliqueSearchStats& stats,
        CliqueBatchOutput& output);
    // Whether a clique of the graph as it is now has no vertex extending it.
    [[nodiscard]] bool maximalInGraph(
        std::span<const Vertex> candidate,
        SubsumedCliqueSearchStats& stats);
    // The candidate's adjacency lists into `memberLists_`, returning the
    // least degree; and the check on them.
    std::size_t readMemberLists(std::span<const Vertex> candidate);
    // `skipped` holds the candidate and is skipped in the scan: the
    // candidate itself, or the whole destroyed clique it came from.
    [[nodiscard]] bool maximalOnMemberLists(
        std::span<const Vertex> candidate,
        std::span<const Vertex> skipped,
        SubsumedCliqueSearchStats& stats);
    // Recovering: whether the batch already accepted an equal candidate,
    // and the record of one it accepts; every candidate the batch accepts
    // is in this table until the end of the batch, when they all enter the
    // index.
    [[nodiscard]] bool acceptedInBatch(
        std::span<const Vertex> candidate,
        std::uint64_t setHash) const noexcept;
    void recordAccepted(std::span<const Vertex> candidate, std::uint64_t setHash);
    void insertAcceptedInBatch(MaximalCliqueIndex& cliques);
    // Residual-kernel recovery: each vertex outside a destroyed clique that
    // extends its core (the clique positions `core`), as the kernel endpoints
    // it misses (bits of `endpoints`, T's numbering, at most 64 of them),
    // handed to `onExtender`, which returns true to end the scan.
    template <typename OnExtender>
    void scanCoreExtenders(
        const Clique& clique,
        std::span<const std::uint32_t> core,
        std::span<const std::uint32_t> endpoints,
        SubsumedCliqueSearchStats& stats,
        OnExtender&& onExtender);
    // The complement search's masks: into `extensionMasks_`, the minimal
    // ones only.  False when some vertex misses none: it extends every
    // candidate.
    [[nodiscard]] bool collectExtensionMasks(
        const Clique& clique,
        std::span<const std::uint32_t> core,
        std::span<const std::uint32_t> endpoints,
        SubsumedCliqueSearchStats& stats);

    // A recovery checks a candidate on the spot unless every member has a
    // high degree; those wait for the end of the batch, where each distinct
    // one is checked once.  Each then gets the member the most of them hold
    // -- as a rule a vertex that many destroyed cliques share -- as its
    // centre, and every vertex extending a candidate is a neighbour of its
    // centre.  A large group of one centre is checked on bit rows over the
    // centre's neighbourhood, each member's row read once for the whole
    // group, when that reads less than checking its candidates one by one;
    // the rest are checked one by one in the order they were made.  It only
    // needs the candidates to be subsets of destroyed cliques, so every
    // algorithm has it since 2026-10-14 (residual-kernel alone before).  A
    // batch that gathers more than `gatherLimit` vertices has them checked
    // then and there, and gathers on: coauthor's window of batch 200000
    // gathers 5.2 billion in one batch (2026-10-12).
    struct GatheredCandidate {
        std::uint64_t setHash = 0U;
        std::size_t offset = 0U;
        std::uint32_t size = 0U;
        Vertex centre = 0;
        bool onRows = false;
    };
    void checkGatheredCandidates(
        SubsumedCliqueSearchStats& stats,
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
        SubsumedCliqueSearchStats& stats);

    const IncrementalGraph& graph_;
    // Set for the length of `recoverBatch`.
    bool recovering_ = false;
    std::vector<Vertex> recoveredScratch_;
    // A maximality check's members by increasing degree: (degree, position
    // in the candidate), with their adjacency lists and the galloping
    // cursors into them.
    std::vector<std::pair<std::size_t, std::uint32_t>> memberOrder_;
    std::vector<std::span<const AdjacencyEntry>> memberLists_;
    std::vector<std::size_t> memberCursors_;
    // The destroyed clique a recovery is consuming (residual-kernel's checks
    // on the spot skip it).
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
    // Set for the length of a residual-kernel `recoverBatch`.  Every
    // algorithm's recovery checks its candidates the same way: the deferred
    // ones grouped on bit rows or one by one at the end of the batch, the
    // rest on the spot by a galloping scan that stops at the first extender,
    // the duplicates found in the batch's own table, and the index filled
    // at the end of the batch.  Residual-kernel's alone: its scans on the
    // spot skip the whole destroyed clique (its candidates are the clique's
    // maximal subsets free of deleted edges; the other algorithms' splits
    // also leave subsets of other candidates, which a left-out vertex does
    // extend), and its complement searches run under extension masks.
    // (Until 2026-10-14 the other algorithms found every vertex extending a
    // candidate before deciding, on the spot, with the index probed and
    // filled one candidate at a time: the user's decision of 2026-10-11,
    // which the user then reversed for all but these two.)
    bool kernelRecovery_ = false;
    // Set while a complement search runs under extension masks: what it
    // emits is maximal, so the candidate is accepted unless a duplicate.
    bool knownMaximal_ = false;
    // Whether the batch accepted a candidate that a deferred one can equal
    // (under extension masks, or in an earlier round of the deferred ones),
    // so that those are looked up in the batch's table.
    bool filterGathered_ = false;
    // `collectExtensionMasks`'s lists: the core's by increasing length, and
    // the kernel endpoints', each with a galloping cursor.
    std::vector<std::pair<std::size_t, std::uint32_t>> coreOrder_;
    std::vector<std::span<const AdjacencyEntry>> coreLists_;
    std::vector<std::size_t> coreCursors_;
    std::vector<std::span<const AdjacencyEntry>> endpointLists_;
    std::vector<std::size_t> endpointCursors_;
    std::vector<std::uint64_t> extensionMasks_;
    // A destroyed clique's core and kernel endpoints as clique positions,
    // and a wide clique's kernel degrees and candidate mask.
    std::vector<std::uint32_t> corePositions_;
    std::vector<std::uint32_t> endpointPositions_;
    std::vector<std::uint32_t> wideDegree_;
    std::vector<std::uint64_t> wideMask_;
    // A recovering complement search's candidates, held (kept endpoints and
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
    SubsumedCliqueAlgorithm preparedAlgorithm_ =
        SubsumedCliqueAlgorithm::Optimized;

    // The head of a vertex's batch-edge incidence list, with the epoch that
    // says whether the head belongs to this batch.  One array rather than two
    // parallel ones: the pair is always read together, and at |V| entries each
    // the separate arrays put two cache misses where there can be one.
    struct BatchAnchor {
        std::uint32_t epoch = 0U;
        EdgeId head = invalidEdgeId;
    };
    std::vector<BatchAnchor> batchAnchorByVertex_;
    std::uint32_t batchIndexEpoch_ = 0U;
    std::vector<Vertex> batchOtherByIncidence_;
    std::vector<EdgeId> nextBatchIncidence_;
    std::vector<InternalBatchEdge> internalEdgeBuffer_;
    std::vector<std::pair<EdgeId, InternalBatchEdge>> rankedEdgeBuffer_;
    bool batchAnchorsPrepared_ = false;
    // The bitmap extraction's batch structures.  The bitmap has one bit per
    // vertex (a few hundred kilobytes on the largest graph, so it stays in
    // L2), and the edge set holds the batch edges as (min << 32 | max) keys in
    // an open-addressed table at most a quarter full.
    std::vector<std::uint64_t> endpointBits_;
    std::vector<Vertex> endpointList_;
    std::vector<std::uint64_t> batchEdgeSet_;
    // Parallel to `batchEdgeSet_`: each stored edge's position in the batch,
    // kept for optimized-bitmap only.
    std::vector<std::uint32_t> batchEdgeSetOrder_;
    std::vector<std::uint32_t> batchEdgeSlots_;
    unsigned batchEdgeSetShift_ = 64U;
    bool bitmapExtractionPrepared_ = false;
    // A clique's kernel: S as clique positions, the kernel edges, and their
    // endpoints' clique positions.
    std::vector<std::uint32_t> kernelPositions_;
    std::vector<InternalBatchEdge> kernelEdgeBuffer_;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> kernelEdgePositions_;
    // optimized-bitmap scratch: kernel edges with their batch positions.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> kernelBatchOrder_;
    // Below this many live cliques a probe has no latency worth hiding and
    // the buffering is pure overhead: 4,096 live cliques is a 64 KB table,
    // the point where it stops fitting in L1, and that is where hiding the
    // record and arena accesses of the staged flush starts to pay.  Measured
    // over 1,024 to 16,384 the full-maintenance time is flat to within the
    // noise on both datasets (dblp-cite 7.76-8.10 ms, facebook 63.2-64.4 ms)
    // while 65,536 costs 2.7% and 3.9% respectively, so the value is chosen
    // from the cache argument rather than fitted.
    std::size_t deferralMinimumIndexSize_ = std::size_t{1} << 12U;
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
        SubsumedCliqueSearchStats& stats,
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

} // namespace immce
