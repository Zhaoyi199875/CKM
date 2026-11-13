#pragma once

#include "types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace ckm {

using CliqueId = std::size_t;

// A commutative set hash: the accumulator is the XOR of one scrambled word
// per vertex, so dropping a vertex from a clique costs one XOR instead of a
// full rehash.  The final mix folds in the member count.
[[nodiscard]] inline std::uint64_t vertexSetWord(Vertex vertex) noexcept {
    std::uint64_t value =
        static_cast<std::uint64_t>(static_cast<std::uint32_t>(vertex)) +
        0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] inline std::uint64_t finishSetHash(
    std::uint64_t accumulator,
    std::size_t memberCount) noexcept {
    std::uint64_t value = accumulator +
        0x9E3779B97F4A7C15ULL * (static_cast<std::uint64_t>(memberCount) + 1U);
    value ^= value >> 33U;
    value *= 0xFF51AFD7ED558CCDULL;
    value ^= value >> 33U;
    value *= 0xC4CEB9FE1A85EC53ULL;
    return value ^ (value >> 33U);
}

[[nodiscard]] inline std::uint64_t setHashOf(std::span<const Vertex> clique) {
    std::uint64_t accumulator = 0U;
    for (Vertex vertex : clique) {
        accumulator ^= vertexSetWord(vertex);
    }
    return finishSetHash(accumulator, clique.size());
}

// Flat, open-addressed store of the current maximal cliques.
//
// Each table slot carries the full 64-bit set hash next to the record id, so
// a probe that misses is decided by one cache line and never follows a
// pointer.  Clique vertices live packed in one arena instead of a vector per
// clique, which removes two allocations per insertion.
//
// Spans handed out by `cliqueVertices` stay valid until the next insertion;
// insertions may compact the arena.
class MaximalCliqueIndex {
public:
    explicit MaximalCliqueIndex(bool trackSizes = false);

    [[nodiscard]] bool contains(std::span<const Vertex> clique) const;
    [[nodiscard]] bool containsCanonical(const Clique& clique) const;
    // For a caller that already holds the canonical clique's set hash.
    [[nodiscard]] bool containsHashed(
        std::span<const Vertex> clique,
        std::uint64_t setHash) const noexcept;
    bool insert(Clique clique);
    bool insert(std::span<const Vertex> clique);
    // For a caller that already holds the canonical clique's set hash.  The
    // subsumed search computes it anyway to derive its candidates, so handing
    // it over saves a second pass over the vertices here.
    bool insertHashed(std::span<const Vertex> clique, std::uint64_t setHash);
    bool eraseCanonical(const Clique& clique);
    // Looks up the subset of `vertices` that `membershipWords` selects, and
    // erases it if it is there.  A subsumed candidate is its new clique under
    // a mask, so this is what lets the search probe for one without building
    // the set first.  The caller supplies the masked subset's set hash, which
    // it derives from the producing clique's own by one XOR per dropped
    // vertex rather than by hashing the subset from scratch.
    [[nodiscard]] std::optional<Clique> extractMaskedHashed(
        std::span<const Vertex> vertices,
        std::span<const std::uint64_t> membershipWords,
        std::size_t memberCount,
        std::uint64_t setHash);

    // Start the memory fetch for a slot that is about to be read or written.
    // The subsumed search knows both a clique ahead of using them.
    void prefetch(std::uint64_t setHash) const noexcept;
    void prefetchForInsert(std::uint64_t setHash) const noexcept;

    // A probe that misses is decided by the table alone -- the tag rejects it
    // before anything else is read -- but a probe that hits then follows two
    // more pointers, into the record and into the arena, and each of those
    // depends on the one before it.  Prefetching cannot hide a dependent
    // chain, so a caller with many probes in hand splits the chain into
    // stages instead: resolve every probe against the table, then fetch every
    // surviving record, then verify against the arena.  These three entry
    // points are what let it do that.
    static constexpr std::size_t noSlot = static_cast<std::size_t>(-1);

    // Stage one: the first slot of this hash's run whose tag matches, decided
    // without reading the record or the arena.  A tag collision is possible,
    // so the slot is a candidate, not yet a hit.
    [[nodiscard]] std::size_t firstTagMatch(
        std::uint64_t setHash) const noexcept;
    [[nodiscard]] std::uint32_t recordAtSlot(std::size_t slot) const noexcept;
    void prefetchRecord(std::uint32_t record) const noexcept;
    void prefetchRecordVertices(std::uint32_t record) const noexcept;
    // Stage three: verify the candidate at `slot` and erase it if it is the
    // one being looked for, resuming the run from there if it is not.  Safe
    // against the slot having been erased since stage one -- a tombstone does
    // not end a run -- which is what a duplicate candidate produces.
    [[nodiscard]] std::optional<Clique> extractMaskedAtSlot(
        std::size_t slot,
        std::span<const Vertex> vertices,
        std::span<const std::uint64_t> membershipWords,
        std::size_t memberCount,
        std::uint64_t setHash);

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;

    [[nodiscard]] bool tracksSizes() const noexcept { return trackSizes_; }

    static constexpr std::size_t trackedSizeLimit = 64U;

    [[nodiscard]] std::size_t minimumLiveCliqueSize() const noexcept {
        return trackSizes_ ? minimumLiveSize_ : 0U;
    }
    [[nodiscard]] bool hasLiveCliqueOfSize(std::size_t size) const noexcept {
        if (!trackSizes_) {
            return true;
        }
        if (size >= trackedSizeLimit) {
            return oversizeLiveCount_ != 0U;
        }
        return liveCountBySize_[size] != 0U;
    }
    [[nodiscard]] std::size_t liveCliqueCountOfSize(
        std::size_t size) const noexcept {
        return size < trackedSizeLimit ? liveCountBySize_[size] : 0U;
    }
    [[nodiscard]] std::size_t oversizeLiveCliqueCount() const noexcept {
        return oversizeLiveCount_;
    }

private:
    // Eight bytes per slot: eight slots share a cache line, so a linear
    // probe run costs one miss, and the whole table is half the size it
    // would be with the full hash inline.  The tag only has to make a
    // vertex-by-vertex verification rare, not to decide membership.
    struct Slot {
        std::uint32_t tag = 0U;
        std::uint32_t record = 0xFFFFFFFFU;
    };

    [[nodiscard]] static std::uint32_t slotTag(
        std::uint64_t setHash) noexcept {
        return static_cast<std::uint32_t>(setHash >> 32U);
    }

    static constexpr std::uint32_t emptyRecord = 0xFFFFFFFFU;
    static constexpr std::uint32_t tombstoneRecord = 0xFFFFFFFEU;
    static constexpr std::size_t notFound =
        static_cast<std::size_t>(-1);

    [[nodiscard]] std::size_t findExact(
        std::uint64_t setHash,
        std::span<const Vertex> clique) const noexcept;
    // One probe run that both rejects a duplicate and picks the slot to
    // write, instead of walking the run twice.
    [[nodiscard]] std::size_t findOrClaim(
        std::uint64_t setHash,
        std::span<const Vertex> clique,
        bool& duplicate) noexcept;
    [[nodiscard]] std::size_t findMasked(
        std::uint64_t setHash,
        std::span<const Vertex> vertices,
        std::span<const std::uint64_t> membershipWords,
        std::size_t memberCount) const noexcept;
    [[nodiscard]] bool recordEquals(
        std::uint32_t record,
        std::span<const Vertex> clique) const noexcept;
    [[nodiscard]] bool recordEqualsMasked(
        std::uint32_t record,
        std::span<const Vertex> vertices,
        std::span<const std::uint64_t> membershipWords,
        std::size_t memberCount) const noexcept;
    void rehash(std::size_t slotCount);
    // Drops the records of erased cliques and renumbers the survivors, and
    // reports the renumbering so the caller can carry the table across it.
    // Only legal while nothing outside holds a clique id, which nothing does:
    // the index hands its cliques out by value.
    void compactRecords(std::vector<std::uint32_t>& remap);
    void reserveSlot();
    void compactArena();
    [[nodiscard]] bool isAlive(CliqueId cliqueId) const noexcept;
    [[nodiscard]] std::span<const Vertex> cliqueVertices(
        CliqueId cliqueId) const;
    [[nodiscard]] Clique cliqueAt(std::uint32_t record) const;
    void eraseRecord(std::size_t slot, std::uint32_t record);

    std::vector<Slot> table_;
    std::size_t tableMask_ = 0U;
    std::size_t occupiedSlots_ = 0U;
    std::size_t liveCount_ = 0U;

    // One record per live clique, twelve bytes.
    //
    // The third word used to be the slot the clique sits in, and the set hash
    // was recomputed from the arena whenever a rehash needed it.  That made a
    // rehash read the vertices of every live clique -- measured at 27% of all
    // index-update time on both datasets, because it is the one pass that
    // touches the whole arena.  Keeping the low half of the hash instead
    // removes those reads: a rehash walks the old table, which already holds
    // the high half as its tag, and pairs it with this.
    //
    // The slot is not missed.  Both callers of `eraseRecord` have just probed
    // for the record and so know where it lives, and a rehash writes each slot
    // where it puts it.
    struct Record {
        std::uint32_t offset = 0U;
        std::uint32_t size = 0U;
        std::uint32_t hashLow = 0U;
    };
    std::vector<Record> records_;
    std::size_t deadRecords_ = 0U;
    std::vector<Vertex> arena_;
    std::size_t deadArenaVertices_ = 0U;

    void noteInserted(std::size_t size);
    void noteErased(std::size_t size);
    void recomputeMinimumLiveSize(std::size_t from) noexcept;

    bool trackSizes_ = false;
    std::array<std::uint32_t, trackedSizeLimit> liveCountBySize_{};
    std::size_t oversizeLiveCount_ = 0U;
    std::size_t minimumLiveSize_ = static_cast<std::size_t>(-1);
};

} // namespace ckm
