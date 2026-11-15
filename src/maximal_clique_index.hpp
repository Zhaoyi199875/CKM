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

class MaximalCliqueIndex {
public:
    explicit MaximalCliqueIndex(bool trackSizes = false);

    [[nodiscard]] bool contains(std::span<const Vertex> clique) const;
    [[nodiscard]] bool containsCanonical(const Clique& clique) const;
    [[nodiscard]] bool containsHashed(
        std::span<const Vertex> clique,
        std::uint64_t setHash) const noexcept;
    bool insert(Clique clique);
    bool insert(std::span<const Vertex> clique);
    bool insertHashed(std::span<const Vertex> clique, std::uint64_t setHash);
    bool eraseCanonical(const Clique& clique);
    [[nodiscard]] std::optional<Clique> extractMaskedHashed(
        std::span<const Vertex> vertices,
        std::span<const std::uint64_t> membershipWords,
        std::size_t memberCount,
        std::uint64_t setHash);

    void prefetch(std::uint64_t setHash) const noexcept;
    void prefetchForInsert(std::uint64_t setHash) const noexcept;

    static constexpr std::size_t noSlot = static_cast<std::size_t>(-1);

    [[nodiscard]] std::size_t firstTagMatch(
        std::uint64_t setHash) const noexcept;
    [[nodiscard]] std::uint32_t recordAtSlot(std::size_t slot) const noexcept;
    void prefetchRecord(std::uint32_t record) const noexcept;
    void prefetchRecordVertices(std::uint32_t record) const noexcept;
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

}
