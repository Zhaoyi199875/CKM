#include "maximal_clique_index.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace ckm {
namespace {

constexpr std::size_t minimumTableSlots = 1U << 10U;
constexpr std::size_t minimumArenaCompaction = 1U << 16U;

}

MaximalCliqueIndex::MaximalCliqueIndex(bool trackSizes)
    : trackSizes_(trackSizes) {
    rehash(minimumTableSlots);
}

bool MaximalCliqueIndex::recordEquals(
    std::uint32_t record,
    std::span<const Vertex> clique) const noexcept {
    if (records_[record].size != clique.size()) {
        return false;
    }
    const Vertex* stored = arena_.data() + records_[record].offset;
    for (std::size_t index = 0; index < clique.size(); ++index) {
        if (stored[index] != clique[index]) {
            return false;
        }
    }
    return true;
}

bool MaximalCliqueIndex::recordEqualsMasked(
    std::uint32_t record,
    std::span<const Vertex> vertices,
    std::span<const std::uint64_t> membershipWords,
    std::size_t memberCount) const noexcept {
    if (records_[record].size != memberCount) {
        return false;
    }
    const Vertex* stored = arena_.data() + records_[record].offset;
    std::size_t storedIndex = 0U;
    for (std::size_t index = 0; index < vertices.size(); ++index) {
        const std::uint64_t bit = std::uint64_t{1} << (index % 64U);
        if ((membershipWords[index / 64U] & bit) == 0U) {
            continue;
        }
        if (stored[storedIndex] != vertices[index]) {
            return false;
        }
        ++storedIndex;
    }
    return true;
}

std::size_t MaximalCliqueIndex::findExact(
    std::uint64_t setHash,
    std::span<const Vertex> clique) const noexcept {
    const std::uint32_t tag = slotTag(setHash);
    std::size_t slot = static_cast<std::size_t>(setHash) & tableMask_;
    while (true) {
        const Slot& entry = table_[slot];
        if (entry.record == emptyRecord) {
            return notFound;
        }
        if (entry.record != tombstoneRecord && entry.tag == tag &&
            recordEquals(entry.record, clique)) {
            return slot;
        }
        slot = (slot + 1U) & tableMask_;
    }
}

std::size_t MaximalCliqueIndex::findMasked(
    std::uint64_t setHash,
    std::span<const Vertex> vertices,
    std::span<const std::uint64_t> membershipWords,
    std::size_t memberCount) const noexcept {
    const std::uint32_t tag = slotTag(setHash);
    std::size_t slot = static_cast<std::size_t>(setHash) & tableMask_;
    while (true) {
        const Slot& entry = table_[slot];
        if (entry.record == emptyRecord) {
            return notFound;
        }
        if (entry.record != tombstoneRecord && entry.tag == tag &&
            recordEqualsMasked(
                entry.record, vertices, membershipWords, memberCount)) {
            return slot;
        }
        slot = (slot + 1U) & tableMask_;
    }
}

void MaximalCliqueIndex::rehash(std::size_t slotCount) {
    if (slotCount > (std::size_t{1} << 32U)) {
        throw std::length_error("maximal-clique index is too large");
    }
    std::vector<Slot> previous;
    previous.swap(table_);
    std::vector<std::uint32_t> remap;
    if (deadRecords_ * 2U > records_.size()) {
        compactRecords(remap);
    }
    table_.assign(slotCount, Slot{});
    tableMask_ = slotCount - 1U;
    occupiedSlots_ = 0U;
    for (const Slot& entry : previous) {
        if (entry.record == emptyRecord || entry.record == tombstoneRecord) {
            continue;
        }
        const std::uint32_t record =
            remap.empty() ? entry.record : remap[entry.record];
        if (record == emptyRecord) {
            continue;
        }
        std::size_t slot =
            static_cast<std::size_t>(records_[record].hashLow) & tableMask_;
        while (table_[slot].record != emptyRecord) {
            slot = (slot + 1U) & tableMask_;
        }
        table_[slot] = Slot{entry.tag, record};
        ++occupiedSlots_;
    }
}

std::size_t MaximalCliqueIndex::findOrClaim(
    std::uint64_t setHash,
    std::span<const Vertex> clique,
    bool& duplicate) noexcept {
    const std::uint32_t tag = slotTag(setHash);
    std::size_t slot = static_cast<std::size_t>(setHash) & tableMask_;
    std::size_t target = notFound;
    while (true) {
        const Slot& entry = table_[slot];
        if (entry.record == emptyRecord) {
            if (target == notFound) {
                target = slot;
                ++occupiedSlots_;
            }
            duplicate = false;
            return target;
        }
        if (entry.record == tombstoneRecord) {
            if (target == notFound) {
                target = slot;
            }
        } else if (entry.tag == tag && recordEquals(entry.record, clique)) {
            duplicate = true;
            return slot;
        }
        slot = (slot + 1U) & tableMask_;
    }
}

void MaximalCliqueIndex::reserveSlot() {
    if ((occupiedSlots_ + 1U) * 4U < table_.size() * 3U) {
        return;
    }
    std::size_t slotCount = minimumTableSlots;
    while (slotCount < (liveCount_ + 1U) * 2U) {
        slotCount *= 2U;
    }
    rehash(slotCount);
}

void MaximalCliqueIndex::compactArena() {
    std::vector<Vertex> compacted;
    compacted.reserve(arena_.size() - deadArenaVertices_);
    for (Record& entry : records_) {
        if (entry.size == 0U) {
            continue;
        }
        const Vertex* stored = arena_.data() + entry.offset;
        entry.offset = static_cast<std::uint32_t>(compacted.size());
        compacted.insert(compacted.end(), stored, stored + entry.size);
    }
    arena_.swap(compacted);
    deadArenaVertices_ = 0U;
}

bool MaximalCliqueIndex::insert(Clique clique) {
    clique = canonicalClique(std::move(clique));
    return insert(std::span<const Vertex>(clique));
}

bool MaximalCliqueIndex::insert(std::span<const Vertex> clique) {
    return insertHashed(clique, setHashOf(clique));
}

bool MaximalCliqueIndex::insertHashed(
    std::span<const Vertex> clique,
    std::uint64_t setHash) {
    if (arena_.size() >= minimumArenaCompaction &&
        deadArenaVertices_ * 2U > arena_.size()) {
        compactArena();
    }
    if (records_.size() >= tombstoneRecord ||
        arena_.size() + clique.size() > std::numeric_limits<std::uint32_t>::max())
        [[unlikely]] {
        throw std::length_error(
            "the maximal-clique index outgrew its 32-bit offsets");
    }
    reserveSlot();

    bool duplicate = false;
    const std::size_t target = findOrClaim(setHash, clique, duplicate);
    if (duplicate) {
        return false;
    }

    const auto record = static_cast<std::uint32_t>(records_.size());
    records_.push_back(Record{
        static_cast<std::uint32_t>(arena_.size()),
        static_cast<std::uint32_t>(clique.size()),
        static_cast<std::uint32_t>(setHash)});
    arena_.insert(arena_.end(), clique.begin(), clique.end());
    table_[target] = Slot{slotTag(setHash), record};
    ++liveCount_;
    noteInserted(clique.size());
    return true;
}

void MaximalCliqueIndex::noteInserted(std::size_t size) {
    if (!trackSizes_) {
        return;
    }
    if (size >= trackedSizeLimit) [[unlikely]] {
        ++oversizeLiveCount_;
        if (trackedSizeLimit < minimumLiveSize_) {
            minimumLiveSize_ = trackedSizeLimit;
        }
        return;
    }
    ++liveCountBySize_[size];
    if (size < minimumLiveSize_) {
        minimumLiveSize_ = size;
    }
}

void MaximalCliqueIndex::noteErased(std::size_t size) {
    if (!trackSizes_) {
        return;
    }
    if (size >= trackedSizeLimit) [[unlikely]] {
        if (oversizeLiveCount_ == 0U) {
            throw std::logic_error("live clique size histogram underflow");
        }
        --oversizeLiveCount_;
        if (oversizeLiveCount_ == 0U && minimumLiveSize_ == trackedSizeLimit) {
            recomputeMinimumLiveSize(0U);
        }
        return;
    }
    std::uint32_t& count = liveCountBySize_[size];
    if (count == 0U) {
        throw std::logic_error("live clique size histogram underflow");
    }
    --count;
    if (count != 0U || size != minimumLiveSize_) {
        return;
    }
    recomputeMinimumLiveSize(size + 1U);
}

void MaximalCliqueIndex::recomputeMinimumLiveSize(std::size_t from) noexcept {
    for (std::size_t candidate = from; candidate < trackedSizeLimit;
         ++candidate) {
        if (liveCountBySize_[candidate] != 0U) {
            minimumLiveSize_ = candidate;
            return;
        }
    }
    minimumLiveSize_ = oversizeLiveCount_ != 0U
        ? trackedSizeLimit
        : static_cast<std::size_t>(-1);
}

Clique MaximalCliqueIndex::cliqueAt(std::uint32_t record) const {
    const std::span<const Vertex> vertices = cliqueVertices(record);
    return Clique(vertices.begin(), vertices.end());
}

void MaximalCliqueIndex::eraseRecord(std::size_t slot, std::uint32_t record) {
    Record& entry = records_[record];
    noteErased(entry.size);
    table_[slot].record = tombstoneRecord;
    deadArenaVertices_ += entry.size;
    entry.size = 0U;
    ++deadRecords_;
    --liveCount_;
}

bool MaximalCliqueIndex::contains(std::span<const Vertex> clique) const {
    const Clique key = canonicalClique(Clique(clique.begin(), clique.end()));
    return containsCanonical(key);
}

bool MaximalCliqueIndex::containsCanonical(const Clique& clique) const {
    return findExact(setHashOf(clique), clique) != notFound;
}

bool MaximalCliqueIndex::containsHashed(
    std::span<const Vertex> clique,
    std::uint64_t setHash) const noexcept {
    return findExact(setHash, clique) != notFound;
}

bool MaximalCliqueIndex::eraseCanonical(const Clique& clique) {
    const std::size_t slot = findExact(setHashOf(clique), clique);
    if (slot == notFound) {
        return false;
    }
    eraseRecord(slot, table_[slot].record);
    return true;
}

std::optional<Clique> MaximalCliqueIndex::extractMaskedHashed(
    std::span<const Vertex> vertices,
    std::span<const std::uint64_t> membershipWords,
    std::size_t memberCount,
    std::uint64_t setHash) {
    const std::size_t slot =
        findMasked(setHash, vertices, membershipWords, memberCount);
    if (slot == notFound) {
        return std::nullopt;
    }
    const std::uint32_t record = table_[slot].record;
    Clique result = cliqueAt(record);
    eraseRecord(slot, record);
    return result;
}

bool MaximalCliqueIndex::isAlive(CliqueId cliqueId) const noexcept {
    return cliqueId < records_.size() && records_[cliqueId].size != 0U;
}

void MaximalCliqueIndex::compactRecords(std::vector<std::uint32_t>& remap) {
    remap.assign(records_.size(), emptyRecord);
    std::uint32_t next = 0U;
    for (std::uint32_t record = 0;
         record < static_cast<std::uint32_t>(records_.size()); ++record) {
        if (records_[record].size == 0U) {
            continue;
        }
        remap[record] = next;
        records_[next] = records_[record];
        ++next;
    }
    records_.resize(next);
    records_.shrink_to_fit();
    deadRecords_ = 0U;
}

std::span<const Vertex> MaximalCliqueIndex::cliqueVertices(
    CliqueId cliqueId) const {
    if (!isAlive(cliqueId)) {
        throw std::out_of_range("clique ID is not alive");
    }
    return std::span<const Vertex>(
        arena_.data() + records_[cliqueId].offset, records_[cliqueId].size);
}

std::size_t MaximalCliqueIndex::size() const noexcept {
    return liveCount_;
}

bool MaximalCliqueIndex::empty() const noexcept {
    return liveCount_ == 0U;
}

std::size_t MaximalCliqueIndex::firstTagMatch(
    std::uint64_t setHash) const noexcept {
    const std::uint32_t tag = slotTag(setHash);
    std::size_t slot = static_cast<std::size_t>(setHash) & tableMask_;
    while (true) {
        const Slot& entry = table_[slot];
        if (entry.record == emptyRecord) {
            return noSlot;
        }
        if (entry.record != tombstoneRecord && entry.tag == tag) {
            return slot;
        }
        slot = (slot + 1U) & tableMask_;
    }
}

std::uint32_t MaximalCliqueIndex::recordAtSlot(
    std::size_t slot) const noexcept {
    return table_[slot].record;
}

void MaximalCliqueIndex::prefetchRecord(
    std::uint32_t record) const noexcept {
    __builtin_prefetch(records_.data() + record, 0, 1);
}

void MaximalCliqueIndex::prefetchRecordVertices(
    std::uint32_t record) const noexcept {
    __builtin_prefetch(arena_.data() + records_[record].offset, 0, 1);
}

std::optional<Clique> MaximalCliqueIndex::extractMaskedAtSlot(
    std::size_t slot,
    std::span<const Vertex> vertices,
    std::span<const std::uint64_t> membershipWords,
    std::size_t memberCount,
    std::uint64_t setHash) {
    const std::uint32_t tag = slotTag(setHash);
    while (true) {
        const Slot& entry = table_[slot];
        if (entry.record == emptyRecord) {
            return std::nullopt;
        }
        if (entry.record != tombstoneRecord && entry.tag == tag &&
            recordEqualsMasked(
                entry.record, vertices, membershipWords, memberCount)) {
            const std::uint32_t record = entry.record;
            Clique result = cliqueAt(record);
            eraseRecord(slot, record);
            return result;
        }
        slot = (slot + 1U) & tableMask_;
    }
}

void MaximalCliqueIndex::prefetch(std::uint64_t setHash) const noexcept {
    __builtin_prefetch(
        table_.data() + (static_cast<std::size_t>(setHash) & tableMask_), 0, 1);
}

void MaximalCliqueIndex::prefetchForInsert(
    std::uint64_t setHash) const noexcept {
    __builtin_prefetch(
        table_.data() + (static_cast<std::size_t>(setHash) & tableMask_), 1, 1);
}

}
