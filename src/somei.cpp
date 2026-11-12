#include "somei.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace immce {
namespace {

using Clock = std::chrono::steady_clock;

[[nodiscard]] double millisecondsBetween(
    Clock::time_point begin,
    Clock::time_point end) noexcept {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

// A vertex's adjacency list is sorted, so its small neighbour list N^s (the
// neighbours below it) and big neighbour list N^b (those above it) are the
// two sides of its own position -- the SNL and BNL of [42]'s vertex
// structure.
[[nodiscard]] std::span<const AdjacencyEntry> bigNeighbours(
    std::span<const AdjacencyEntry> list,
    Vertex owner) {
    const auto first = std::upper_bound(
        list.begin(), list.end(), owner,
        [](Vertex value, const AdjacencyEntry& entry) {
            return value < entry.neighbor;
        });
    return {first, list.end()};
}

[[nodiscard]] std::span<const AdjacencyEntry> smallNeighbours(
    std::span<const AdjacencyEntry> list,
    Vertex owner) {
    const auto last = std::lower_bound(
        list.begin(), list.end(), owner,
        [](const AdjacencyEntry& entry, Vertex value) {
            return entry.neighbor < value;
        });
    return {list.begin(), last};
}

[[nodiscard]] std::span<const AdjacencyEntry>::iterator firstAtLeast(
    std::span<const AdjacencyEntry> list,
    Vertex value) {
    return std::lower_bound(
        list.begin(), list.end(), value,
        [](const AdjacencyEntry& entry, Vertex wanted) {
            return entry.neighbor < wanted;
        });
}

// The hash of a body is folded over its vertices in order, so appending a
// vertex extends it.
constexpr std::uint64_t emptyBodyHash = 0x84222325cbf29ce4ULL;

} // namespace

// ----------------------------------------------------------------------------
// ICTree child table

SomeiIndex::ChildTable::ChildTable() : slots_(1024U), mask_(1023U) {}

std::uint64_t SomeiIndex::ChildTable::key(
    NodeId parent,
    Vertex vertex) noexcept {
    return (std::uint64_t{parent} << 32U) |
        std::uint64_t{static_cast<std::uint32_t>(vertex)};
}

std::size_t SomeiIndex::ChildTable::home(std::uint64_t value) const noexcept {
    value ^= value >> 30U;
    value *= 0xbf58476d1ce4e5b9ULL;
    value ^= value >> 27U;
    value *= 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return static_cast<std::size_t>(value) & mask_;
}

SomeiIndex::NodeId SomeiIndex::ChildTable::find(
    NodeId parent,
    Vertex vertex) const noexcept {
    const std::uint64_t wanted = key(parent, vertex);
    for (std::size_t slot = home(wanted);; slot = (slot + 1U) & mask_) {
        if (slots_[slot].key == wanted) {
            return slots_[slot].child;
        }
        if (slots_[slot].key == emptyKey) {
            return noNode;
        }
    }
}

void SomeiIndex::ChildTable::insert(NodeId parent, Vertex vertex, NodeId child) {
    if ((size_ + 1U) * 10U > slots_.size() * 7U) {
        grow();
    }
    const std::uint64_t value = key(parent, vertex);
    std::size_t slot = home(value);
    while (slots_[slot].key != emptyKey) {
        slot = (slot + 1U) & mask_;
    }
    slots_[slot] = Slot{value, child};
    ++size_;
}

void SomeiIndex::ChildTable::erase(NodeId parent, Vertex vertex) {
    const std::uint64_t value = key(parent, vertex);
    std::size_t hole = home(value);
    while (slots_[hole].key != value) {
        if (slots_[hole].key == emptyKey) {
            throw std::logic_error("SOMEi: erasing an ICTree child that is not there");
        }
        hole = (hole + 1U) & mask_;
    }
    // Backward-shift deletion: pull later members of the probe run into the
    // hole unless their home slot lies cyclically in (hole, next].
    std::size_t next = hole;
    while (true) {
        next = (next + 1U) & mask_;
        if (slots_[next].key == emptyKey) {
            break;
        }
        const std::size_t ideal = home(slots_[next].key);
        const bool stays = hole <= next
            ? (hole < ideal && ideal <= next)
            : (hole < ideal || ideal <= next);
        if (!stays) {
            slots_[hole] = slots_[next];
            hole = next;
        }
    }
    slots_[hole] = Slot{};
    --size_;
}

void SomeiIndex::ChildTable::grow() {
    std::vector<Slot> old(slots_.size() * 2U);
    old.swap(slots_);
    mask_ = slots_.size() - 1U;
    for (const Slot& entry : old) {
        if (entry.key == emptyKey) {
            continue;
        }
        std::size_t slot = home(entry.key);
        while (slots_[slot].key != emptyKey) {
            slot = (slot + 1U) & mask_;
        }
        slots_[slot] = entry;
    }
}

// ----------------------------------------------------------------------------
// SOMEi

SomeiIndex::SomeiIndex(
    const IncrementalGraph& graph,
    SomeiConfiguration configuration)
    : graph_(graph),
      configuration_(configuration),
      headers_(static_cast<std::size_t>(graph.vertexCount())),
      nodes_(1U),
      rebuildFrom_(static_cast<std::size_t>(graph.vertexCount()), -1) {}

std::size_t SomeiIndex::treeNodeCount() const noexcept {
    return nodes_.size() - 1U - freeNodes_.size();
}

void SomeiIndex::applyBatch(
    std::span<const EdgeId> batch,
    SomeiStats& stats,
    std::vector<Clique>& newCliques,
    std::vector<Clique>& deletedCliques) {
    beginBatch(batch, stats);
    endBatch(stats, newCliques, deletedCliques);
}

void SomeiIndex::beginBatch(std::span<const EdgeId> batch, SomeiStats& stats) {
    ++stats.batches;
    stats.batchEdges += static_cast<long long>(batch.size());
    events_.clear();
    eventVertices_.clear();

    // Line 1: BNI <- GetUpdateScope(G, H).
    const auto scopeBegin = Clock::now();
    computeUpdateScope(batch, stats);
    stats.scopeMilliseconds += millisecondsBetween(scopeBegin, Clock::now());
}

void SomeiIndex::endBatch(
    SomeiStats& stats,
    std::vector<Clique>& newCliques,
    std::vector<Clique>& deletedCliques) {
    const auto rebuildBegin = Clock::now();
    // Lines 2-6, header by header.
    lastScope_.clear();
    for (const Vertex header : scopeHeaders_) {
        Vertex& entry = rebuildFrom_[static_cast<std::size_t>(header)];
        const Vertex from = entry;
        entry = -1;
        lastScope_.emplace_back(header, from);
        if (configuration_.profilePhases) {
            const auto begin = Clock::now();
            rollBack(header, from, stats);
            const auto rolledBack = Clock::now();
            redo(header, from, stats);
            const auto redone = Clock::now();
            addAllToTree(header, stats);
            const auto added = Clock::now();
            stats.rollbackMilliseconds += millisecondsBetween(begin, rolledBack);
            stats.redoMilliseconds += millisecondsBetween(rolledBack, redone);
            stats.treeAddMilliseconds += millisecondsBetween(redone, added);
        } else {
            rollBack(header, from, stats);
            redo(header, from, stats);
            addAllToTree(header, stats);
        }
    }
    const auto changeBegin = Clock::now();
    stats.rebuildMilliseconds += millisecondsBetween(rebuildBegin, changeBegin);

    if (configuration_.trackChanges) {
        netChanges(newCliques, deletedCliques, stats);
    }
    stats.changeMilliseconds +=
        millisecondsBetween(changeBegin, Clock::now());
}

void SomeiIndex::computeUpdateScope(
    std::span<const EdgeId> batch,
    SomeiStats& stats) {
    // Lines 2-5: every change e(v_x, v_y), v_x < v_y, under its big endpoint.
    // Sorted, the changes at one v_y are a run whose last small endpoint is
    // SNI[v_y] and whose small endpoints are updateNeighborMap[v_y]; SNI is
    // walked in vertex order, as an ordered map would be.
    changes_.clear();
    changes_.reserve(batch.size());
    for (const EdgeId edgeId : batch) {
        const Edge& edge = graph_.edge(edgeId);
        changes_.emplace_back(edge.v, edge.u);
    }
    std::ranges::sort(changes_);
    scopeHeaders_.clear();

    std::size_t begin = 0U;
    while (begin < changes_.size()) {
        const Vertex big = changes_[begin].first;
        std::size_t end = begin + 1U;
        while (end < changes_.size() && changes_[end].first == big) {
            ++end;
        }
        const Vertex boundary = changes_[end - 1U].second;
        changedSmall_.clear();
        for (std::size_t index = begin; index < end; ++index) {
            changedSmall_.push_back(changes_[index].second);
        }
        ++stats.bigVertices;

        // Lines 6-9: each small neighbour v_k <= SNI[v_y] of v_y takes v_y
        // as its rollback vertex unless it already has an earlier one.
        for (const AdjacencyEntry& entry :
             smallNeighbours(graph_.neighbors(big), big)) {
            const Vertex candidate = entry.neighbor;
            if (candidate > boundary) {
                break;
            }
            ++stats.scopeVisits;
            Vertex& from = rebuildFrom_[static_cast<std::size_t>(candidate)];
            if (from != -1 && from <= big) {
                continue;
            }
            if (configuration_.pruning &&
                !touchesHeader(candidate, changedSmall_)) {
                ++stats.prunedVisits;
                continue;
            }
            if (from == -1) {
                scopeHeaders_.push_back(candidate);
            }
            from = big;
        }
        begin = end;
    }
    // Algorithm 3 walks BNI's keys; in vertex order here.
    std::ranges::sort(scopeHeaders_);
    stats.rebuiltHeaders += static_cast<long long>(scopeHeaders_.size());
}

bool SomeiIndex::touchesHeader(
    Vertex header,
    std::span<const Vertex> changed) const {
    for (const Vertex smallEndpoint : changed) {
        if (smallEndpoint == header) {
            return true;
        }
        if (smallEndpoint > header && graph_.hasEdge(header, smallEndpoint)) {
            return true;
        }
    }
    return false;
}

void SomeiIndex::rollBack(Vertex header, Vertex from, SomeiStats& stats) {
    std::vector<Body>& bodies = headers_[static_cast<std::size_t>(header)];
    const auto count = static_cast<std::uint32_t>(bodies.size());
    for (std::uint32_t index = 0; index < count; ++index) {
        ++stats.rolledBackBodies;
        // Lines 2-3: delete the link of l and its branch.
        treeDelete(header, index, stats);
        // Lines 6-7: a body created by T_{j-1} goes back to its version at
        // T_{j-1}, its vertices below v_j.  One created later (lines 4-5) is
        // removed below, once every body of the header is unlinked: until
        // then a later body's index is what its marked ICTree node points at.
        Body& body = bodies[index];
        if (body.createdAt >= from) {
            continue;
        }
        const auto cut = std::ranges::lower_bound(body.vertices, from);
        if (cut != body.vertices.end()) {
            ++stats.truncatedBodies;
            body.vertices.erase(cut, body.vertices.end());
        }
    }
    const std::size_t removed = std::erase_if(
        bodies, [from](const Body& body) { return body.createdAt >= from; });
    stats.removedBodies += static_cast<long long>(removed);
    bodyCount_ -= removed;
}

void SomeiIndex::redo(Vertex header, Vertex from, SomeiStats& stats) {
    // L <- Gamma_r(v_j, N^b_{v_i}) u {v_j}: the big neighbours from v_j on,
    // popped from the front.
    const std::span<const AdjacencyEntry> big =
        bigNeighbours(graph_.neighbors(header), header);
    // The untimed build's steps by enumeration (see somei.hpp).
    const bool enumerate = building_ && configuration_.enumerateBuild;
    if (enumerate) {
        indexBodies(header);
    }
    for (auto position = firstAtLeast(big, from); position != big.end();
         ++position) {
        const Vertex renovator = position->neighbor;
        ++stats.redoSteps;
        // Line 3: R <- N^b_{v_i} n N^s_{v_r}, merged as [42] does with the
        // part of N^s_{v_r} after v_i (its "left list"; v_i is in N^s_{v_r}
        // because v_r is in N^b_{v_i}).
        const std::span<const AdjacencyEntry> small =
            smallNeighbours(graph_.neighbors(renovator), renovator);
        auto left = std::upper_bound(
            small.begin(), small.end(), header,
            [](Vertex value, const AdjacencyEntry& entry) {
                return value < entry.neighbor;
            });
        auto right = big.begin();
        renovation_.clear();
        while (left != small.end() && right != big.end()) {
            if (left->neighbor < right->neighbor) {
                ++left;
            } else if (right->neighbor < left->neighbor) {
                ++right;
            } else {
                renovation_.push_back(left->neighbor);
                ++left;
                ++right;
            }
        }
        if (enumerate) {
            updateCliqueByEnumeration(header, renovator, stats);
        } else {
            updateClique(header, renovator, stats);
        }
    }
    if (enumerate) {
        releaseBodyIndex();
    }
}

void SomeiIndex::updateClique(
    Vertex header,
    Vertex renovator,
    SomeiStats& stats) {
    std::vector<Body>& bodies = headers_[static_cast<std::size_t>(header)];
    // Lines 4-5.
    if (bodies.empty() || renovation_.empty()) {
        ++stats.emptyRenovations;
        Body created;
        created.vertices.push_back(renovator);
        created.createdAt = renovator;
        bodies.push_back(std::move(created));
        ++stats.createdBodies;
        ++bodyCount_;
        return;
    }

    // Lines 7-17.  One two-way merge per body gives l_j n R, and its size
    // against |l_j| and |R| decides which of the four cases holds.
    partialVertices_.clear();
    partialRanges_.clear();
    fullVertices_.clear();
    fullRanges_.clear();
    bool matched = false;
    const std::size_t existing = bodies.size();
    for (std::size_t index = 0; index < existing; ++index) {
        ++stats.bodyComparisons;
        std::vector<Vertex>& body = bodies[index].vertices;
        intersection_.clear();
        std::ranges::set_intersection(
            body, renovation_, std::back_inserter(intersection_));
        const std::size_t common = intersection_.size();
        if (common == body.size() && common == renovation_.size()) {
            // Lines 9-10: l_j == R.
            body.push_back(renovator);
            ++stats.extendedBodies;
            matched = true;
            break;
        }
        if (common == renovation_.size()) {
            // Lines 11-13: R is a proper subset of l_j.
            Body created;
            created.vertices.reserve(renovation_.size() + 1U);
            created.vertices.assign(renovation_.begin(), renovation_.end());
            created.vertices.push_back(renovator);
            created.createdAt = renovator;
            bodies.push_back(std::move(created));
            ++stats.createdBodies;
            ++bodyCount_;
            matched = true;
            break;
        }
        if (common == body.size()) {
            // Lines 14-15: l_j is a proper subset of R.
            body.push_back(renovator);
            ++stats.extendedBodies;
            fullRanges_.emplace_back(fullVertices_.size(), common);
            fullVertices_.insert(
                fullVertices_.end(), intersection_.begin(), intersection_.end());
        }
        if (common != 0U) {
            // Lines 16-17.
            partialRanges_.emplace_back(partialVertices_.size(), common);
            partialVertices_.insert(
                partialVertices_.end(), intersection_.begin(),
                intersection_.end());
        }
    }
    if (matched) {
        return;
    }

    // Lines 18-22.  pSet is a set of vertex lists, so a repeated
    // intersection counts once; a list then survives unless it is contained
    // in a list of fSet or in another list of pSet.
    const auto partial = [this](const std::pair<std::size_t, std::size_t>& range) {
        return std::span<const Vertex>(
            partialVertices_.data() + range.first, range.second);
    };
    const auto full = [this](const std::pair<std::size_t, std::size_t>& range) {
        return std::span<const Vertex>(
            fullVertices_.data() + range.first, range.second);
    };
    std::ranges::sort(partialRanges_, [&](const auto& lhs, const auto& rhs) {
        return std::ranges::lexicographical_compare(partial(lhs), partial(rhs));
    });
    partialRanges_.erase(
        std::unique(
            partialRanges_.begin(), partialRanges_.end(),
            [&](const auto& lhs, const auto& rhs) {
                return std::ranges::equal(partial(lhs), partial(rhs));
            }),
        partialRanges_.end());
    stats.partialIntersections += static_cast<long long>(partialRanges_.size());
    for (std::size_t index = 0; index < partialRanges_.size(); ++index) {
        const std::span<const Vertex> candidate = partial(partialRanges_[index]);
        bool contained = false;
        for (std::size_t other = 0; other < fullRanges_.size() && !contained;
             ++other) {
            contained = std::ranges::includes(full(fullRanges_[other]), candidate);
        }
        for (std::size_t other = 0;
             other < partialRanges_.size() && !contained; ++other) {
            contained = other != index &&
                std::ranges::includes(partial(partialRanges_[other]), candidate);
        }
        if (contained) {
            continue;
        }
        Body created;
        created.vertices.reserve(candidate.size() + 1U);
        created.vertices.assign(candidate.begin(), candidate.end());
        created.vertices.push_back(renovator);
        created.createdAt = renovator;
        bodies.push_back(std::move(created));
        ++stats.createdBodies;
        ++bodyCount_;
    }
}

std::uint64_t SomeiIndex::extendHash(std::uint64_t hash, Vertex vertex) noexcept {
    hash = (hash + static_cast<std::uint64_t>(vertex) + 1U) * 0x9e3779b97f4a7c15ULL;
    return hash ^ (hash >> 29U);
}

void SomeiIndex::releaseBodyIndex() {
    // clear() keeps every bucket and zeroes them all, so after a header with
    // a hundred thousand bodies it would tax every header that follows.
    if (bodyByHash_.bucket_count() > 1024U) {
        bodyByHash_ = decltype(bodyByHash_){};
    } else {
        bodyByHash_.clear();
    }
}

void SomeiIndex::indexBodies(Vertex header) {
    releaseBodyIndex();
    const std::vector<Body>& bodies = headers_[static_cast<std::size_t>(header)];
    bodyByHash_.reserve(bodies.size());
    for (std::uint32_t index = 0;
         index < static_cast<std::uint32_t>(bodies.size()); ++index) {
        std::uint64_t hash = emptyBodyHash;
        for (const Vertex vertex : bodies[index].vertices) {
            hash = extendHash(hash, vertex);
        }
        bodyByHash_.emplace(hash, index);
    }
}

void SomeiIndex::enumerateRenovation(std::size_t depth) {
    const std::size_t words = renovationWords_;
    std::uint64_t* const candidates = enumerationStack_.data() + 2U * depth * words;
    std::uint64_t* const excluded = candidates + words;
    bool candidatesLeft = false;
    bool excludedLeft = false;
    for (std::size_t word = 0; word < words; ++word) {
        candidatesLeft = candidatesLeft || candidates[word] != 0U;
        excludedLeft = excludedLeft || excluded[word] != 0U;
    }
    if (!candidatesLeft) {
        if (!excludedLeft) {
            const std::size_t offset = foundPositions_.size();
            foundPositions_.insert(
                foundPositions_.end(), enumerationClique_.begin(),
                enumerationClique_.end());
            std::sort(foundPositions_.begin() + static_cast<std::ptrdiff_t>(offset),
                      foundPositions_.end());
            foundRanges_.emplace_back(offset, enumerationClique_.size());
        }
        return;
    }
    // Tomita's pivot: the vertex of P u X with the most neighbours in P.
    std::size_t pivot = 0;
    long best = -1;
    for (std::size_t word = 0; word < words; ++word) {
        std::uint64_t pool = candidates[word] | excluded[word];
        while (pool != 0U) {
            const std::size_t vertex =
                word * 64U + static_cast<std::size_t>(std::countr_zero(pool));
            pool &= pool - 1U;
            const std::uint64_t* const row = renovationRows_.data() + vertex * words;
            long count = 0;
            for (std::size_t other = 0; other < words; ++other) {
                count += std::popcount(candidates[other] & row[other]);
            }
            if (count > best) {
                best = count;
                pivot = vertex;
            }
        }
    }
    const std::uint64_t* const pivotRow = renovationRows_.data() + pivot * words;
    std::uint64_t* const nextCandidates = candidates + 2U * words;
    std::uint64_t* const nextExcluded = nextCandidates + words;
    for (std::size_t word = 0; word < words; ++word) {
        std::uint64_t branches = candidates[word] & ~pivotRow[word];
        while (branches != 0U) {
            const auto bit = static_cast<unsigned>(std::countr_zero(branches));
            branches &= branches - 1U;
            const std::size_t vertex = word * 64U + bit;
            const std::uint64_t* const row = renovationRows_.data() + vertex * words;
            for (std::size_t other = 0; other < words; ++other) {
                nextCandidates[other] = candidates[other] & row[other];
                nextExcluded[other] = excluded[other] & row[other];
            }
            enumerationClique_.push_back(static_cast<std::uint32_t>(vertex));
            enumerateRenovation(depth + 1U);
            enumerationClique_.pop_back();
            candidates[word] &= ~(std::uint64_t{1} << bit);
            excluded[word] |= std::uint64_t{1} << bit;
        }
    }
}

void SomeiIndex::updateCliqueByEnumeration(
    Vertex header,
    Vertex renovator,
    SomeiStats& stats) {
    std::vector<Body>& bodies = headers_[static_cast<std::size_t>(header)];
    const std::size_t size = renovation_.size();
    // Lines 4-5.
    if (bodies.empty() || size == 0U) {
        ++stats.emptyRenovations;
        bodyByHash_.emplace(
            extendHash(emptyBodyHash, renovator),
            static_cast<std::uint32_t>(bodies.size()));
        Body created;
        created.vertices.push_back(renovator);
        created.createdAt = renovator;
        bodies.push_back(std::move(created));
        ++stats.createdBodies;
        ++bodyCount_;
        return;
    }
    if (size > configuration_.enumerateBuildLimit) {
        updateClique(header, renovator, stats);
        indexBodies(header);
        return;
    }

    // G[R] as bit rows over R's positions (R is ascending).
    const std::size_t words = (size + 63U) / 64U;
    renovationWords_ = words;
    renovationRows_.assign(size * words, 0U);
    const auto connect = [&](std::size_t first, std::size_t second) {
        renovationRows_[first * words + second / 64U] |=
            std::uint64_t{1} << (second % 64U);
        renovationRows_[second * words + first / 64U] |=
            std::uint64_t{1} << (first % 64U);
    };
    const auto before = [](const AdjacencyEntry& entry, Vertex value) {
        return entry.neighbor < value;
    };
    for (std::size_t first = 0; first + 1U < size; ++first) {
        const std::span<const AdjacencyEntry> list =
            graph_.neighbors(renovation_[first]);
        auto position = std::lower_bound(
            list.begin(), list.end(), renovation_[first + 1U], before);
        const auto end = std::lower_bound(
            position, list.end(), renovation_.back() + 1, before);
        // The later members of R in this vertex's list: merged when the list
        // is short there, binary-searched when it is long.
        if (static_cast<std::size_t>(end - position) <= 8U * (size - first - 1U)) {
            std::size_t second = first + 1U;
            while (position != end && second < size) {
                if (position->neighbor < renovation_[second]) {
                    ++position;
                } else if (renovation_[second] < position->neighbor) {
                    ++second;
                } else {
                    connect(first, second);
                    ++position;
                    ++second;
                }
            }
        } else {
            for (std::size_t second = first + 1U; second < size; ++second) {
                position = std::lower_bound(position, end, renovation_[second], before);
                if (position == end) {
                    break;
                }
                if (position->neighbor == renovation_[second]) {
                    connect(first, second);
                }
            }
        }
    }

    // The maximal cliques of G[R]: P every position, X empty.
    const std::size_t stackWords = 2U * words * (size + 2U);
    if (enumerationStack_.size() < stackWords) {
        enumerationStack_.resize(stackWords);
    }
    std::fill(enumerationStack_.begin(),
              enumerationStack_.begin() + static_cast<std::ptrdiff_t>(2U * words),
              std::uint64_t{0});
    for (std::size_t position = 0; position < size; ++position) {
        enumerationStack_[position / 64U] |= std::uint64_t{1} << (position % 64U);
    }
    enumerationClique_.clear();
    foundPositions_.clear();
    foundRanges_.clear();
    enumerateRenovation(0U);

    // A clique that is a body -- one that lies in R -- is extended (lines
    // 9-10 and 14-15); any other, K, becomes the body K u {v_r} (lines 11-13
    // and 18-22).
    toCreate_.clear();
    for (std::uint32_t found = 0;
         found < static_cast<std::uint32_t>(foundRanges_.size()); ++found) {
        const auto [offset, length] = foundRanges_[found];
        std::uint64_t hash = emptyBodyHash;
        for (std::size_t index = offset; index < offset + length; ++index) {
            hash = extendHash(hash, renovation_[foundPositions_[index]]);
        }
        bool extended = false;
        const auto [first, last] = bodyByHash_.equal_range(hash);
        for (auto entry = first; entry != last; ++entry) {
            std::vector<Vertex>& body = bodies[entry->second].vertices;
            bool same = body.size() == length;
            for (std::size_t index = 0; index < length && same; ++index) {
                same = body[index] == renovation_[foundPositions_[offset + index]];
            }
            if (!same) {
                continue;
            }
            const std::uint32_t bodyIndex = entry->second;
            bodyByHash_.erase(entry);
            body.push_back(renovator);
            bodyByHash_.emplace(extendHash(hash, renovator), bodyIndex);
            ++stats.extendedBodies;
            extended = true;
            break;
        }
        if (!extended) {
            toCreate_.push_back(found);
        }
    }
    // Created in the order the paper's lines 18-22 create them, pSet's
    // lexicographic one; positions in R compare as the vertices do.
    std::ranges::sort(toCreate_, [&](std::uint32_t lhs, std::uint32_t rhs) {
        const auto [lhsOffset, lhsLength] = foundRanges_[lhs];
        const auto [rhsOffset, rhsLength] = foundRanges_[rhs];
        return std::lexicographical_compare(
            foundPositions_.begin() + static_cast<std::ptrdiff_t>(lhsOffset),
            foundPositions_.begin() + static_cast<std::ptrdiff_t>(lhsOffset + lhsLength),
            foundPositions_.begin() + static_cast<std::ptrdiff_t>(rhsOffset),
            foundPositions_.begin() + static_cast<std::ptrdiff_t>(rhsOffset + rhsLength));
    });
    for (const std::uint32_t found : toCreate_) {
        const auto [offset, length] = foundRanges_[found];
        Body created;
        created.vertices.reserve(length + 1U);
        std::uint64_t hash = emptyBodyHash;
        for (std::size_t index = offset; index < offset + length; ++index) {
            created.vertices.push_back(renovation_[foundPositions_[index]]);
            hash = extendHash(hash, created.vertices.back());
        }
        created.vertices.push_back(renovator);
        created.createdAt = renovator;
        bodyByHash_.emplace(
            extendHash(hash, renovator), static_cast<std::uint32_t>(bodies.size()));
        bodies.push_back(std::move(created));
        ++stats.createdBodies;
        ++bodyCount_;
    }
}

void SomeiIndex::addAllToTree(Vertex header, SomeiStats& stats) {
    const auto count = static_cast<std::uint32_t>(
        headers_[static_cast<std::size_t>(header)].size());
    for (std::uint32_t index = 0; index < count; ++index) {
        treeAdd(header, index, stats);
    }
}

SomeiIndex::NodeId SomeiIndex::newChild(
    NodeId parent,
    Vertex vertex,
    SomeiStats& stats) {
    // A marked leaf that gains a child stops being a maximal clique.
    if (parent != rootNode && nodes_[parent].children == 0U &&
        nodes_[parent].marked) {
        recordLeaf(parent, -1);
    }
    NodeId node = noNode;
    if (!freeNodes_.empty()) {
        node = freeNodes_.back();
        freeNodes_.pop_back();
    } else {
        if (nodes_.size() >= static_cast<std::size_t>(noNode)) {
            throw std::length_error("SOMEi: too many ICTree nodes");
        }
        node = static_cast<NodeId>(nodes_.size());
        nodes_.emplace_back();
    }
    nodes_[node] = TreeNode{vertex, parent, 0U, 0U, false};
    ++nodes_[parent].children;
    children_.insert(parent, vertex, node);
    ++stats.treeNodesCreated;
    return node;
}

void SomeiIndex::deleteNode(NodeId node, SomeiStats& stats) {
    const TreeNode removed = nodes_[node];
    children_.erase(removed.parent, removed.vertex);
    --nodes_[removed.parent].children;
    nodes_[node] = TreeNode{};
    freeNodes_.push_back(node);
    ++stats.treeNodesDeleted;
}

void SomeiIndex::treeAdd(
    Vertex header,
    std::uint32_t bodyIndex,
    SomeiStats& stats) {
    ++stats.treeAdds;
    Body& body = headers_[static_cast<std::size_t>(header)][bodyIndex];
    // Lines 1-6: down from the root along the body in reverse.
    NodeId node = rootNode;
    for (auto vertex = body.vertices.rbegin(); vertex != body.vertices.rend();
         ++vertex) {
        NodeId child = children_.find(node, *vertex);
        if (child == noNode) {
            child = newChild(node, *vertex, stats);
        }
        node = child;
    }
    // Lines 7-9: the header, marked and linked.
    NodeId target = children_.find(node, header);
    if (target == noNode) {
        target = newChild(node, header, stats);
    }
    TreeNode& linked = nodes_[target];
    if (linked.marked) {
        throw std::logic_error("SOMEi: a v-header clique was added twice");
    }
    linked.marked = true;
    linked.body = bodyIndex;
    body.link = target;
    if (linked.children == 0U) {
        recordLeaf(header, body, +1);
    }
}

void SomeiIndex::treeDelete(
    Vertex header,
    std::uint32_t bodyIndex,
    SomeiStats& stats) {
    ++stats.treeDeletes;
    Body& body = headers_[static_cast<std::size_t>(header)][bodyIndex];
    NodeId node = body.link;
    if (node == noNode) {
        throw std::logic_error("SOMEi: deleting a clique that has no ICTree link");
    }
    body.link = noNode;
    nodes_[node].marked = false;
    // Lines 2-3: a node with children only loses its mark.
    if (nodes_[node].children != 0U) {
        return;
    }
    // Lines 5-9: the branch it ended was a maximal clique; delete upwards
    // while the nodes are unmarked and childless.
    recordLeaf(header, body, -1);
    NodeId parent = nodes_[node].parent;
    deleteNode(node, stats);
    node = parent;
    while (node != rootNode && !nodes_[node].marked &&
           nodes_[node].children == 0U) {
        parent = nodes_[node].parent;
        deleteNode(node, stats);
        node = parent;
    }
    // A marked node left without children is now a leaf: its v-header clique
    // has become a maximal clique.
    if (node != rootNode && nodes_[node].marked && nodes_[node].children == 0U) {
        recordLeaf(node, +1);
    }
}

void SomeiIndex::recordLeaf(NodeId node, int sign) {
    const TreeNode& leaf = nodes_[node];
    recordLeaf(
        leaf.vertex,
        headers_[static_cast<std::size_t>(leaf.vertex)][leaf.body],
        sign);
}

void SomeiIndex::recordLeaf(Vertex header, const Body& body, int sign) {
    if (sign > 0) {
        ++leafCount_;
    } else {
        --leafCount_;
    }
    if (!configuration_.trackChanges) {
        return;
    }
    // The header is the clique's smallest vertex, so header + body is the
    // clique in ascending order.
    events_.push_back(LeafEvent{
        eventVertices_.size(),
        static_cast<std::uint32_t>(body.vertices.size() + 1U),
        sign});
    eventVertices_.push_back(header);
    eventVertices_.insert(
        eventVertices_.end(), body.vertices.begin(), body.vertices.end());
}

void SomeiIndex::netChanges(
    std::vector<Clique>& newCliques,
    std::vector<Clique>& deletedCliques,
    SomeiStats& stats) {
    stats.leafEvents += static_cast<long long>(events_.size());
    const auto clique = [this](std::uint32_t index) {
        const LeafEvent& event = events_[index];
        return std::span<const Vertex>(
            eventVertices_.data() + event.offset, event.length);
    };
    eventOrder_.resize(events_.size());
    std::iota(eventOrder_.begin(), eventOrder_.end(), std::uint32_t{0});
    std::ranges::sort(eventOrder_, [&](std::uint32_t lhs, std::uint32_t rhs) {
        return std::ranges::lexicographical_compare(clique(lhs), clique(rhs));
    });
    // A clique's leaf comes and goes in turn, so the signs of one clique sum
    // to its net change: +1 a new maximal clique, -1 a deleted one.
    std::size_t begin = 0U;
    while (begin < eventOrder_.size()) {
        const std::span<const Vertex> members = clique(eventOrder_[begin]);
        int net = events_[eventOrder_[begin]].sign;
        std::size_t end = begin + 1U;
        while (end < eventOrder_.size() &&
               std::ranges::equal(members, clique(eventOrder_[end]))) {
            net += events_[eventOrder_[end]].sign;
            ++end;
        }
        if (net == 1) {
            newCliques.emplace_back(members.begin(), members.end());
            ++stats.newCliques;
        } else if (net == -1) {
            deletedCliques.emplace_back(members.begin(), members.end());
            ++stats.deletedCliques;
        } else if (net != 0) {
            throw std::logic_error("SOMEi: inconsistent ICTree leaf changes");
        }
        begin = end;
    }
}

std::vector<SomeiIndex::BodyView> SomeiIndex::bodies(Vertex header) const {
    std::vector<BodyView> result;
    for (const Body& body : headers_.at(static_cast<std::size_t>(header))) {
        result.push_back(BodyView{body.vertices, body.createdAt});
    }
    return result;
}

std::vector<Clique> SomeiIndex::maximalCliques() const {
    std::vector<Clique> result;
    for (std::size_t index = 1U; index < nodes_.size(); ++index) {
        const TreeNode& node = nodes_[index];
        if (node.vertex < 0 || node.children != 0U) {
            continue;
        }
        if (!node.marked) {
            throw std::logic_error("SOMEi: an unmarked ICTree leaf");
        }
        Clique clique;
        for (NodeId walk = static_cast<NodeId>(index); walk != rootNode;
             walk = nodes_[walk].parent) {
            clique.push_back(nodes_[walk].vertex);
        }
        std::ranges::sort(clique);
        result.push_back(std::move(clique));
    }
    return result;
}

} // namespace immce
