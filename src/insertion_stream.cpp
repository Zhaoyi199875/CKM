#include "insertion_stream.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ckm {
namespace {

struct Row {
    int a = 0;
    int b = 0;
    int timestamp = 0;
    std::uint32_t line = 0U;
};

constexpr std::size_t readBufferBytes = 1U << 20U;
constexpr std::size_t idChunkLimit = 1U << 22U;

class FileScanner {
public:
    explicit FileScanner(const std::string& path)
        : file_(std::fopen(path.c_str(), "rb")), buffer_(readBufferBytes) {
        if (file_ == nullptr) {
            throw std::runtime_error("cannot open input file: " + path);
        }
    }

    ~FileScanner() {
        if (file_ != nullptr) {
            std::fclose(file_);
        }
    }

    FileScanner(const FileScanner&) = delete;
    FileScanner& operator=(const FileScanner&) = delete;

    [[nodiscard]] int peek() {
        if (position_ == filled_ && !refill()) {
            return -1;
        }
        return static_cast<unsigned char>(buffer_[position_]);
    }

    void advance() noexcept { ++position_; }

private:
    [[nodiscard]] bool refill() {
        filled_ = std::fread(buffer_.data(), 1U, buffer_.size(), file_);
        position_ = 0U;
        return filled_ != 0U;
    }

    std::FILE* file_ = nullptr;
    std::vector<char> buffer_;
    std::size_t position_ = 0U;
    std::size_t filled_ = 0U;
};

[[nodiscard]] bool readRow(FileScanner& scanner, Row& row, int& lineNumber) {
    while (true) {
        int character = scanner.peek();
        while (character == ' ' || character == '\t' || character == '\r' ||
               character == '\n' || character == '#') {
            if (character == '\n') {
                ++lineNumber;
            }
            if (character == '#') {
                while (character != -1 && character != '\n') {
                    scanner.advance();
                    character = scanner.peek();
                }
                continue;
            }
            scanner.advance();
            character = scanner.peek();
        }
        if (character == -1) {
            return false;
        }
        ++lineNumber;

        int* const fields[3] = {&row.a, &row.b, &row.timestamp};
        for (std::size_t index = 0; index < 3U; ++index) {
            character = scanner.peek();
            while (character == ' ' || character == '\t') {
                scanner.advance();
                character = scanner.peek();
            }
            bool negative = false;
            if (character == '-' || character == '+') {
                negative = character == '-';
                scanner.advance();
                character = scanner.peek();
            }
            if (character < '0' || character > '9') {
                throw std::runtime_error(
                    "line " + std::to_string(lineNumber) +
                    " must contain exactly 3 integer columns: u v timestamp");
            }
            long long value = 0;
            while (character >= '0' && character <= '9') {
                value = value * 10 + (character - '0');
                if (value > std::numeric_limits<int>::max()) {
                    throw std::runtime_error(
                        "line " + std::to_string(lineNumber) +
                        " has a column that does not fit in a 32-bit integer");
                }
                scanner.advance();
                character = scanner.peek();
            }
            *fields[index] = static_cast<int>(negative ? -value : value);
        }

        character = scanner.peek();
        while (character == ' ' || character == '\t' || character == '\r') {
            scanner.advance();
            character = scanner.peek();
        }
        if (character == '#') {
            while (character != -1 && character != '\n') {
                scanner.advance();
                character = scanner.peek();
            }
        }
        if (character != -1 && character != '\n') {
            throw std::runtime_error(
                "line " + std::to_string(lineNumber) +
                " must contain exactly 3 integer columns: u v timestamp");
        }
        return true;
    }
}

void foldUniqueIds(std::vector<int>& pending, std::vector<int>& unique) {
    std::ranges::sort(pending);
    pending.erase(std::unique(pending.begin(), pending.end()), pending.end());
    if (unique.empty()) {
        unique.swap(pending);
        pending.clear();
        return;
    }
    std::vector<int> merged;
    merged.reserve(unique.size() + pending.size());
    std::ranges::set_union(unique, pending, std::back_inserter(merged));
    unique.swap(merged);
    pending.clear();
}

}

InsertionStream InsertionStream::readEdgeList(const std::string& path) {
    FileScanner scanner(path);
    std::vector<Row> rows;
    std::vector<int> pendingIds;
    std::vector<int> uniqueIds;
    pendingIds.reserve(idChunkLimit);

    Row row;
    int lineNumber = 0;
    while (readRow(scanner, row, lineNumber)) {
        if (row.a == row.b) {
            continue;
        }
        if (rows.size() >= std::numeric_limits<std::uint32_t>::max()) {
            throw std::length_error("too many input rows for 32-bit positions");
        }
        row.line = static_cast<std::uint32_t>(rows.size());
        rows.push_back(row);
        pendingIds.push_back(row.a);
        pendingIds.push_back(row.b);
        if (pendingIds.size() >= idChunkLimit) {
            foldUniqueIds(pendingIds, uniqueIds);
        }
    }
    if (!pendingIds.empty()) {
        foldUniqueIds(pendingIds, uniqueIds);
    }
    pendingIds.shrink_to_fit();
    if (rows.empty()) {
        throw std::runtime_error("input graph has no non-self edges");
    }
    if (uniqueIds.size() >
        static_cast<std::size_t>(std::numeric_limits<Vertex>::max())) {
        throw std::length_error("too many vertices for 32-bit dense IDs");
    }

    for (Row& entry : rows) {
        const auto denseA = static_cast<Vertex>(
            std::ranges::lower_bound(uniqueIds, entry.a) - uniqueIds.begin());
        const auto denseB = static_cast<Vertex>(
            std::ranges::lower_bound(uniqueIds, entry.b) - uniqueIds.begin());
        const Edge edge = normalizedEdge(denseA, denseB);
        entry.a = edge.u;
        entry.b = edge.v;
    }

    std::ranges::sort(rows, [](const Row& lhs, const Row& rhs) {
        if (lhs.a != rhs.a) {
            return lhs.a < rhs.a;
        }
        if (lhs.b != rhs.b) {
            return lhs.b < rhs.b;
        }
        if (lhs.timestamp != rhs.timestamp) {
            return lhs.timestamp < rhs.timestamp;
        }
        return lhs.line < rhs.line;
    });
    rows.erase(
        std::unique(
            rows.begin(),
            rows.end(),
            [](const Row& lhs, const Row& rhs) {
                return lhs.a == rhs.a && lhs.b == rhs.b;
            }),
        rows.end());
    rows.shrink_to_fit();
    if (rows.size() >
        static_cast<std::size_t>(std::numeric_limits<EdgeId>::max())) {
        throw std::length_error("too many edges for 32-bit edge IDs");
    }

    std::ranges::sort(rows, [](const Row& lhs, const Row& rhs) {
        if (lhs.timestamp != rhs.timestamp) {
            return lhs.timestamp < rhs.timestamp;
        }
        return lhs.line < rhs.line;
    });

    InsertionStream stream;
    stream.originalVertexIds_.swap(uniqueIds);
    stream.edges_.reserve(rows.size());
    stream.insertionTimestamps_.reserve(rows.size());
    for (const Row& entry : rows) {
        stream.edges_.push_back(Edge{entry.a, entry.b});
        stream.insertionTimestamps_.push_back(entry.timestamp);
    }
    return stream;
}

int InsertionStream::vertexCount() const noexcept {
    return static_cast<int>(originalVertexIds_.size());
}

int InsertionStream::edgeCount() const noexcept {
    return static_cast<int>(edges_.size());
}

const Edge& InsertionStream::edge(EdgeId edgeId) const {
    return edges_.at(static_cast<std::size_t>(edgeId));
}

std::span<const Edge> InsertionStream::edges() const noexcept {
    return edges_;
}

std::vector<Edge> InsertionStream::takeEdges() noexcept {
    return std::move(edges_);
}

int InsertionStream::originalVertexId(Vertex denseVertex) const {
    return originalVertexIds_.at(static_cast<std::size_t>(denseVertex));
}

std::span<const int> InsertionStream::insertionTimestamps() const noexcept {
    return insertionTimestamps_;
}

std::size_t InsertionStream::batchCount(
    std::optional<std::size_t> batchSize) const {
    if (batchSize.has_value()) {
        if (*batchSize == 0U) {
            throw std::invalid_argument("batch size must be positive");
        }
        return (insertionTimestamps_.size() + *batchSize - 1U) / *batchSize;
    }
    std::size_t count = 0U;
    for (std::size_t index = 0; index < insertionTimestamps_.size(); ++index) {
        if (index == 0U ||
            insertionTimestamps_[index] != insertionTimestamps_[index - 1U]) {
            ++count;
        }
    }
    return count;
}

std::vector<EdgeBatch> InsertionStream::makeBatches(
    std::optional<std::size_t> batchSize,
    std::optional<std::size_t> maxBatches,
    std::size_t firstEdge,
    std::size_t endEdge) const {
    if (batchSize.has_value() && *batchSize == 0U) {
        throw std::invalid_argument("batch size must be positive");
    }
    const bool wholeStream =
        firstEdge == 0U && endEdge == std::numeric_limits<std::size_t>::max();
    if (!batchSize.has_value() && !wholeStream) {
        throw std::invalid_argument("an edge range needs a fixed batch size");
    }
    const std::size_t limit =
        maxBatches.value_or(std::numeric_limits<std::size_t>::max());

    std::vector<EdgeBatch> batches;
    if (batchSize.has_value()) {
        const std::size_t stop = std::min(insertionTimestamps_.size(), endEdge);
        for (std::size_t begin = firstEdge; begin < stop;
             begin += *batchSize) {
            if (batches.size() >= limit) {
                break;
            }
            const std::size_t end = std::min(stop, begin + *batchSize);
            EdgeBatch batch;
            batch.timestamp = insertionTimestamps_[end - 1U];
            batch.edgeIds.reserve(end - begin);
            for (std::size_t index = begin; index < end; ++index) {
                batch.edgeIds.push_back(static_cast<EdgeId>(index));
            }
            batches.push_back(std::move(batch));
        }
        return batches;
    }

    for (std::size_t index = 0; index < insertionTimestamps_.size(); ++index) {
        const int timestamp = insertionTimestamps_[index];
        if (batches.empty() || batches.back().timestamp != timestamp) {
            if (batches.size() >= limit) {
                break;
            }
            batches.push_back(EdgeBatch{timestamp, {}});
        }
        batches.back().edgeIds.push_back(static_cast<EdgeId>(index));
    }
    return batches;
}

}
