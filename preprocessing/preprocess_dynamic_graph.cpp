#include <algorithm>
#include <charconv>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr std::size_t ioBufferSize = 1U << 20U;
constexpr std::string_view outputPrefix = "processed_";

struct Options {
    fs::path inputPath;
    fs::path outputPath;
    bool force = false;
};

struct Record {
    int source = 0;
    int destination = 0;
    std::int64_t timestamp = 0;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error(message);
}

void printUsage(const char* program) {
    std::cerr
        << "usage: " << program << " INPUT [--force]\n"
        << "\n"
        << "Extract the first two columns and the last column from a temporal\n"
        << "edge list. Vertex IDs must be signed 32-bit integers; timestamps\n"
        << "must be signed 64-bit integers. Middle columns are ignored.\n"
        << "Timestamp 0 means unknown and is skipped. Distinct timestamps\n"
        << "are mapped in ascending order\n"
        << "to consecutive IDs 0, 1, 2, ... . Lines beginning with % or # are\n"
        << "ignored. The output is written beside INPUT with the prefix \""
        << outputPrefix << "\". Existing output files are not overwritten unless\n"
        << "--force is specified.\n";
}

Options parseOptions(const int argc, char* argv[]) {
    if (argc == 2 && (std::string_view(argv[1]) == "--help" ||
                      std::string_view(argv[1]) == "-h")) {
        printUsage(argv[0]);
        std::exit(0);
    }
    if (argc != 2 && argc != 3) {
        printUsage(argv[0]);
        fail("expected one INPUT path");
    }

    const fs::path inputPath(argv[1]);
    const fs::path outputPath = inputPath.parent_path() /
        (std::string(outputPrefix) + inputPath.filename().string());
    Options options{inputPath, outputPath, false};
    if (argc == 3) {
        if (std::string_view(argv[2]) != "--force") {
            printUsage(argv[0]);
            fail("unknown option: " + std::string(argv[2]));
        }
        options.force = true;
    }
    return options;
}

void validatePaths(const Options& options) {
    if (!fs::is_regular_file(options.inputPath)) {
        fail("input is not a regular file: " + options.inputPath.string());
    }

    if (fs::exists(options.outputPath)) {
        if (fs::equivalent(options.inputPath, options.outputPath)) {
            fail("INPUT and OUTPUT must be different files");
        }
        if (!options.force) {
            fail("output already exists (use --force to overwrite): " +
                 options.outputPath.string());
        }
    }
}

std::string_view nextToken(std::string_view& line) {
    constexpr std::string_view whitespace = " \t\r\n\f\v";
    const auto begin = line.find_first_not_of(whitespace);
    if (begin == std::string_view::npos) {
        line = {};
        return {};
    }
    line.remove_prefix(begin);
    const auto token = line.substr(0, line.find_first_of(whitespace));
    line.remove_prefix(token.size());
    return token;
}

template <std::integral Integer>
Integer parseInteger(
    std::string_view token,
    const char* field,
    const std::size_t lineNumber) {
    const auto originalToken = token;
    const bool hasPlus = token.starts_with('+');
    if (hasPlus) {
        token.remove_prefix(1);
    }
    Integer value = 0;
    const char* const end = token.data() + token.size();
    const auto [parsedEnd, error] = std::from_chars(token.data(), end, value);
    if (error != std::errc{} || parsedEnd != end ||
        (hasPlus && token.starts_with('-'))) {
        fail("line " + std::to_string(lineNumber) +
             " has an invalid or out-of-range " + field + ": " +
             std::string(originalToken));
    }
    return value;
}

std::optional<Record> parseRecord(
    std::string_view line,
    const std::size_t lineNumber) {
    const auto source = nextToken(line);
    if (source.empty() || source.front() == '%' || source.front() == '#') {
        return std::nullopt;
    }
    const auto destination = nextToken(line);
    auto timestamp = nextToken(line);
    if (destination.empty() || timestamp.empty()) {
        fail("line " + std::to_string(lineNumber) +
             " must contain at least three columns: u v ... timestamp");
    }
    for (auto token = nextToken(line); !token.empty(); token = nextToken(line)) {
        timestamp = token;
    }

    const Record record{
        parseInteger<int>(source, "source vertex ID", lineNumber),
        parseInteger<int>(destination, "destination vertex ID", lineNumber),
        parseInteger<std::int64_t>(timestamp, "timestamp", lineNumber)};
    if (record.timestamp == 0) {
        return std::nullopt;
    }
    return record;
}

template <typename Callback>
std::size_t forEachRecord(const fs::path& inputPath, Callback&& callback) {
    std::vector<char> buffer(ioBufferSize);
    std::filebuf fileBuffer;
    fileBuffer.pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    if (fileBuffer.open(inputPath, std::ios::in | std::ios::binary) == nullptr) {
        fail("cannot open input: " + inputPath.string());
    }

    std::istream input(&fileBuffer);
    std::string line;
    std::size_t lineNumber = 0;
    std::size_t recordCount = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        const std::optional<Record> record = parseRecord(line, lineNumber);
        if (record.has_value()) {
            callback(*record);
            ++recordCount;
        }
    }
    if (!input.eof()) {
        fail("I/O error while reading input: " + inputPath.string());
    }
    return recordCount;
}

std::vector<std::int64_t> collectSortedTimestamps(
    const fs::path& inputPath,
    std::size_t& recordCount,
    bool& inputWasTimeSorted) {
    std::vector<std::int64_t> timestamps;
    std::unordered_set<std::int64_t> seenTimestamps;
    std::optional<std::int64_t> previousTimestamp;
    inputWasTimeSorted = true;

    recordCount = forEachRecord(inputPath, [&](const Record& record) {
        if (previousTimestamp.has_value() && record.timestamp < *previousTimestamp) {
            inputWasTimeSorted = false;
        }

        if (seenTimestamps.insert(record.timestamp).second) {
            if (timestamps.size() >
                static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                fail("too many distinct timestamps for 32-bit timestamp IDs");
            }
            timestamps.push_back(record.timestamp);
        }
        previousTimestamp = record.timestamp;
    });

    std::ranges::sort(timestamps);
    return timestamps;
}

void writeOutput(
    const Options& options,
    const std::vector<std::int64_t>& timestamps,
    const std::size_t expectedRecordCount) {
    std::vector<char> buffer(ioBufferSize);
    std::filebuf fileBuffer;
    fileBuffer.pubsetbuf(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    if (fileBuffer.open(
            options.outputPath,
            std::ios::out | std::ios::binary | std::ios::trunc) == nullptr) {
        fail("cannot open output: " + options.outputPath.string());
    }

    std::ostream output(&fileBuffer);
    const std::size_t writtenRecordCount = forEachRecord(
        options.inputPath,
        [&](const Record& record) {
            const auto timestampPosition =
                std::lower_bound(timestamps.begin(), timestamps.end(), record.timestamp);
            if (timestampPosition == timestamps.end() ||
                *timestampPosition != record.timestamp) {
                fail("internal error: timestamp mapping is incomplete");
            }
            const auto timestampId = static_cast<int>(timestampPosition - timestamps.begin());
            output << record.source << ' ' << record.destination << ' '
                   << timestampId << '\n';
            if (!output) {
                fail("I/O error while writing output: " +
                     options.outputPath.string());
            }
        });

    if (writtenRecordCount != expectedRecordCount) {
        fail("input changed while it was being processed");
    }
    output.flush();
    if (!output) {
        fail("I/O error while finalizing output: " + options.outputPath.string());
    }
}

}

int main(int argc, char* argv[]) {
    try {
        const Options options = parseOptions(argc, argv);
        validatePaths(options);

        std::size_t recordCount = 0;
        bool inputWasTimeSorted = false;
        const std::vector<std::int64_t> timestamps = collectSortedTimestamps(
            options.inputPath, recordCount, inputWasTimeSorted);
        writeOutput(options, timestamps, recordCount);

        std::cerr << "processed: output=" << options.outputPath.string()
                  << ", edges=" << recordCount
                  << ", distinct-timestamps=" << timestamps.size()
                  << ", timestamp-ids=";
        if (timestamps.empty()) {
            std::cerr << "empty";
        } else {
            std::cerr << "0.." << (timestamps.size() - 1U);
        }
        std::cerr << ", input-time-sorted="
                  << (inputWasTimeSorted ? "yes" : "no") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
