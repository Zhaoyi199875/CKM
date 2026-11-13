#include <algorithm>
#include <charconv>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
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
    std::int64_t source = 0;
    std::int64_t destination = 0;
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
        << "edge list. Distinct UNIX timestamps are mapped in ascending order\n"
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

std::string pathForMessage(const fs::path& path) {
    return path.string();
}

void validatePaths(const Options& options) {
    std::error_code error;
    if (!fs::is_regular_file(options.inputPath, error)) {
        fail("input is not a regular file: " + pathForMessage(options.inputPath));
    }

    const fs::path inputAbsolute = fs::absolute(options.inputPath).lexically_normal();
    const fs::path outputAbsolute = fs::absolute(options.outputPath).lexically_normal();
    if (inputAbsolute == outputAbsolute) {
        fail("INPUT and OUTPUT must be different files");
    }

    const bool outputExists = fs::exists(options.outputPath, error);
    if (outputExists) {
        if (fs::equivalent(options.inputPath, options.outputPath, error)) {
            fail("INPUT and OUTPUT must be different files");
        }
        if (!options.force) {
            fail("output already exists (use --force to overwrite): " +
                 pathForMessage(options.outputPath));
        }
    }

    const fs::path parent = options.outputPath.parent_path();
    if (!parent.empty() && !fs::is_directory(parent, error)) {
        fail("output directory does not exist: " + pathForMessage(parent));
    }
}

std::optional<Record> parseRecord(
    const std::string_view line,
    const std::size_t lineNumber) {
    std::size_t position = 0;
    while (position < line.size() &&
           std::isspace(static_cast<unsigned char>(line[position])) != 0) {
        ++position;
    }
    if (position == line.size() || line[position] == '%' || line[position] == '#') {
        return std::nullopt;
    }

    Record record;
    std::size_t column = 0;
    while (position < line.size()) {
        while (position < line.size() &&
               std::isspace(static_cast<unsigned char>(line[position])) != 0) {
            ++position;
        }
        if (position == line.size()) {
            break;
        }

        const std::size_t tokenBegin = position;
        while (position < line.size() &&
               std::isspace(static_cast<unsigned char>(line[position])) == 0) {
            ++position;
        }

        std::int64_t value = 0;
        const char* const begin = line.data() + tokenBegin;
        const char* const end = line.data() + position;
        const auto [parsedEnd, parseError] = std::from_chars(begin, end, value);
        if (parseError != std::errc{} || parsedEnd != end) {
            fail("line " + std::to_string(lineNumber) +
                 " contains a non-integer column: " +
                 std::string(begin, end));
        }

        if (column == 0U) {
            record.source = value;
        } else if (column == 1U) {
            record.destination = value;
        }
        record.timestamp = value;
        ++column;
    }

    if (column < 3U) {
        fail("line " + std::to_string(lineNumber) +
             " must contain at least three integer columns");
    }
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
        fail("cannot open input: " + pathForMessage(inputPath));
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
        fail("I/O error while reading input: " + pathForMessage(inputPath));
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
        fail("cannot open output: " + pathForMessage(options.outputPath));
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
            const auto timestampId = static_cast<std::size_t>(
                std::distance(timestamps.begin(), timestampPosition));
            output << record.source << ' ' << record.destination << ' '
                   << timestampId << '\n';
            if (!output) {
                fail("I/O error while writing output: " +
                     pathForMessage(options.outputPath));
            }
        });

    if (writtenRecordCount != expectedRecordCount) {
        fail("input changed while it was being processed");
    }
    output.flush();
    if (!output) {
        fail("I/O error while finalizing output: " + pathForMessage(options.outputPath));
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

        std::cerr << "processed: output=" << pathForMessage(options.outputPath)
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
