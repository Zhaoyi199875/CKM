#include "insertion_stream.hpp"
#include "maintenance_engine.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <time.h>
#include <unistd.h>

static timespec measuredBeginTime{};

extern "C" void stopAtMeasuredTimeLimit(int) {
    timespec now{};
    ::clock_gettime(CLOCK_MONOTONIC, &now);
    long long milliseconds =
        (static_cast<long long>(now.tv_sec) - measuredBeginTime.tv_sec) * 1000LL +
        (static_cast<long long>(now.tv_nsec) - measuredBeginTime.tv_nsec) / 1000000LL;
    if (milliseconds < 0) {
        milliseconds = 0;
    }
    static constexpr char head[] = "progress: measured-time-limit, measured-wall-ms=";
    char line[sizeof head + 24U];
    std::size_t length = 0U;
    for (std::size_t index = 0U; index + 1U < sizeof head; ++index) {
        line[length++] = head[index];
    }
    char digits[24];
    std::size_t count = 0U;
    do {
        digits[count++] = static_cast<char>('0' + milliseconds % 10);
        milliseconds /= 10;
    } while (milliseconds > 0);
    while (count > 0U) {
        line[length++] = digits[--count];
    }
    line[length++] = '\n';
    const ssize_t written = ::write(STDERR_FILENO, line, length);
    static_cast<void>(written);
    ::_exit(124);
}

namespace ckm {
namespace {

constexpr std::size_t defaultBatchSize = 1000U;

[[nodiscard]] double ratioOrNaN(double numerator, double denominator) noexcept {
    return denominator == 0.0
        ? std::numeric_limits<double>::quiet_NaN()
        : numerator / denominator;
}

enum class UpdateModeOption {
    Insertion,
    Deletion,
    Hybrid,
};

enum class ChangeOutputOption {
    New,
    Removed,
    Both,
};

struct Options {
    std::string inputPath;
    std::optional<std::size_t> batchSize = defaultBatchSize;
    std::optional<std::size_t> maxBatches;
    std::optional<std::size_t> maxEdges;
    UpdateModeOption updateMode = UpdateModeOption::Insertion;
    bool printProgress = false;
    std::optional<std::size_t> measuredTimeLimitSeconds;
    ChangeOutputOption changeOutput = ChangeOutputOption::Both;
    std::optional<std::string> outputPath;
    bool suppressOutput = false;
    bool printStats = false;
};

[[nodiscard]] constexpr std::string_view usageText() noexcept {
    return
        "usage: ckm INPUT [--batch-size=N|--batch-by-timestamp] "
        "[--max-batches=N] [--max-edges=N] "
        "[--update-mode=insertion|deletion|hybrid] "
        "[--change-output=new|removed|both] "
        "[--output=FILE|--no-output] [--stats] [--progress] "
        "[--measured-time-limit=SECONDS] "
        "(default: --update-mode=insertion --batch-size=1000 "
        "--change-output=both)";
}

[[noreturn]] void usageError(const std::string& message) {
    throw std::invalid_argument(message + "\n" + std::string(usageText()));
}

[[nodiscard]] std::size_t parsePositiveSize(
    std::string_view text,
    std::string_view optionName) {
    std::size_t value = 0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() || value == 0) {
        usageError(std::string(optionName) + " must be a positive integer");
    }
    return value;
}

[[nodiscard]] const char* updateModeOptionName(
    UpdateModeOption updateMode) noexcept {
    switch (updateMode) {
    case UpdateModeOption::Deletion:
        return "deletion";
    case UpdateModeOption::Hybrid:
        return "hybrid";
    case UpdateModeOption::Insertion:
        break;
    }
    return "insertion";
}

[[nodiscard]] const char* changeOutputOptionName(
    ChangeOutputOption changeOutput) noexcept {
    switch (changeOutput) {
    case ChangeOutputOption::New:
        return "new";
    case ChangeOutputOption::Removed:
        return "removed";
    case ChangeOutputOption::Both:
        break;
    }
    return "both";
}

[[nodiscard]] Options parseOptions(int argc, char** argv) {
    Options options;
    bool batchModeSeen = false;
    for (int index = 1; index < argc; ++index) {
        const std::string_view argument(argv[index]);
        if (argument.starts_with("--batch-size=")) {
            if (batchModeSeen) {
                usageError("batch mode may be specified only once");
            }
            batchModeSeen = true;
            options.batchSize =
                parsePositiveSize(argument.substr(13), "--batch-size");
        } else if (argument == "--batch-by-timestamp") {
            if (batchModeSeen) {
                usageError("batch mode may be specified only once");
            }
            batchModeSeen = true;
            options.batchSize.reset();
        } else if (argument.starts_with("--max-batches=")) {
            options.maxBatches =
                parsePositiveSize(argument.substr(14), "--max-batches");
        } else if (argument.starts_with("--max-edges=")) {
            options.maxEdges =
                parsePositiveSize(argument.substr(12), "--max-edges");
        } else if (argument.starts_with("--update-mode=")) {
            const std::string_view value = argument.substr(14);
            if (value == "insertion") {
                options.updateMode = UpdateModeOption::Insertion;
            } else if (value == "deletion") {
                options.updateMode = UpdateModeOption::Deletion;
            } else if (value == "hybrid") {
                options.updateMode = UpdateModeOption::Hybrid;
            } else {
                usageError("unknown update mode");
            }
        } else if (argument.starts_with("--change-output=")) {
            const auto value = argument.substr(16);
            if (value == "new") {
                options.changeOutput = ChangeOutputOption::New;
            } else if (value == "removed") {
                options.changeOutput = ChangeOutputOption::Removed;
            } else if (value == "both") {
                options.changeOutput = ChangeOutputOption::Both;
            } else {
                usageError("unknown change output mode");
            }
        } else if (argument.starts_with("--output=")) {
            options.outputPath = std::string(argument.substr(9));
        } else if (argument == "--no-output") {
            options.suppressOutput = true;
        } else if (argument == "--stats") {
            options.printStats = true;
        } else if (argument == "--progress") {
            options.printProgress = true;
        } else if (argument.starts_with("--measured-time-limit=")) {
            const std::size_t seconds = parsePositiveSize(
                argument.substr(22), "--measured-time-limit");
            if (seconds > std::numeric_limits<unsigned>::max()) {
                usageError("--measured-time-limit is too large");
            }
            options.measuredTimeLimitSeconds = seconds;
        } else if (argument == "--help" || argument == "-h") {
            std::cout << usageText() << '\n';
            std::exit(0);
        } else if (argument.starts_with('-')) {
            usageError("unknown option: " + std::string(argument));
        } else if (options.inputPath.empty()) {
            options.inputPath = std::string(argument);
        } else {
            usageError("only one input file is allowed");
        }
    }
    if (options.inputPath.empty()) {
        usageError("missing input file");
    }
    if (options.suppressOutput && options.outputPath.has_value()) {
        usageError("--output and --no-output cannot be combined");
    }
    if (options.maxEdges.has_value() && !options.batchSize.has_value()) {
        usageError("--max-edges needs --batch-size: the stream is cut by edge "
                   "count");
    }
    if (options.updateMode == UpdateModeOption::Hybrid &&
        !options.batchSize.has_value()) {
        usageError("--update-mode=hybrid needs --batch-size: the window "
                   "slides by edge count");
    }
    return options;
}

template <typename Callable>
[[nodiscard]] auto measurePhase(
    Callable&& operation,
    double& elapsedMilliseconds) {
    const auto begin = std::chrono::steady_clock::now();
    auto result = operation();
    elapsedMilliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - begin).count();
    return result;
}

void printLocalGraphStatistics(const LocalGraphStats& stats) {
    std::cerr
        << ", local-builds=" << stats.builds
        << ", local-vertices=" << stats.vertexOccurrences
        << ", local-edges=" << stats.edgeOccurrences
        << ", max-local-vertices=" << stats.maximumVertices
        << ", peak-local-bytes=" << stats.peakBytes;
}

void printCoverSeedStatistics(
    const CoverSeedStats& stats,
    double elapsedMilliseconds) {
    std::cerr
        << "coverseed: time-ms=" << elapsedMilliseconds
        << ", batches=" << stats.batches
        << ", batch-edges=" << stats.batchEdges
        << ", cover-vertices=" << stats.coverVertices
        << ", cover-both-endpoints=" << stats.coverBothEndpoints
        << ", subproblems=" << stats.seeds
        << ", single-edge-seeds=" << stats.singleEdgeSeeds
        << ", anchor-root-seeds=" << stats.anchorRootSeeds
        << ", pivot-root-seeds=" << stats.pivotRootSeeds
        << ", edges-per-seed="
        << ratioOrNaN(
               static_cast<double>(stats.assignedEdges),
               static_cast<double>(stats.seeds))
        << ", max-seed-edges=" << stats.maximumSeedEdges
        << ", root-anchor-branches=" << stats.rootAnchorBranches
        << ", root-pivot-branches=" << stats.rootPivotBranches
        << ", root-chosen-branches=" << stats.rootChosenBranches
        << ", root-degree-decisions=" << stats.rootDegreeDecisions
        << ", root-pivot-scans=" << stats.rootPivotScans
        << ", unanchored-anchor-choices=" << stats.unanchoredAnchorChoices
        << ", single-partner-anchors=" << stats.singlePartnerAnchors
        << ", domain-vertices=" << stats.domainVertexOccurrences
        << ", domain-vertices-per-seed="
        << ratioOrNaN(
               static_cast<double>(stats.domainVertexOccurrences),
               static_cast<double>(stats.seeds))
        << ", max-domain-vertices=" << stats.maximumDomainVertices
        << ", assign-ms=" << stats.assignMilliseconds
        << ", recursive-calls=" << stats.recursiveCalls
        << ", pivot-calls=" << stats.pivotCalls
        << ", pivot-candidates=" << stats.pivotCandidatesEvaluated
        << ", leaves=" << stats.leaves
        << ", unanchored-nodes=" << stats.unanchoredNodes
        << ", live-test-prunes=" << stats.liveTestPrunes
        << ", emitted=" << stats.emittedCliques
        << ", ownership-moves=" << stats.ownershipMoves
        << ", ownership-probes=" << stats.ownershipProbes
        << ", reindexed-trees=" << stats.reindexedTrees
        << ", root-child-reindexed-subtrees=" << stats.rootChildReindexedSubtrees
        << ", deep-reindexed-subtrees=" << stats.deepReindexedSubtrees
        << ", node-words=" << stats.nodeWords
        << ", pivot-words=" << stats.pivotWords;
    printLocalGraphStatistics(stats.localGraph);
    std::cerr << '\n';
}

void printSubsumedStatistics(
    const SubsumedSearchStats& stats,
    double elapsedMilliseconds) {
    std::cerr
        << "subsumed: time-ms=" << elapsedMilliseconds
        << ", batches=" << stats.batches
        << ", new-cliques-processed=" << stats.newCliquesProcessed
        << ", internal-batch-edges=" << stats.internalBatchEdges
        << ", one-edge-kernels=" << stats.oneEdgeKernels
        << ", index-lookups=" << stats.indexLookups
        << ", size-gate-skipped-cliques=" << stats.sizeGateSkippedCliques
        << ", size-gate-skipped-probes=" << stats.sizeGateSkippedProbes
        << ", kernel-candidates=" << stats.kernelCandidates
        << ", kernel-searches=" << stats.kernelSearches
        << ", kernel-search-nodes=" << stats.kernelSearchNodes
        << ", emitted=" << stats.emittedCliques
        << '\n';
}

}
}

int main(int argc, char** argv) {
    using namespace ckm;
    try {
        const auto programBegin = std::chrono::steady_clock::now();
        const Options options = parseOptions(argc, argv);
        double inputPreprocessMilliseconds = 0.0;
        double batchBuildMilliseconds = 0.0;
        double graphInitializationMilliseconds = 0.0;
        InsertionStream stream = measurePhase(
            [&] {
                return InsertionStream::readEdgeList(options.inputPath);
            },
            inputPreprocessMilliseconds);
        const std::size_t streamEdgeCount = stream.insertionTimestamps().size();
        if (options.maxEdges.value_or(0U) > streamEdgeCount) {
            throw std::invalid_argument(
                "--max-edges is larger than the stream's " +
                std::to_string(streamEdgeCount) + " edges");
        }
        const std::size_t runEdgeCount =
            options.maxEdges.value_or(streamEdgeCount);
        const bool deletion = options.updateMode == UpdateModeOption::Deletion;
        const bool hybrid = options.updateMode == UpdateModeOption::Hybrid;
        const std::size_t windowEdges = hybrid ? runEdgeCount / 2U : 0U;
        if (hybrid && windowEdges < *options.batchSize) {
            throw std::invalid_argument(
                "the " + std::to_string(windowEdges) +
                "-edge window is smaller than a batch, which would delete "
                "edges it has not inserted yet");
        }
        std::size_t windowBatchCount = 0U;
        std::vector<EdgeBatch> batches = measurePhase(
            [&] {
                std::vector<EdgeBatch> all;
                if (windowEdges != 0U) {
                    all = stream.makeBatches(
                        options.batchSize, std::nullopt, 0U, windowEdges);
                }
                windowBatchCount = all.size();
                std::vector<EdgeBatch> measured = stream.makeBatches(
                    options.batchSize, options.maxBatches, windowEdges,
                    options.maxEdges.value_or(
                        std::numeric_limits<std::size_t>::max()));
                all.insert(all.end(), std::make_move_iterator(measured.begin()),
                           std::make_move_iterator(measured.end()));
                return all;
            },
            batchBuildMilliseconds);
        const std::size_t availableBatchCount = options.batchSize.has_value()
            ? (runEdgeCount - windowEdges + *options.batchSize - 1U) /
                  *options.batchSize
            : stream.batchCount(options.batchSize);
        std::size_t processedEdgeCount = 0;
        for (std::size_t index = windowBatchCount; index < batches.size();
             ++index) {
            processedEdgeCount += batches[index].edgeIds.size();
        }
        IncrementalGraph graph = measurePhase(
            [&] {
                return IncrementalGraph(
                    stream.vertexCount(), stream.takeEdges());
            },
            graphInitializationMilliseconds);

        MaintenanceConfiguration configuration;
        configuration.builds = deletion || hybrid;
        configuration.sizeHistogram = !deletion;
        MaintenanceEngine engine(graph, configuration);

        std::ofstream outputFile;
        std::ostream* output = &std::cout;
        if (options.outputPath.has_value()) {
            outputFile.open(*options.outputPath);
            if (!outputFile) {
                throw std::runtime_error(
                    "cannot open output file: " + *options.outputPath);
            }
            output = &outputFile;
        }

        enum class StepKind { Build, Insert, Delete, Slide };
        struct Step {
            StepKind kind = StepKind::Insert;
            std::size_t batch = 0;
        };
        std::vector<Step> steps;
        steps.reserve(batches.size() * (deletion ? 2U : 1U));
        for (std::size_t index = 0; index < batches.size(); ++index) {
            StepKind kind = deletion ? StepKind::Build : StepKind::Insert;
            if (hybrid) {
                kind = index < windowBatchCount ? StepKind::Build
                                                : StepKind::Slide;
            }
            steps.push_back(Step{kind, index});
        }
        if (deletion) {
            for (std::size_t index = batches.size(); index > 0U; --index) {
                steps.push_back(Step{StepKind::Delete, index - 1U});
            }
        }
        const std::size_t measuredBegin =
            deletion ? batches.size() : windowBatchCount;
        const std::size_t buildEdges = deletion ? processedEdgeCount : windowEdges;

        BatchChanges changes;
        const auto writeChanges = [&](int timestamp) {
            if (options.suppressOutput) {
                return;
            }
            std::ranges::sort(changes.created, CliqueOutputOrder{});
            std::ranges::sort(changes.removed, CliqueOutputOrder{});
            const auto writeClique = [&](std::string_view change,
                                         const Clique& clique) {
                *output << timestamp << '\t' << change << '\t';
                for (std::size_t index = 0; index < clique.size(); ++index) {
                    if (index != 0) {
                        *output << ' ';
                    }
                    *output << stream.originalVertexId(clique[index]);
                }
                *output << '\n';
            };
            if (options.changeOutput == ChangeOutputOption::New ||
                options.changeOutput == ChangeOutputOption::Both) {
                for (const Clique& clique : changes.created) {
                    writeClique("new", clique);
                }
            }
            if (options.changeOutput == ChangeOutputOption::Removed ||
                options.changeOutput == ChangeOutputOption::Both) {
                for (const Clique& clique : changes.removed) {
                    writeClique("del", clique);
                }
            }
        };
        std::vector<EdgeId> expired;
        const auto searchAndMaintenanceMilliseconds = [&] {
            const MaintenanceTimes& times = engine.times();
            return times.coverSeed + times.subsumed + times.cliqueIndexUpdate +
                times.transientNetting + times.singleVertex;
        };
        const auto wallMilliseconds = [&] {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - programBegin).count();
        };
        double buildMilliseconds = 0.0;
        const auto buildBegin = std::chrono::steady_clock::now();

        for (std::size_t stepIndex = 0; stepIndex < steps.size(); ++stepIndex) {
            if (stepIndex == measuredBegin) {
                if (measuredBegin != 0U) {
                    buildMilliseconds = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - buildBegin).count();
                    engine.resetMeasurements();
                }
                if (options.printProgress) {
                    std::cerr << "progress: measured-begin, measured-batches=0/"
                              << steps.size() - measuredBegin
                              << ", wall-ms=" << wallMilliseconds() << std::endl;
                }
                if (options.measuredTimeLimitSeconds.has_value()) {
                    ::clock_gettime(CLOCK_MONOTONIC, &measuredBeginTime);
                    std::signal(SIGALRM, stopAtMeasuredTimeLimit);
                    ::alarm(static_cast<unsigned>(
                        *options.measuredTimeLimitSeconds));
                }
            }
            const Step& step = steps[stepIndex];
            const EdgeBatch& batch = batches[step.batch];
            switch (step.kind) {
            case StepKind::Build:
                if (!batch.edgeIds.empty()) {
                    engine.build(batch.edgeIds);
                }
                break;
            case StepKind::Delete:
                engine.updateBatch({}, batch.edgeIds, changes);
                writeChanges(batch.timestamp);
                break;
            case StepKind::Slide: {
                const auto first = static_cast<std::size_t>(batch.edgeIds.front());
                expired.clear();
                for (std::size_t offset = 0; offset < batch.edgeIds.size();
                     ++offset) {
                    expired.push_back(
                        static_cast<EdgeId>(first - windowEdges + offset));
                }
                engine.updateBatch(batch.edgeIds, expired, changes);
                writeChanges(batch.timestamp);
                break;
            }
            case StepKind::Insert:
                engine.updateBatch(batch.edgeIds, {}, changes);
                writeChanges(batch.timestamp);
                break;
            }
            if (options.printProgress && stepIndex >= measuredBegin) {
                std::cerr << "progress: measured-batches="
                          << stepIndex + 1U - measuredBegin << '/'
                          << steps.size() - measuredBegin
                          << ", search-and-maintenance-ms="
                          << searchAndMaintenanceMilliseconds()
                          << ", wall-ms=" << wallMilliseconds() << std::endl;
            }
        }
        if (options.measuredTimeLimitSeconds.has_value()) {
            ::alarm(0U);
        }
        if (measuredBegin != 0U && measuredBegin == steps.size()) {
            buildMilliseconds = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - buildBegin).count();
            engine.resetMeasurements();
        }

        if (options.printStats) {
            const MaintenanceTimes& times = engine.times();
            std::cerr
                << "graph: vertices=" << graph.vertexCount()
                << ", unique-edges=" << graph.edgeCount()
                << ", max-edges=";
            if (options.maxEdges.has_value()) {
                std::cerr << *options.maxEdges;
            } else {
                std::cerr << "all";
            }
            std::cerr
                << ", processed-edges=" << processedEdgeCount
                << ", batches=" << batches.size() - windowBatchCount
                << ", available-batches=" << availableBatchCount
                << ", max-batches=";
            if (options.maxBatches.has_value()) {
                std::cerr << *options.maxBatches;
            } else {
                std::cerr << "all";
            }
            std::cerr
                << ", build-edges=" << buildEdges
                << ", build-ms=" << buildMilliseconds;
            const UpdateCounts& counts = engine.updateCounts();
            std::cerr
                << ", new-maximal-cliques=" << engine.totalCreated()
                << ", removed-maximal-cliques=" << engine.totalRemoved()
                << ", current-maximal-cliques=" << engine.currentCliqueCount()
                << ", new-single-vertex-cliques=" << counts.singleVertexCreated
                << ", removed-single-vertex-cliques="
                << counts.singleVertexRemoved
                << ", isolated-vertices=" << engine.isolatedVertexCount()
                << ", batch-size=";
            if (options.batchSize.has_value()) {
                std::cerr << *options.batchSize;
            } else {
                std::cerr << "timestamp";
            }
            std::cerr
                << ", update-mode=" << updateModeOptionName(options.updateMode)
                << ", change-output="
                << changeOutputOptionName(options.changeOutput)
                << ", measured-time-limit=";
            if (options.measuredTimeLimitSeconds.has_value()) {
                std::cerr << *options.measuredTimeLimitSeconds;
            } else {
                std::cerr << "none";
            }
            std::cerr << '\n';
            const double excludedMilliseconds = inputPreprocessMilliseconds +
                batchBuildMilliseconds + graphInitializationMilliseconds +
                times.graphUpdate;
            const double maintenanceMilliseconds = times.subsumed +
                times.cliqueIndexUpdate + times.transientNetting +
                times.singleVertex;
            std::cerr
                << "timing: input-preprocess-ms=" << inputPreprocessMilliseconds
                << ", batch-build-ms=" << batchBuildMilliseconds
                << ", graph-init-ms=" << graphInitializationMilliseconds
                << ", graph-update-ms=" << times.graphUpdate
                << ", clique-index-update-ms=" << times.cliqueIndexUpdate
                << ", full-maintenance-ms=" << maintenanceMilliseconds
                << ", search-and-maintenance-ms="
                << searchAndMaintenanceMilliseconds()
                << ", single-vertex-ms=" << times.singleVertex
                << ", excluded-ms=" << excludedMilliseconds
                << '\n';
            printCoverSeedStatistics(engine.coverSeedStats(), times.coverSeed);
            const MaximalCliqueIndex& maximalCliques = engine.index();
            std::cerr << "index: live=" << maximalCliques.size()
                      << ", size-histogram="
                      << (maximalCliques.tracksSizes() ? "on" : "off");
            if (maximalCliques.tracksSizes()) {
                std::cerr << ", min-live-size="
                          << maximalCliques.minimumLiveCliqueSize()
                          << ", sizes=";
                for (std::size_t k = 0;
                     k < MaximalCliqueIndex::trackedSizeLimit; ++k) {
                    const std::size_t count =
                        maximalCliques.liveCliqueCountOfSize(k);
                    if (count != 0U) {
                        std::cerr << k << ":" << count << " ";
                    }
                }
                if (maximalCliques.oversizeLiveCliqueCount() != 0U) {
                    std::cerr << "oversize:"
                              << maximalCliques.oversizeLiveCliqueCount()
                              << " ";
                }
            }
            std::cerr << '\n';
            printSubsumedStatistics(engine.subsumedStats(), times.subsumed);
            if (deletion || hybrid) {
                const SubsumedSearchStats::Recovery& recovery =
                    engine.subsumedStats().recovery;
                std::cerr << "recovery: destroyed=" << counts.destroyed
                          << ", created=" << counts.recovered
                          << ", candidates=" << recovery.candidates
                          << ", below-two=" << recovery.belowTwo
                          << ", duplicates=" << recovery.duplicates
                          << ", checks=" << recovery.checks
                          << ", rejected=" << recovery.rejected
                          << ", accepted=" << recovery.accepted
                          << '\n';
            }
            if (hybrid) {
                std::cerr << "hybrid: window-edges=" << windowEdges
                          << ", batches=" << counts.batches
                          << ", inserted-edges=" << counts.insertedEdges
                          << ", deleted-edges=" << counts.deletedEdges
                          << ", destroyed=" << counts.destroyed
                          << ", recovered=" << counts.recovered
                          << ", new=" << counts.holdingNew
                          << ", subsumed=" << counts.subsumed
                          << ", transient=" << counts.transient
                          << ", transient-single-vertex="
                          << counts.transientSingleVertex
                          << ", destroyed-search-ms=" << times.destroyedSearch
                          << ", new-search-ms=" << times.newSearch
                          << ", subsumed-ms=" << times.subsumed - times.recovery
                          << ", recovery-ms=" << times.recovery
                          << ", transient-netting-ms=" << times.transientNetting
                          << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
