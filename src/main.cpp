#include "insertion_stream.hpp"
#include "maintenance_engine.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <time.h>
#include <unistd.h>

// --measured-time-limit and --measured-time-extend: where the measured
// batches began (CLOCK_MONOTONIC), set before the alarm is armed.
static timespec immceMeasuredBegin{};

// The alarm ends the run where it stands, as an external `timeout` would,
// with the same exit status (124), after a line that says how long the
// measured batches ran.  clock_gettime(2), write(2) and _exit(2) are
// async-signal-safe; the progress lines before it are flushed.  Called
// directly too, when the batches stop between two batches.
extern "C" void immceStopAtMeasuredTimeLimit(int /*signal*/) {
    timespec now{};
    ::clock_gettime(CLOCK_MONOTONIC, &now);
    long long milliseconds =
        (static_cast<long long>(now.tv_sec) - immceMeasuredBegin.tv_sec) * 1000LL +
        (static_cast<long long>(now.tv_nsec) - immceMeasuredBegin.tv_nsec) / 1000000LL;
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

namespace immce {
namespace {

constexpr std::size_t defaultBatchSize = 1000U;

[[nodiscard]] double ratioOrNaN(double numerator, double denominator) noexcept {
    return denominator == 0.0
        ? std::numeric_limits<double>::quiet_NaN()
        : numerator / denominator;
}

// The new-clique search: star (the method), edge or vertex (its baselines),
// or star and edge (both) or all three (all), compared batch by batch.  SOMEi
// (somei.hpp) is a baseline of its own: it replaces the new-clique search,
// the subsumed search and the clique index together, so it has no
// producer/consumer split.
enum class FrameworkOption {
    Star,
    Edge,
    Vertex,
    Both,
    All,
    Somei,
};

enum class MaintenanceOption {
    NewOnly,
    Full,
};

// Insert: the stream's batches go in, in order.  Reverse delete: every batch
// goes in first, untimed (the warm-up and the measured window alike), and the
// measured batches are then deleted again, last first -- the decremental
// workload of Das et al. (start from the graph, delete batches in decreasing
// timestamp order).  Random labels: their fully dynamic workload -- every
// edge of the measured window is labelled at random, insert or delete; the
// ones to delete go in first, untimed (with the warm-up), and the window's
// batches then insert their insert-labelled edges and delete their
// delete-labelled ones at once.  Sliding window: the graph keeps the newest W
// edges -- the first W go in untimed, and every batch then inserts its edges
// and deletes as many of the oldest, at once.
enum class UpdateModeOption {
    Insert,
    ReverseDelete,
    RandomLabels,
    SlidingWindow,
};

enum class ChangeOutputOption {
    New,
    Subsumed,
    Both,
};

enum class MaintenancePipelineOption {
    // Streaming wherever the configuration allows it, batched otherwise.
    Auto,
    Batched,
    Streaming,
};

struct Options {
    std::string inputPath;
    std::optional<std::size_t> batchSize = defaultBatchSize;
    std::optional<std::size_t> maxBatches;
    // Only the stream's first `maxEdges` edges are used (unset: all).
    std::optional<std::size_t> maxEdges;
    // Warm-up: the stream's first `warmupEdges` edges go in before anything is
    // measured, in batches of `warmupBatchSize` (the batch size when unset),
    // and then every timer and counter starts again from zero, so the
    // statistics cover only the batches after them; `--max-batches` counts
    // those alone.  The state the warm-up leaves -- the graph and its maximal
    // cliques -- does not depend on how it was batched.  Its changes are
    // still written, so a warm-up cut into the same batches as the run leaves
    // the change stream exactly as it is without one.
    std::size_t warmupEdges = 0U;
    std::optional<std::size_t> warmupBatchSize;
    // Off: in new-only mode the warm-up only inserts its edges. The searches
    // it would run leave no state behind there, only output, so a warm-up
    // into the middle of a large stream need not pay for them.
    bool warmupSearch = true;
    FrameworkOption framework = FrameworkOption::Star;
    // Which of the two star methods: the branch-count rule at the centre
    // only, or at every state that holds no partner yet.
    StarBranchScope starBranchScope = StarBranchScope::Unbound;
    // The rule, or for the branching ablation one of its two partitions fixed.
    StarBranching starBranching = StarBranching::Rule;
    // Off only for the ablation of the star's domain.
    bool starRootFilter = true;
    // How the vertex baseline branches before its clique holds a batch edge;
    // anything but the pivot is the branching ablation.
    VertexBranching vertexBranching = VertexBranching::Pivot;
    // Which batch endpoints seed the vertex baseline; only the cover for the
    // seed-cover ablation.
    VertexSeeds vertexSeeds = VertexSeeds::Endpoints;
    MaintenanceOption maintenance = MaintenanceOption::NewOnly;
    UpdateModeOption updateMode = UpdateModeOption::Insert;
    // Random labels: the chance that an edge is labelled insert, and the seed
    // of the 64-bit Mersenne Twister drawing the labels in stream order.
    std::optional<double> insertFraction;
    std::optional<std::uint64_t> labelSeed;
    // Sliding window: how many of the newest edges the graph keeps (unset:
    // half of the edges the run uses).
    std::optional<std::size_t> windowEdges;
    // Batches that both delete and insert (random labels, sliding window).
    std::optional<MixedStrategy> mixedStrategy;
    MaintenancePipelineOption maintenancePipeline =
        MaintenancePipelineOption::Auto;
    SubsumedCliqueAlgorithm subsumedAlgorithm =
        SubsumedCliqueAlgorithm::Optimized;
    // Whether the index keeps the live clique-size histogram, which is the
    // residual-kernel size gate: off switches the gate off.  Unset means "on
    // exactly when the algorithm reads it"; the controls never do, and on
    // only buys them the size distribution on the `index:` line.
    std::optional<bool> subsumedSizeHistogram;
    bool someiPruning = true;
    bool someiTrackChanges = true;
    bool someiProfilePhases = false;
    // SOMEi only: the untimed build (a deletion's whole graph, a sliding
    // window's first window, the delete-labelled edges) in one batch rather
    // than the run's own batches -- one construction instead of hundreds of
    // rollbacks, so that a time limit falls on the measured batches.
    bool someiOneBatchBuild = false;
    // SOMEi only: the untimed build's steps find the maximal cliques of G[R]
    // by Bron-Kerbosch (bk, the default) or by the paper's lines 7-22 of
    // Algorithm 6 (paper).  The same state either way (somei.hpp); the
    // measured batches always run the paper's steps.
    std::optional<bool> someiEnumerateBuild;
    // A line per measured batch on stderr, flushed: how far a run that hits
    // its time limit got.
    bool printProgress = false;
    // Wall seconds the measured batches may take, counted from where they
    // begin -- after the reading and any untimed build or warm-up.  When
    // they are up the run stops at once with exit status 124, as a
    // `timeout` would (2026-11-06: a time limit on a whole run let a long
    // build eat into SOMEi's measured time).
    std::optional<std::size_t> measuredTimeLimitSeconds;
    // With the limit: measured batches that look set to finish within this
    // many seconds -- the time so far plus the last batch's time for every
    // batch left, judged after each batch -- may run until then instead
    // (2026-11-06, so that a run a little over the limit still finishes).
    std::optional<std::size_t> measuredTimeExtendSeconds;
    ChangeOutputOption changeOutput = ChangeOutputOption::New;
    std::optional<std::string> outputPath;
    bool suppressOutput = false;
    bool printStats = false;
};

[[nodiscard]] constexpr std::string_view usageText() noexcept {
    return
        "usage: immce_star INPUT [--batch-size=N|--batch-by-timestamp] "
        "[--max-batches=N] [--max-edges=N] "
        "[--warmup-edges=N] [--warmup-batch-size=N] "
        "[--warmup-search=on|off] "
        "[--tie-order=file] "
        "[--framework=star|edge|vertex|both|all|somei] "
        "[--star-branch-scope=unbound|root] "
        "[--star-branching=rule|anchor|pivot] [--star-root-filter=on|off] "
        "[--vertex-branching=pivot|new-edges|rule-root|rule-unbound] "
        "[--vertex-seeds=endpoints|cover] "
        "[--maintenance=new-only|full] "
        "[--update-mode=insert|reverse-delete|random-labels|sliding-window] "
        "[--insert-fraction=P] [--label-seed=N] [--window-edges=N] "
        "[--mixed-strategy=two-phase|single-pass] "
        "[--maintenance-pipeline=auto|batched|streaming] "
        "[--subsumed-algorithm=optimized|optimized-bitmap|residual-kernel|full-bk|kernel-bk] "
        "[--subsumed-size-histogram=on|off] "
        "[--somei-pruning=on|off] [--somei-changes=on|off] "
        "[--somei-profile=on|off] [--somei-build=batches|one-batch] "
        "[--somei-build-update=bk|paper] "
        "[--change-output=new|subsumed|both] "
        "[--output=FILE|--no-output] [--stats] [--progress] "
        "[--measured-time-limit=SECONDS] [--measured-time-extend=SECONDS] "
        "(default: --framework=star --star-branch-scope=unbound "
        "--batch-size=1000 --tie-order=file)";
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

[[nodiscard]] std::size_t parseNonnegativeSize(
    std::string_view text,
    std::string_view optionName) {
    std::size_t value = 0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        usageError(std::string(optionName) + " must be a nonnegative integer");
    }
    return value;
}

[[nodiscard]] double parseFraction(
    std::string_view text,
    std::string_view optionName) {
    double value = 0.0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() ||
        !(value >= 0.0 && value <= 1.0)) {
        usageError(std::string(optionName) + " must be a number from 0 to 1");
    }
    return value;
}

[[nodiscard]] std::uint64_t parseUnsigned64(
    std::string_view text,
    std::string_view optionName) {
    std::uint64_t value = 0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        usageError(std::string(optionName) + " must be a nonnegative integer");
    }
    return value;
}

[[nodiscard]] bool parseOnOff(
    std::string_view text,
    std::string_view optionName) {
    if (text == "on") {
        return true;
    }
    if (text == "off") {
        return false;
    }
    usageError(std::string(optionName) + " must be on or off");
}

[[nodiscard]] const char* updateModeOptionName(
    UpdateModeOption updateMode) noexcept {
    switch (updateMode) {
    case UpdateModeOption::ReverseDelete:
        return "reverse-delete";
    case UpdateModeOption::RandomLabels:
        return "random-labels";
    case UpdateModeOption::SlidingWindow:
        return "sliding-window";
    case UpdateModeOption::Insert:
        break;
    }
    return "insert";
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
            if (options.maxBatches.has_value()) {
                usageError("--max-batches may be specified only once");
            }
            options.maxBatches =
                parsePositiveSize(argument.substr(14), "--max-batches");
        } else if (argument.starts_with("--max-edges=")) {
            options.maxEdges =
                parsePositiveSize(argument.substr(12), "--max-edges");
        } else if (argument.starts_with("--warmup-edges=")) {
            options.warmupEdges =
                parseNonnegativeSize(argument.substr(15), "--warmup-edges");
        } else if (argument.starts_with("--warmup-batch-size=")) {
            options.warmupBatchSize =
                parsePositiveSize(argument.substr(20), "--warmup-batch-size");
        } else if (argument.starts_with("--warmup-search=")) {
            options.warmupSearch =
                parseOnOff(argument.substr(16), "--warmup-search");
        } else if (argument.starts_with("--tie-order=")) {
            // Edges sharing a timestamp keep the input file's order, the only
            // order there is; the option stays so that scripts can state it.
            if (argument.substr(12) != "file") {
                usageError(
                    "--tie-order takes only file: edges sharing a timestamp "
                    "keep the input file's order (vertex and shuffle were "
                    "removed on 2026-10-08)");
            }
        } else if (argument.starts_with("--framework=")) {
            const auto value = argument.substr(12);
            if (value == "star") {
                options.framework = FrameworkOption::Star;
            } else if (value == "edge") {
                options.framework = FrameworkOption::Edge;
            } else if (value == "vertex") {
                options.framework = FrameworkOption::Vertex;
            } else if (value == "both") {
                options.framework = FrameworkOption::Both;
            } else if (value == "all") {
                options.framework = FrameworkOption::All;
            } else if (value == "somei") {
                options.framework = FrameworkOption::Somei;
            } else {
                usageError("unknown framework");
            }
        } else if (argument.starts_with("--star-branch-scope=")) {
            const auto value = argument.substr(20);
            if (value == "root") {
                options.starBranchScope = StarBranchScope::Root;
            } else if (value == "unbound") {
                options.starBranchScope = StarBranchScope::Unbound;
            } else {
                usageError("unknown star branch scope");
            }
        } else if (argument.starts_with("--star-branching=")) {
            const auto value = argument.substr(17);
            if (value == "rule") {
                options.starBranching = StarBranching::Rule;
            } else if (value == "anchor") {
                options.starBranching = StarBranching::Anchor;
            } else if (value == "pivot") {
                options.starBranching = StarBranching::Pivot;
            } else {
                usageError("unknown star branching");
            }
        } else if (argument.starts_with("--star-root-filter=")) {
            options.starRootFilter =
                parseOnOff(argument.substr(19), "--star-root-filter");
        } else if (argument.starts_with("--vertex-branching=")) {
            const auto value = argument.substr(19);
            if (value == "pivot") {
                options.vertexBranching = VertexBranching::Pivot;
            } else if (value == "new-edges") {
                options.vertexBranching = VertexBranching::NewEdges;
            } else if (value == "rule-root") {
                options.vertexBranching = VertexBranching::RuleRoot;
            } else if (value == "rule-unbound") {
                options.vertexBranching = VertexBranching::RuleUnbound;
            } else {
                usageError("unknown vertex branching");
            }
        } else if (argument.starts_with("--vertex-seeds=")) {
            const auto value = argument.substr(15);
            if (value == "endpoints") {
                options.vertexSeeds = VertexSeeds::Endpoints;
            } else if (value == "cover") {
                options.vertexSeeds = VertexSeeds::Cover;
            } else {
                usageError("unknown vertex seeds");
            }
        } else if (argument.starts_with("--maintenance=")) {
            const auto value = argument.substr(14);
            if (value == "new-only") {
                options.maintenance = MaintenanceOption::NewOnly;
            } else if (value == "full") {
                options.maintenance = MaintenanceOption::Full;
            } else {
                usageError("unknown maintenance mode");
            }
        } else if (argument.starts_with("--update-mode=")) {
            const std::string_view value = argument.substr(14);
            if (value == "insert") {
                options.updateMode = UpdateModeOption::Insert;
            } else if (value == "reverse-delete") {
                options.updateMode = UpdateModeOption::ReverseDelete;
            } else if (value == "random-labels") {
                options.updateMode = UpdateModeOption::RandomLabels;
            } else if (value == "sliding-window") {
                options.updateMode = UpdateModeOption::SlidingWindow;
            } else {
                usageError("unknown update mode");
            }
        } else if (argument.starts_with("--insert-fraction=")) {
            options.insertFraction =
                parseFraction(argument.substr(18), "--insert-fraction");
        } else if (argument.starts_with("--label-seed=")) {
            options.labelSeed =
                parseUnsigned64(argument.substr(13), "--label-seed");
        } else if (argument.starts_with("--window-edges=")) {
            options.windowEdges =
                parsePositiveSize(argument.substr(15), "--window-edges");
        } else if (argument.starts_with("--mixed-strategy=")) {
            const std::string_view value = argument.substr(17);
            if (value == "single-pass") {
                options.mixedStrategy = MixedStrategy::SinglePass;
            } else if (value == "two-phase") {
                options.mixedStrategy = MixedStrategy::TwoPhase;
            } else {
                usageError("unknown mixed strategy");
            }
        } else if (argument.starts_with("--maintenance-pipeline=")) {
            const std::string_view value = argument.substr(23);
            if (value == "auto") {
                options.maintenancePipeline = MaintenancePipelineOption::Auto;
            } else if (value == "batched") {
                options.maintenancePipeline =
                    MaintenancePipelineOption::Batched;
            } else if (value == "streaming") {
                options.maintenancePipeline =
                    MaintenancePipelineOption::Streaming;
            } else {
                usageError("unknown maintenance pipeline");
            }
        } else if (argument.starts_with("--subsumed-algorithm=")) {
            const auto value = argument.substr(21);
            if (value == "optimized") {
                options.subsumedAlgorithm = SubsumedCliqueAlgorithm::Optimized;
            } else if (value == "residual-kernel") {
                options.subsumedAlgorithm = SubsumedCliqueAlgorithm::ResidualKernel;
            } else if (value == "optimized-bitmap") {
                options.subsumedAlgorithm = SubsumedCliqueAlgorithm::OptimizedBitmap;
            } else if (value == "full-bk") {
                options.subsumedAlgorithm = SubsumedCliqueAlgorithm::FullBK;
            } else if (value == "kernel-bk") {
                options.subsumedAlgorithm = SubsumedCliqueAlgorithm::KernelBK;
            } else {
                usageError("unknown subsumed clique algorithm");
            }
        } else if (argument.starts_with("--subsumed-size-histogram=")) {
            options.subsumedSizeHistogram =
                parseOnOff(argument.substr(26), "subsumed size histogram");
        } else if (argument.starts_with("--somei-pruning=")) {
            options.someiPruning =
                parseOnOff(argument.substr(16), "SOMEi pruning");
        } else if (argument.starts_with("--somei-changes=")) {
            options.someiTrackChanges =
                parseOnOff(argument.substr(16), "SOMEi change tracking");
        } else if (argument.starts_with("--somei-profile=")) {
            options.someiProfilePhases =
                parseOnOff(argument.substr(16), "SOMEi profiling");
        } else if (argument.starts_with("--somei-build-update=")) {
            const auto value = argument.substr(21);
            if (value == "bk") {
                options.someiEnumerateBuild = true;
            } else if (value == "paper") {
                options.someiEnumerateBuild = false;
            } else {
                usageError("unknown SOMEi build update: " + std::string(value));
            }
        } else if (argument.starts_with("--somei-build=")) {
            const auto value = argument.substr(14);
            if (value == "batches") {
                options.someiOneBatchBuild = false;
            } else if (value == "one-batch") {
                options.someiOneBatchBuild = true;
            } else {
                usageError("unknown SOMEi build: " + std::string(value));
            }
        } else if (argument.starts_with("--change-output=")) {
            const auto value = argument.substr(16);
            if (value == "new") {
                options.changeOutput = ChangeOutputOption::New;
            } else if (value == "subsumed" || value == "deleted") {
                options.changeOutput = ChangeOutputOption::Subsumed;
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
            // alarm(2) takes an unsigned count of seconds.
            if (seconds > std::numeric_limits<unsigned>::max()) {
                usageError("--measured-time-limit is too large");
            }
            options.measuredTimeLimitSeconds = seconds;
        } else if (argument.starts_with("--measured-time-extend=")) {
            const std::size_t seconds = parsePositiveSize(
                argument.substr(23), "--measured-time-extend");
            if (seconds > std::numeric_limits<unsigned>::max()) {
                usageError("--measured-time-extend is too large");
            }
            options.measuredTimeExtendSeconds = seconds;
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
    if (options.warmupBatchSize.has_value() && options.warmupEdges == 0U) {
        usageError("--warmup-batch-size needs --warmup-edges");
    }
    if (options.warmupEdges != 0U && !options.batchSize.has_value()) {
        usageError("--warmup-edges needs --batch-size: the warm-up is cut by "
                   "edge count");
    }
    if (options.measuredTimeExtendSeconds.has_value() &&
        (!options.measuredTimeLimitSeconds.has_value() ||
         *options.measuredTimeExtendSeconds < *options.measuredTimeLimitSeconds)) {
        usageError("--measured-time-extend needs a --measured-time-limit no "
                   "larger than it");
    }
    if (options.maxEdges.has_value() && !options.batchSize.has_value()) {
        usageError("--max-edges needs --batch-size: the stream is cut by edge "
                   "count");
    }
    // Full maintenance (SOMEi included) builds the maximal-clique state the
    // measured batches start from out of the warm-up's own searches.
    if (!options.warmupSearch &&
        options.maintenance != MaintenanceOption::NewOnly) {
        usageError("--warmup-search=off requires --maintenance=new-only");
    }
    if (options.maintenance == MaintenanceOption::NewOnly &&
        options.changeOutput != ChangeOutputOption::New) {
        usageError("subsumed output requires --maintenance=full");
    }
    // SOMEi keeps every maximal clique by construction and updates a batch
    // as a whole, so it runs full maintenance in the batched pipeline only.
    if (options.framework == FrameworkOption::Somei) {
        if (options.maintenance != MaintenanceOption::Full) {
            usageError("--framework=somei maintains every maximal clique, so "
                       "it needs --maintenance=full");
        }
        if (options.maintenancePipeline ==
            MaintenancePipelineOption::Streaming) {
            usageError("--framework=somei updates whole batches, so it has no "
                       "streaming pipeline");
        }
        if (!options.someiTrackChanges && !options.suppressOutput) {
            usageError("--somei-changes=off lists no changes, so it needs "
                       "--no-output");
        }
    } else if (options.someiOneBatchBuild) {
        usageError("--somei-build is SOMEi's: the other frameworks build their "
                   "cliques batch by batch");
    } else if (options.someiEnumerateBuild.has_value()) {
        usageError("--somei-build-update is SOMEi's: the other frameworks "
                   "build their cliques by their own search");
    }
    // Deleting, the cliques the deletions create can only be told apart by
    // the index, and they are recovered batch by batch.
    if (options.updateMode != UpdateModeOption::Insert) {
        if (isBKComparison(options.subsumedAlgorithm)) {
            usageError("full-bk and kernel-bk require --update-mode=insert");
        }
        const std::string mode = std::string("--update-mode=") +
            updateModeOptionName(options.updateMode);
        if (options.maintenance != MaintenanceOption::Full) {
            usageError(mode + " maintains the cliques the deletions create, "
                       "so it needs --maintenance=full");
        }
        if (options.maintenancePipeline ==
            MaintenancePipelineOption::Streaming) {
            usageError(mode + " runs the batched pipeline only");
        }
    }
    if (options.updateMode != UpdateModeOption::RandomLabels &&
        (options.insertFraction.has_value() || options.labelSeed.has_value())) {
        usageError("--insert-fraction and --label-seed need "
                   "--update-mode=random-labels");
    }
    if (options.updateMode != UpdateModeOption::RandomLabels &&
        options.updateMode != UpdateModeOption::SlidingWindow &&
        options.mixedStrategy.has_value()) {
        usageError("--mixed-strategy needs --update-mode=random-labels or "
                   "sliding-window: only their batches both delete and insert");
    }
    if (options.updateMode != UpdateModeOption::SlidingWindow &&
        options.windowEdges.has_value()) {
        usageError("--window-edges needs --update-mode=sliding-window");
    }
    // The window slides by edge count, and its first W edges are its own
    // warm-up.
    if (options.updateMode == UpdateModeOption::SlidingWindow) {
        if (!options.batchSize.has_value()) {
            usageError("--update-mode=sliding-window needs --batch-size: the "
                       "window slides by edge count");
        }
        if (options.warmupEdges != 0U) {
            usageError("--update-mode=sliding-window builds its first window "
                       "itself, so it takes no --warmup-edges");
        }
    }
    const bool streamingIsPossible =
        options.maintenance == MaintenanceOption::Full &&
        options.updateMode == UpdateModeOption::Insert &&
        options.framework != FrameworkOption::Both &&
        options.framework != FrameworkOption::All &&
        options.framework != FrameworkOption::Somei;
    if (options.maintenancePipeline == MaintenancePipelineOption::Streaming &&
        !streamingIsPossible) {
        if (options.maintenance != MaintenanceOption::Full) {
            usageError(
                "--maintenance-pipeline=streaming requires --maintenance=full");
        }
        usageError(
            "--maintenance-pipeline=streaming needs a single framework, "
            "because cross-framework comparison has to materialize both "
            "result sets");
    }
    if (options.maintenancePipeline == MaintenancePipelineOption::Auto) {
        options.maintenancePipeline = streamingIsPossible
            ? MaintenancePipelineOption::Streaming
            : MaintenancePipelineOption::Batched;
    }
    return options;
}

[[nodiscard]] const char* frameworkOptionName(
    FrameworkOption framework) noexcept {
    switch (framework) {
    case FrameworkOption::Star:
        return "star";
    case FrameworkOption::Edge:
        return "edge";
    case FrameworkOption::Vertex:
        return "vertex";
    case FrameworkOption::All:
        return "all";
    case FrameworkOption::Somei:
        return "somei";
    case FrameworkOption::Both:
        break;
    }
    return "both";
}

[[nodiscard]] const char* maintenanceOptionName(
    MaintenanceOption maintenance) noexcept {
    return maintenance == MaintenanceOption::Full ? "full" : "new-only";
}

[[nodiscard]] const char* maintenancePipelineOptionName(
    MaintenancePipelineOption pipeline) noexcept {
    switch (pipeline) {
    case MaintenancePipelineOption::Batched:
        return "batched";
    case MaintenancePipelineOption::Streaming:
        return "streaming";
    case MaintenancePipelineOption::Auto:
        break;
    }
    return "auto";
}

[[nodiscard]] const char* changeOutputOptionName(
    ChangeOutputOption changeOutput) noexcept {
    switch (changeOutput) {
    case ChangeOutputOption::New:
        return "new";
    case ChangeOutputOption::Subsumed:
        return "subsumed";
    case ChangeOutputOption::Both:
        break;
    }
    return "both";
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

void printEdgeStatistics(
    const NewCliqueSearchStats& stats,
    double elapsedMilliseconds) {
    std::cerr
        << "edge: time-ms=" << elapsedMilliseconds
        << ", seeds=" << stats.seeds
        << ", recursive-calls=" << stats.recursiveCalls
        << ", pivot-calls=" << stats.pivotCalls
        << ", pivot-candidates=" << stats.pivotCandidatesEvaluated
        << ", leaves=" << stats.leaves
        << ", emitted=" << stats.emittedCliques
        << ", ownership-moves=" << stats.ownershipMoves
        << ", local-builds=" << stats.localGraph.builds
        << ", local-vertices=" << stats.localGraph.vertexOccurrences
        << ", local-edges=" << stats.localGraph.edgeOccurrences
        << ", max-local-vertices=" << stats.localGraph.maximumVertices
        << ", peak-local-bytes=" << stats.localGraph.peakBytes
        << ", snapshot-edge-lookups=" << stats.localGraph.snapshotEdgeIdLookups
        << '\n';
}

void printVertexStatistics(
    const VertexSearchStats& stats,
    double elapsedMilliseconds,
    VertexBranching branching,
    VertexSeeds seeds) {
    std::cerr
        << "vertex: time-ms=" << elapsedMilliseconds
        << ", branching=" << vertexBranchingName(branching)
        << ", seeds=" << vertexSeedsName(seeds)
        << ", root-edge-choices=" << stats.rootEdgeChoices
        << ", root-pivot-choices=" << stats.rootPivotChoices
        << ", unbound-edge-choices=" << stats.unboundEdgeChoices
        << ", unbound-pivot-choices=" << stats.unboundPivotChoices
        << ", degree-decisions=" << stats.degreeDecisions
        << ", unit-bindings=" << stats.unitBindings
        << ", edge-children=" << stats.edgeChildren
        << ", skipped-edge-children=" << stats.skippedEdgeChildren
        << ", exclusion-moves=" << stats.exclusionMoves
        << ", batches=" << stats.batches
        << ", batch-edges=" << stats.batchEdges
        << ", roots=" << stats.roots
        << ", pruned-roots=" << stats.prunedRoots
        << ", recursive-calls=" << stats.recursiveCalls
        << ", pivot-calls=" << stats.pivotCalls
        << ", pivot-candidates=" << stats.pivotCandidatesEvaluated
        << ", new-edge-prunes=" << stats.newEdgePrunes
        << ", batch-edge-checks=" << stats.batchEdgeChecks
        << ", leaves=" << stats.leaves
        << ", emitted=" << stats.emittedCliques
        << ", local-builds=" << stats.localGraph.builds
        << ", local-vertices=" << stats.localGraph.vertexOccurrences
        << ", local-edges=" << stats.localGraph.edgeOccurrences
        << ", max-local-vertices=" << stats.localGraph.maximumVertices
        << ", peak-local-bytes=" << stats.localGraph.peakBytes
        << '\n';
}

void printStarStatistics(
    const StarSearchStats& stats,
    double elapsedMilliseconds,
    StarBranchScope branchScope,
    StarBranching branching,
    bool rootFilter) {
    std::cerr
        << "star: time-ms=" << elapsedMilliseconds
        << ", branch-scope=" << starBranchScopeName(branchScope)
        << ", branching=" << starBranchingName(branching)
        << ", root-filter=" << (rootFilter ? "on" : "off")
        << ", reindexed-trees=" << stats.reindexedTrees
        << ", merged-reindexed-subtrees=" << stats.mergedReindexedSubtrees
        << ", deep-reindexed-subtrees=" << stats.deepReindexedSubtrees
        << ", node-words=" << stats.nodeWords
        << ", pivot-words=" << stats.pivotWords
        << ", branch-root-edge-stars=" << stats.branchRootEdgeStars
        << ", branch-root-pivot-stars=" << stats.branchRootPivotStars
        << ", branch-root-edge-branches=" << stats.branchRootEdgeBranches
        << ", branch-root-pivot-branches=" << stats.branchRootPivotBranches
        << ", branch-root-chosen-branches=" << stats.branchRootChosenBranches
        << ", branch-root-degree-decisions=" << stats.branchRootDegreeDecisions
        << ", branch-root-scans=" << stats.branchRootScans
        << ", branch-unbound-bindings=" << stats.branchUnboundBindings
        << ", branch-unit-bindings=" << stats.branchUnitBindings
        << ", singleton-stars=" << stats.singletonStars
        << ", merged-star-edges=" << stats.mergedStarEdges
        << ", per-partner-star-edges=" << stats.perPartnerStarEdges
        << ", singleton-star-edges=" << stats.singletonStarEdges
        << ", batches=" << stats.batches
        << ", batch-edges=" << stats.batchEdges
        << ", cover-vertices=" << stats.coverVertices
        << ", cover-both-endpoints=" << stats.coverBothEndpoints
        << ", stars=" << stats.stars
        << ", star-edges-per-star="
        << ratioOrNaN(
               static_cast<double>(stats.starEdgeOccurrences),
               static_cast<double>(stats.stars))
        << ", max-star-edges=" << stats.maximumStarEdges
        << ", domain-vertices=" << stats.domainVertexOccurrences
        << ", domain-vertices-per-star="
        << ratioOrNaN(
               static_cast<double>(stats.domainVertexOccurrences),
               static_cast<double>(stats.stars))
        << ", max-domain-vertices=" << stats.maximumDomainVertices
        << ", assign-ms=" << stats.assignMilliseconds
        << ", domain-scan-reads=" << stats.domainScanReads
        << ", domain-mark-writes=" << stats.domainMarkWrites
        << ", recursive-calls=" << stats.recursiveCalls
        << ", pivot-calls=" << stats.pivotCalls
        << ", pivot-candidates=" << stats.pivotCandidatesEvaluated
        << ", leaves=" << stats.leaves
        << ", partner-prunes=" << stats.partnerPrunes
        << ", partnerless-nodes=" << stats.partnerlessNodes
        << ", emitted=" << stats.emittedCliques
        << ", ownership-moves=" << stats.ownershipMoves
        << ", ownership-probes=" << stats.ownershipProbes
        << ", ownership-walks=" << stats.ownershipWalks
        << ", local-builds=" << stats.localGraph.builds
        << ", local-vertices=" << stats.localGraph.vertexOccurrences
        << ", local-edges=" << stats.localGraph.edgeOccurrences
        << ", max-local-vertices=" << stats.localGraph.maximumVertices
        << ", peak-local-bytes=" << stats.localGraph.peakBytes
        << '\n';
}

void printSubsumedStatistics(
    const SubsumedCliqueSearchStats& stats,
    double elapsedMilliseconds,
    SubsumedCliqueAlgorithm algorithm,
    MaintenancePipelineOption pipeline) {
    std::cerr
        << "subsumed: time-ms=" << elapsedMilliseconds
        // Streaming leaves this at zero: its subsumed work runs inside the
        // search timer, so `search-and-maintenance-ms` is the comparable one.
        << ", time-scope="
        << (pipeline == MaintenancePipelineOption::Streaming
                ? "streamed-into-search"
                : "subsumed")
        << ", algorithm=" << subsumedCliqueAlgorithmName(algorithm)
        << ", batches=" << stats.batches
        << ", new-cliques-processed=" << stats.newCliquesProcessed
        << ", internal-batch-edges=" << stats.internalBatchEdges
        << ", split-rounds=" << stats.splitRounds
        << ", candidate-occurrences=" << stats.candidateOccurrences
        << ", duplicate-candidates-removed="
        << stats.duplicateCandidatesRemoved
        << ", index-lookups=" << stats.indexLookups
        << ", batch-edge-visits=" << stats.batchEdgeVisits
        << ", bitmask-cliques=" << stats.bitmaskCliques
        << ", fallback-cliques=" << stats.fallbackCliques
        << ", size-gate-skipped-cliques=" << stats.sizeGateSkippedCliques
        << ", size-gate-skipped-probes=" << stats.sizeGateSkippedProbes
        << ", fast-path-cliques=" << stats.fastPathCliques
        << ", deferred-probes=" << stats.deferredProbes
        << ", inline-probes=" << stats.inlineProbes
        << ", structural-star-cliques=" << stats.structuralStarCliques
        << ", structural-complete-cliques=" << stats.structuralCompleteCliques
        << ", kernel-candidates=" << stats.kernelCandidates
        << ", max-deferred-probes=" << stats.maximumDeferredProbes
        << ", max-candidates=" << stats.maximumCandidates
        << ", emitted=" << stats.emittedCliques
        << ", residual-kernel-searches=" << stats.residualKernelSearches
        << ", residual-kernel-nodes=" << stats.residualKernelNodes
        << ", phase-prepare-ms=" << stats.phasePrepareMilliseconds
        << ", phase-cliques-ms=" << stats.phaseCliquesMilliseconds
        << ", phase-flush-ms=" << stats.phaseFlushMilliseconds
        << ", phase-finish-ms=" << stats.phaseFinishMilliseconds
        << '\n';
    // Shared BK engine counters by search width: <= 64, 65-256, and > 256.
    const wideresidual::Counters& wide = stats.wide;
    std::cerr << "wideresidual: nodes=" << wide.nodes
              << ", leaves=" << wide.leaves
              << ", kernel-calls=" << wide.kernelCalls
              << ", kernel-narrow=" << wide.kernelNarrow
              << ", kernel-wide=" << wide.kernelWide
              << ", kernel-dynamic=" << wide.kernelDynamic
              << ", dynamic-nodes=" << wide.dynamicNodes << '\n';
}

void printSomeiStatistics(
    const SomeiStats& stats,
    double elapsedMilliseconds,
    const SomeiIndex& index,
    bool oneBatchBuild) {
    const SomeiConfiguration& configuration = index.configuration();
    std::cerr
        << "somei: time-ms=" << elapsedMilliseconds
        << ", pruning=" << (configuration.pruning ? "on" : "off")
        << ", changes=" << (configuration.trackChanges ? "on" : "off")
        << ", profile=" << (configuration.profilePhases ? "on" : "off")
        << ", build=" << (oneBatchBuild ? "one-batch" : "batches")
        << ", build-update=" << (configuration.enumerateBuild ? "bk" : "paper")
        << ", batches=" << stats.batches
        << ", batch-edges=" << stats.batchEdges
        << ", scope-ms=" << stats.scopeMilliseconds
        << ", rebuild-ms=" << stats.rebuildMilliseconds
        << ", change-ms=" << stats.changeMilliseconds
        << ", rollback-ms=" << stats.rollbackMilliseconds
        << ", redo-ms=" << stats.redoMilliseconds
        << ", tree-add-ms=" << stats.treeAddMilliseconds
        << ", big-vertices=" << stats.bigVertices
        << ", scope-visits=" << stats.scopeVisits
        << ", pruned-visits=" << stats.prunedVisits
        << ", rebuilt-headers=" << stats.rebuiltHeaders
        << ", rolled-back-bodies=" << stats.rolledBackBodies
        << ", removed-bodies=" << stats.removedBodies
        << ", truncated-bodies=" << stats.truncatedBodies
        << ", redo-steps=" << stats.redoSteps
        << ", empty-renovations=" << stats.emptyRenovations
        << ", body-comparisons=" << stats.bodyComparisons
        << ", extended-bodies=" << stats.extendedBodies
        << ", created-bodies=" << stats.createdBodies
        << ", partial-intersections=" << stats.partialIntersections
        << ", tree-adds=" << stats.treeAdds
        << ", tree-deletes=" << stats.treeDeletes
        << ", tree-nodes-created=" << stats.treeNodesCreated
        << ", tree-nodes-deleted=" << stats.treeNodesDeleted
        << ", leaf-events=" << stats.leafEvents
        << ", new-cliques=" << stats.newCliques
        << ", deleted-cliques=" << stats.deletedCliques
        << ", candidate-cliques=" << index.candidateCliqueCount()
        << ", maximal-cliques=" << index.maximalCliqueCount()
        << ", tree-nodes=" << index.treeNodeCount()
        << '\n';
}

} // namespace
} // namespace immce

int main(int argc, char** argv) {
    using namespace immce;
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
        // The run uses the stream's first `runEdgeCount` edges.  The warm-up's
        // batches come first; `warmupBatchCount` is where the measured ones
        // begin.  A sliding window's warm-up is its first window, cut into the
        // run's own batches.
        const std::size_t streamEdgeCount = stream.insertionTimestamps().size();
        if (options.maxEdges.value_or(0U) > streamEdgeCount) {
            throw std::invalid_argument(
                "--max-edges is larger than the stream's " +
                std::to_string(streamEdgeCount) + " edges");
        }
        const std::size_t runEdgeCount =
            options.maxEdges.value_or(streamEdgeCount);
        const bool slidingWindow =
            options.updateMode == UpdateModeOption::SlidingWindow;
        const std::size_t windowEdges = slidingWindow
            ? options.windowEdges.value_or(runEdgeCount / 2U)
            : 0U;
        if (slidingWindow) {
            if (windowEdges < *options.batchSize) {
                throw std::invalid_argument(
                    "the " + std::to_string(windowEdges) +
                    "-edge window is smaller than a batch, which would delete "
                    "edges it has not inserted yet");
            }
            if (windowEdges >= runEdgeCount) {
                throw std::invalid_argument(
                    "the " + std::to_string(windowEdges) +
                    "-edge window leaves none of the run's " +
                    std::to_string(runEdgeCount) + " edges to slide over");
            }
        }
        const std::size_t warmupEdges =
            slidingWindow ? windowEdges : options.warmupEdges;
        const std::size_t warmupBatchSize = slidingWindow
            ? *options.batchSize
            : options.warmupBatchSize.value_or(options.batchSize.value_or(0U));
        if (warmupEdges > runEdgeCount) {
            throw std::invalid_argument(
                "--warmup-edges is larger than the run's " +
                std::to_string(runEdgeCount) + " edges");
        }
        std::size_t warmupBatchCount = 0U;
        std::vector<EdgeBatch> batches = measurePhase(
            [&] {
                std::vector<EdgeBatch> all;
                if (warmupEdges != 0U) {
                    all = stream.makeBatches(
                        warmupBatchSize, std::nullopt, 0U, warmupEdges);
                }
                warmupBatchCount = all.size();
                std::vector<EdgeBatch> measured = stream.makeBatches(
                    options.batchSize, options.maxBatches, warmupEdges,
                    options.maxEdges.value_or(
                        std::numeric_limits<std::size_t>::max()));
                all.insert(all.end(), std::make_move_iterator(measured.begin()),
                           std::make_move_iterator(measured.end()));
                return all;
            },
            batchBuildMilliseconds);
        const std::size_t availableBatchCount = options.batchSize.has_value()
            ? (runEdgeCount - warmupEdges + *options.batchSize - 1U) /
                  *options.batchSize
            : stream.batchCount(options.batchSize);
        std::size_t processedEdgeCount = 0;
        for (std::size_t index = warmupBatchCount; index < batches.size();
             ++index) {
            processedEdgeCount += batches[index].edgeIds.size();
        }
        IncrementalGraph graph = measurePhase(
            [&] {
                return IncrementalGraph(
                    stream.vertexCount(), stream.takeEdges());
            },
            graphInitializationMilliseconds);

        const bool runStar = options.framework == FrameworkOption::Star ||
            options.framework == FrameworkOption::Both ||
            options.framework == FrameworkOption::All;
        const bool runEdge = options.framework == FrameworkOption::Edge ||
            options.framework == FrameworkOption::Both ||
            options.framework == FrameworkOption::All;
        const bool runVertex = options.framework == FrameworkOption::Vertex ||
            options.framework == FrameworkOption::All;
        const bool runSomei = options.framework == FrameworkOption::Somei;
        const bool fullMaintenance =
            options.maintenance == MaintenanceOption::Full;
        const bool streamingMaintenance = fullMaintenance &&
            options.maintenancePipeline ==
                MaintenancePipelineOption::Streaming;
        const bool streamNewOutput = !options.suppressOutput &&
            (options.changeOutput == ChangeOutputOption::New ||
             options.changeOutput == ChangeOutputOption::Both);
        const bool streamSubsumedOutput = !options.suppressOutput &&
            (options.changeOutput == ChangeOutputOption::Subsumed ||
             options.changeOutput == ChangeOutputOption::Both);
        // A single producer: with two or three, the batches are materialized
        // and compared.
        const bool countOnly = !fullMaintenance && options.suppressOutput &&
            (runStar ? 1 : 0) + (runEdge ? 1 : 0) + (runVertex ? 1 : 0) == 1;

        const bool reverseDelete =
            options.updateMode == UpdateModeOption::ReverseDelete;
        const bool randomLabels =
            options.updateMode == UpdateModeOption::RandomLabels;
        const double insertFraction = options.insertFraction.value_or(0.5);
        const std::uint64_t labelSeed = options.labelSeed.value_or(1U);
        // Das et al.'s two phases by default (2026-10-11): star-rk's single
        // pass came out slower.
        const MixedStrategy mixedStrategy =
            options.mixedStrategy.value_or(MixedStrategy::TwoPhase);

        MaintenanceConfiguration configuration;
        configuration.runStar = runStar;
        configuration.runEdge = runEdge;
        configuration.runVertex = runVertex;
        configuration.runSomei = runSomei;
        configuration.fullMaintenance = fullMaintenance;
        configuration.builds = reverseDelete || randomLabels || slidingWindow;
        configuration.starBranchScope = options.starBranchScope;
        configuration.starBranching = options.starBranching;
        configuration.starRootFilter = options.starRootFilter;
        configuration.vertexBranching = options.vertexBranching;
        configuration.vertexSeeds = options.vertexSeeds;
        configuration.subsumedAlgorithm = options.subsumedAlgorithm;
        // Only the framework that reads the histogram pays to maintain it,
        // so the control is not taxed for a structure it never consults.
        // Deleting alone, nothing reads it: the cliques the deletions create
        // are not in the index.
        configuration.sizeHistogram = options.subsumedSizeHistogram.value_or(
            !reverseDelete && usesSubsumedSizeGate(options.subsumedAlgorithm));
        configuration.somei.pruning = options.someiPruning;
        configuration.somei.trackChanges = options.someiTrackChanges;
        configuration.somei.profilePhases = options.someiProfilePhases;
        configuration.somei.enumerateBuild =
            options.someiEnumerateBuild.value_or(true);
        configuration.mixedStrategy = mixedStrategy;
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

        // Random labels: every edge of the measured window is labelled insert
        // with probability `insertFraction`, in stream order, by a 64-bit
        // Mersenne Twister (whose output, unlike the standard distributions',
        // is the same under every standard library).  `insertedOf` and
        // `deletedOf` hold each measured batch's edges by label.
        std::vector<std::vector<EdgeId>> insertedOf;
        std::vector<std::vector<EdgeId>> deletedOf;
        if (randomLabels) {
            std::mt19937_64 generator(labelSeed);
            for (std::size_t index = warmupBatchCount; index < batches.size();
                 ++index) {
                std::vector<EdgeId>& inserted = insertedOf.emplace_back();
                std::vector<EdgeId>& deleted = deletedOf.emplace_back();
                for (const EdgeId id : batches[index].edgeIds) {
                    const double draw =
                        static_cast<double>(generator() >> 11U) * 0x1.0p-53;
                    (draw < insertFraction ? inserted : deleted).push_back(id);
                }
            }
        }

        // The end of the warm-up: whatever it measured is dropped, and its
        // own duration is kept for the record.
        double warmupMilliseconds = 0.0;
        const auto warmupBegin = std::chrono::steady_clock::now();
        const auto endWarmup = [&] {
            warmupMilliseconds = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - warmupBegin).count();
            engine.resetMeasurements();
        };

        // The order the batches run in.  Inserting: each batch once, forward.
        // Reverse delete: every batch as an untimed insertion (the build),
        // then the measured batches deleted, last first.  Random labels: the
        // warm-up and each measured batch's delete-labelled edges as the
        // build, then the measured batches, each inserting and deleting at
        // once.  Sliding window: the first window as the build, then the
        // measured batches, each inserting its edges and deleting as many of
        // the oldest.  The steps from `measuredBegin` on are the measured ones.
        enum class StepKind { Build, Insert, Delete, Update, Slide };
        struct Step {
            StepKind kind = StepKind::Insert;
            std::size_t batch = 0;
        };
        const bool builds = reverseDelete || randomLabels;
        std::vector<Step> steps;
        steps.reserve(batches.size() * (builds ? 2U : 1U));
        for (std::size_t index = 0; index < batches.size(); ++index) {
            StepKind kind = builds ? StepKind::Build : StepKind::Insert;
            if (slidingWindow) {
                kind = index < warmupBatchCount ? StepKind::Build
                                                : StepKind::Slide;
            }
            steps.push_back(Step{kind, index});
        }
        if (reverseDelete) {
            for (std::size_t index = batches.size(); index > warmupBatchCount;
                 --index) {
                steps.push_back(Step{StepKind::Delete, index - 1U});
            }
        } else if (randomLabels) {
            for (std::size_t index = warmupBatchCount; index < batches.size();
                 ++index) {
                steps.push_back(Step{StepKind::Update, index});
            }
        }
        const std::size_t measuredBegin =
            builds ? batches.size() : warmupBatchCount;

        // The batched pipeline's order and format: a batch's new cliques,
        // then the ones it removed, each sorted, so that change streams
        // compare byte for byte between frameworks and pipelines.
        BatchChanges changes;
        const auto writeChanges = [&](int timestamp) {
            if (options.suppressOutput) {
                return;
            }
            std::ranges::sort(changes.created, CliqueOutputOrder{});
            std::ranges::sort(changes.removed, CliqueOutputOrder{});
            const auto writeClique = [&](std::string_view change,
                                         const Clique& clique) {
                *output << timestamp << '\t';
                if (fullMaintenance) {
                    *output << change << '\t';
                }
                for (std::size_t index = 0; index < clique.size(); ++index) {
                    if (index != 0) {
                        *output << ' ';
                    }
                    *output << stream.originalVertexId(clique[index]);
                }
                *output << '\n';
            };
            if (!fullMaintenance ||
                options.changeOutput == ChangeOutputOption::New ||
                options.changeOutput == ChangeOutputOption::Both) {
                for (const Clique& clique : changes.created) {
                    writeClique("new", clique);
                }
            }
            if (fullMaintenance &&
                (options.changeOutput == ChangeOutputOption::Subsumed ||
                 options.changeOutput == ChangeOutputOption::Both)) {
                for (const Clique& clique : changes.removed) {
                    writeClique("del", clique);
                }
            }
        };
        const StreamTarget streamTarget{
            *output, stream, streamNewOutput, streamSubsumedOutput};
        // A sliding batch's expired edges: edge ids are stream positions, so
        // they are the batch's own positions moved back by the window.
        std::vector<EdgeId> expired;
        // `--somei-build=one-batch`: the build steps' edges, gathered for one
        // build after the last of them (every build step precedes the
        // measured ones).
        std::vector<EdgeId> oneBatchBuild;
        // --measured-time-extend: when the last measured batch ended, and the
        // time the batches are now allowed (seconds since they began).
        timespec lastMeasuredStepEnd{};
        double measuredAllowedSeconds =
            static_cast<double>(options.measuredTimeLimitSeconds.value_or(0U));
        const auto secondsBetween = [](const timespec& from, const timespec& to) {
            return static_cast<double>(to.tv_sec - from.tv_sec) +
                static_cast<double>(to.tv_nsec - from.tv_nsec) * 1e-9;
        };
        // The timing line's search-and-maintenance total so far.
        const auto searchAndMaintenanceMilliseconds = [&] {
            const MaintenanceTimes& times = engine.times();
            if (streamingMaintenance) {
                return times.streamOuter + times.singleVertex;
            }
            return times.star + times.edge + times.vertex + times.subsumed +
                times.cliqueIndexUpdate + times.somei + times.transientNetting +
                times.singleVertex;
        };
        const auto wallMilliseconds = [&] {
            return std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - programBegin).count();
        };

        for (std::size_t stepIndex = 0; stepIndex < steps.size(); ++stepIndex) {
            if (stepIndex == measuredBegin && measuredBegin != 0U) {
                endWarmup();
            }
            if (stepIndex == measuredBegin && options.printProgress) {
                std::cerr << "progress: measured-begin, measured-batches=0/"
                          << steps.size() - measuredBegin
                          << ", wall-ms=" << wallMilliseconds() << std::endl;
            }
            if (stepIndex == measuredBegin &&
                options.measuredTimeLimitSeconds.has_value()) {
                ::clock_gettime(CLOCK_MONOTONIC, &immceMeasuredBegin);
                lastMeasuredStepEnd = immceMeasuredBegin;
                std::signal(SIGALRM, immceStopAtMeasuredTimeLimit);
                ::alarm(static_cast<unsigned>(*options.measuredTimeLimitSeconds));
            }
            const Step& step = steps[stepIndex];
            const EdgeBatch& batch = batches[step.batch];
            switch (step.kind) {
            case StepKind::Build: {
                const std::span<const EdgeId> edges =
                    randomLabels && step.batch >= warmupBatchCount
                        ? std::span<const EdgeId>(
                              deletedOf[step.batch - warmupBatchCount])
                        : std::span<const EdgeId>(batch.edgeIds);
                if (options.someiOneBatchBuild) {
                    oneBatchBuild.insert(
                        oneBatchBuild.end(), edges.begin(), edges.end());
                    const bool lastBuild = stepIndex + 1U == steps.size() ||
                        steps[stepIndex + 1U].kind != StepKind::Build;
                    if (lastBuild && !oneBatchBuild.empty()) {
                        engine.build(oneBatchBuild);
                        oneBatchBuild = {};
                    }
                } else if (!edges.empty()) {
                    engine.build(edges);
                }
                break;
            }
            case StepKind::Delete:
                engine.updateBatch(batch.timestamp, {}, batch.edgeIds, changes);
                writeChanges(batch.timestamp);
                break;
            case StepKind::Update: {
                const std::size_t measured = step.batch - warmupBatchCount;
                engine.updateBatch(
                    batch.timestamp, insertedOf[measured], deletedOf[measured],
                    changes);
                writeChanges(batch.timestamp);
                break;
            }
            case StepKind::Slide: {
                const auto first = static_cast<std::size_t>(batch.edgeIds.front());
                expired.clear();
                for (std::size_t offset = 0; offset < batch.edgeIds.size();
                     ++offset) {
                    expired.push_back(
                        static_cast<EdgeId>(first - windowEdges + offset));
                }
                engine.updateBatch(
                    batch.timestamp, batch.edgeIds, expired, changes);
                writeChanges(batch.timestamp);
                break;
            }
            case StepKind::Insert:
                if (step.batch < warmupBatchCount && !options.warmupSearch) {
                    engine.applyOnly(batch.edgeIds);
                } else if (countOnly) {
                    engine.countBatch(batch.edgeIds);
                } else if (streamingMaintenance) {
                    engine.streamBatch(batch, streamTarget);
                } else {
                    engine.updateBatch(
                        batch.timestamp, batch.edgeIds, {}, changes);
                    writeChanges(batch.timestamp);
                }
                break;
            }
            // --measured-time-extend, judged after every measured batch but
            // the last: the batches may run until the extended limit while
            // the time so far plus the last batch's time for every batch
            // left stays within it, and until the plain limit otherwise.
            const bool judging = options.measuredTimeExtendSeconds.has_value() &&
                stepIndex >= measuredBegin && stepIndex + 1U < steps.size();
            double projectedSeconds = 0.0;
            if (judging) {
                timespec now{};
                ::clock_gettime(CLOCK_MONOTONIC, &now);
                const double elapsed = secondsBetween(immceMeasuredBegin, now);
                const double last = secondsBetween(lastMeasuredStepEnd, now);
                lastMeasuredStepEnd = now;
                projectedSeconds = elapsed +
                    last * static_cast<double>(steps.size() - stepIndex - 1U);
                measuredAllowedSeconds = static_cast<double>(
                    projectedSeconds <=
                            static_cast<double>(*options.measuredTimeExtendSeconds)
                        ? *options.measuredTimeExtendSeconds
                        : *options.measuredTimeLimitSeconds);
            }
            if (options.printProgress && stepIndex >= measuredBegin) {
                std::cerr << "progress: measured-batches="
                          << stepIndex + 1U - measuredBegin << '/'
                          << steps.size() - measuredBegin
                          << ", search-and-maintenance-ms="
                          << searchAndMaintenanceMilliseconds()
                          << ", wall-ms=" << wallMilliseconds();
                if (judging) {
                    std::cerr << ", projected-measured-s=" << projectedSeconds
                              << ", allowed-measured-s=" << measuredAllowedSeconds;
                }
                std::cerr << std::endl;
            }
            if (judging) {
                timespec now{};
                ::clock_gettime(CLOCK_MONOTONIC, &now);
                const double elapsed = secondsBetween(immceMeasuredBegin, now);
                if (measuredAllowedSeconds <= elapsed) {
                    immceStopAtMeasuredTimeLimit(0);
                }
                ::alarm(static_cast<unsigned>(
                    std::ceil(measuredAllowedSeconds - elapsed)));
            }
        }
        if (options.measuredTimeLimitSeconds.has_value()) {
            ::alarm(0U);
        }
        // A warm-up with nothing after it still measures nothing.
        if (measuredBegin != 0U && measuredBegin == steps.size()) {
            endWarmup();
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
                << ", batches=" << batches.size() - warmupBatchCount
                << ", available-batches=" << availableBatchCount
                << ", max-batches=";
            if (options.maxBatches.has_value()) {
                std::cerr << *options.maxBatches;
            } else {
                std::cerr << "all";
            }
            // Everything on this line and below covers the measured batches
            // only; the warm-up (a sliding window's first window) is reported
            // here and nowhere else.
            std::cerr
                << ", warmup-edges=" << warmupEdges
                << ", warmup-batches=" << warmupBatchCount
                << ", warmup-batch-size="
                << (warmupEdges == 0U ? std::size_t{0} : warmupBatchSize)
                << ", warmup-ms=" << warmupMilliseconds;
            // The three counts include the single-vertex cliques, which the
            // next three break out (see MaintenanceEngine).
            const UpdateCounts& changeCounts = engine.updateCounts();
            std::cerr
                << ", new-maximal-cliques=" << engine.totalCreated()
                << ", subsumed-maximal-cliques=" << engine.totalRemoved()
                << ", current-maximal-cliques=" << engine.currentCliqueCount()
                << ", new-single-vertex-cliques="
                << changeCounts.singleVertexCreated
                << ", subsumed-single-vertex-cliques="
                << changeCounts.singleVertexRemoved
                << ", isolated-vertices=" << engine.isolatedVertexCount()
                << ", batch-mode="
                << (options.batchSize.has_value() ? "fixed-size" : "timestamp")
                << ", batch-size=";
            if (options.batchSize.has_value()) {
                std::cerr << *options.batchSize;
            } else {
                std::cerr << "timestamp";
            }
            // Stated so that a run can be checked against the batches it means.
            std::cerr << ", tie-order=file";
            std::cerr
                << ", warmup-search=" << (options.warmupSearch ? "on" : "off")
                << ", framework=" << frameworkOptionName(options.framework)
                << ", maintenance="
                << maintenanceOptionName(options.maintenance)
                << ", update-mode=" << updateModeOptionName(options.updateMode)
                << ", maintenance-pipeline="
                << maintenancePipelineOptionName(options.maintenancePipeline)
                << ", subsumed-algorithm="
                << (runSomei
                        ? "none"
                        : subsumedCliqueAlgorithmName(options.subsumedAlgorithm))
                << ", change-output="
                << changeOutputOptionName(options.changeOutput)
                << ", result-mode="
                << (countOnly ? "count-only" : "materialized")
                << ", measured-time-limit=";
            if (options.measuredTimeLimitSeconds.has_value()) {
                std::cerr << *options.measuredTimeLimitSeconds;
            } else {
                std::cerr << "none";
            }
            std::cerr << ", measured-time-extend=";
            if (options.measuredTimeExtendSeconds.has_value()) {
                std::cerr << *options.measuredTimeExtendSeconds;
            } else {
                std::cerr << "none";
            }
            std::cerr << '\n';
            const double excludedFromSearchMilliseconds =
                inputPreprocessMilliseconds + batchBuildMilliseconds +
                graphInitializationMilliseconds + times.graphUpdate;
            const double maintenanceMilliseconds = times.subsumed +
                times.cliqueIndexUpdate + times.somei + times.transientNetting +
                times.singleVertex;
            std::cerr
                << "timing: scope=search-only"
                << ", input-preprocess-ms=" << inputPreprocessMilliseconds
                << ", batch-build-ms=" << batchBuildMilliseconds
                << ", graph-init-ms=" << graphInitializationMilliseconds
                << ", graph-update-ms=" << times.graphUpdate
                << ", clique-index-update-ms=" << times.cliqueIndexUpdate
                // SOMEi has no search apart from its maintenance: its whole
                // time is maintenance, and the search-and-maintenance total.
                << ", full-maintenance-ms=" << maintenanceMilliseconds
                // Streaming folds the subsumed search and the index update
                // into the search timer, so this sum is the only figure that
                // compares the two pipelines directly.  Under streaming it is
                // the outer per-batch timer, which additionally covers the
                // batch preparation and the closing flush -- neither of which
                // the visit timer sees.
                << ", search-and-maintenance-ms="
                << searchAndMaintenanceMilliseconds()
                // Streaming only: what the outer timer covers and the visit
                // timer does not -- `prepareBatch` plus the held-back probe
                // flush in `finishBatch`.
                << ", stream-prepare-flush-ms="
                << (streamingMaintenance
                        ? times.streamOuter - times.star - times.edge -
                              times.vertex
                        : 0.0)
                << ", search-time-scope="
                << (streamingMaintenance || runSomei ? "search+maintenance"
                                                     : "search")
                << ", stream-output-ms=" << times.streamOutput
                // Reading the batch endpoints' degrees for the single-vertex
                // cliques; counted in the two sums above.
                << ", single-vertex-ms=" << times.singleVertex
                << ", excluded-from-search-ms="
                << excludedFromSearchMilliseconds
                << '\n';
            if (runStar) {
                printStarStatistics(
                    engine.starStats(), times.star, options.starBranchScope,
                    options.starBranching, options.starRootFilter);
            }
            if (runEdge) {
                printEdgeStatistics(engine.edgeStats(), times.edge);
            }
            if (runVertex) {
                printVertexStatistics(
                    engine.vertexStats(), times.vertex, options.vertexBranching,
                    options.vertexSeeds);
            }
            if (runSomei) {
                printSomeiStatistics(
                    engine.someiStats(), times.somei, *engine.somei(),
                    options.someiOneBatchBuild);
            }
            const UpdateCounts& counts = engine.updateCounts();
            if (fullMaintenance && !runSomei) {
                // Live clique sizes: the size gate reads this histogram, so a
                // run that wants to know why it fired (or did not) needs it.
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
                std::cerr << char{10};
                printSubsumedStatistics(
                    engine.subsumedStats(),
                    times.subsumed,
                    options.subsumedAlgorithm,
                    options.maintenancePipeline);
                if (reverseDelete || randomLabels || slidingWindow) {
                    // The subsumed search run backwards for the deletions:
                    // deleting alone, its time and counters above are the
                    // recovery's.
                    const SubsumedCliqueSearchStats::Recovery& recovery =
                        engine.subsumedStats().recovery;
                    std::cerr << "recovery: destroyed=" << counts.destroyed
                              << ", created=" << counts.recovered
                              << ", candidates=" << recovery.candidates
                              << ", below-two=" << recovery.belowTwo
                              << ", duplicates=" << recovery.duplicates
                              << ", checks=" << recovery.checks
                              << ", rejected=" << recovery.rejected
                              << ", accepted=" << recovery.accepted
                              << ", neighbor-reads=" << recovery.neighborReads
                              << ", adjacency-tests=" << recovery.adjacencyTests
                              << ", groups=" << recovery.groups
                              << ", grouped-checks=" << recovery.groupedChecks
                              << ", row-reads=" << recovery.rowReads
                              << ", masked-searches=" << recovery.maskedSearches
                              << ", core-extenders=" << recovery.coreExtenders
                              << ", extension-masks=" << recovery.extensionMasks
                              << ", mask-cuts=" << recovery.maskCuts
                              << ", extended-whole=" << recovery.extendedWhole
                              << '\n';
                }
            }
            if (randomLabels || slidingWindow) {
                // Where the batches' changes come from (the labels, or the
                // window), and the times above split by side: the producers
                // on the deleted and on the inserted edges, the subsumed
                // search (insertions) and the recovery (deletions), and
                // two-phase's netting.  SOMEi has no such parts; only its
                // edge counts are filled in.
                std::cerr << "mixed: strategy="
                          << (runSomei ? "somei" : mixedStrategyName(mixedStrategy));
                if (randomLabels) {
                    std::cerr << ", insert-fraction=" << insertFraction
                              << ", label-seed=" << labelSeed;
                } else {
                    std::cerr << ", window-edges=" << windowEdges;
                }
                std::cerr << ", batches=" << counts.batches
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
                          << ", destroyed-side-recursive-calls="
                          << counts.destroyedSideCalls
                          << ", new-side-recursive-calls=" << counts.newSideCalls
                          << '\n';
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
