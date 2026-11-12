#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=immce_edge_vertex_crossover
#SBATCH --time=24:00:00
#SBATCH --output=/home/zzhong/zzy_project/DynamicClique/IMMCE_SubsumedClique/output/edge_vertex_crossover_%j.out
#SBATCH --error=/home/zzhong/zzy_project/DynamicClique/IMMCE_SubsumedClique/output/edge_vertex_crossover_%j.err
# OzSTAR refuses --exclusive; to keep the node to this job, ask for all its
# cores and memory on the sbatch line instead (see NODE below).

# Where the vertex baseline overtakes the edge one, and whether the method,
# star, stays below both: the new-clique search of the three (the paper's
# EdgeSeed, VertexSeed and CoverSeed), each one's running time and the number
# of branches its search takes, over a grid of batch sizes per dataset placed
# around the crossing of the two baselines.  The second version of
# submit_vertex_crossover_times.sh (job 18678150), from which this script was
# made (2026-10-21): there the crossing fell at the last one or two sizes of
# most grids, between coarse steps.  Here each grid is denser near its
# crossing and runs past it.  It calls no other script.
#
# WHY THEY CROSS.  edge roots one search per batch edge, so its time barely
# moves with the batch size.  vertex roots one search per batch endpoint over
# the endpoint's whole neighbourhood; the larger the batch, the more batch
# edges one such graph covers, so its time falls as batches grow.  Below the
# crossing edge is the faster baseline, above it vertex.  star puts the batch
# edges into vertex-cover stars and picks, per state, between branching on
# the star's new edges (as edge does) and on the pivot (as vertex does).
#
# THE ARMS:
#
#   edge        the edge baseline: one search per batch edge over the common
#               neighbourhood of its endpoints
#   vertex      the vertex baseline: one search per batch endpoint over its
#               closed neighbourhood N[s], the classic Tomita BK with the
#               new-edge condition pruning it (time-limited, see TIME LIMIT)
#   star        the method: vertex-cover stars, one local graph per star, the
#               branch-count rule at every state whose clique holds no new
#               edge yet (--star-branch-scope=unbound, the code default)
#
# and, not run by default,
#
#   star-root   the same rule asked at each star's centre only
#               (--star-branch-scope=root)
#
# METHOD_ARM (default star) names the arm the others are compared with; the
# crossover tables set it against both baselines.  Run without star
# (ARMS="edge vertex" METHOD_ARM=edge), they compare the two baselines only.
# Lists take entries separated by spaces or commas.
#
# 1. THE BRANCHES, per method and cell, from the statistics of its runs:
#   branches       recursive calls: every state (R, P, X) the search entered
#                  over the measured batches, roots and nodes cut on entry
#                  included -- the tree size submit_star_branch_tree_times.sh
#                  reports as search_nodes.  Roots count: how a method seeds
#                  its searches is part of how it branches, and edge, which
#                  roots one search per batch edge, would look small without
#                  them next to star, which branches on the same edges below
#                  one root
#   search_roots   the trees: edge one per batch edge; vertex one per root that
#                  passed the new-edge test before its graph was built; star
#                  one per single-edge star, one per partner of a star whose
#                  centre took the new edges, one per star whose centre took
#                  the pivot
#   branches_below_roots  branches - search_roots
#   leaves, emitted, pivot_calls, pivot_candidates, and per method: vertex's
#                  roots cut before their graphs were built (pruned_roots; not
#                  branches, no search was entered) and nodes cut by the
#                  new-edge condition (new_edge_prunes); star's subtrees cut
#                  for holding no partner (partner_prunes)
# The search is deterministic, so the runs of one arm must agree on it;
# tree.tsv says whether they did.
#
# 2. THE RUNNING TIME: `<method>: time-ms` (star / edge / vertex), the
# new-clique search alone on the count-only path (--maintenance=new-only
# --no-output).  Reading, batching and graph updates are outside it.
#
# COMPARISONS (results.tsv, geomeans.tsv): time(from)/time(to) and
# branches(from)/branches(to), above 1 meaning `to` is faster / opens fewer
# branches:
#   edge_vs_star     edge -> star
#   vertex_vs_star   vertex -> star
#   vertex_vs_edge   vertex -> edge: below 1, vertex is the faster baseline
# with star-root run, the first two for it as well and star_root_vs_star.
# EXPERIMENTS (name=from:to entries) replaces the list.  A cell's speedup is
# the median of its paired per-repeat ratios.
#
# TIME LIMIT.  The vertex baseline builds one bitmap over N[s] per root,
# |N[s]|^2/8 bytes, and a high-degree vertex is an endpoint in almost every
# batch: at the largest degree of each whole dataset that is 6.6 MB on
# coauthor, 146 MB on flickr, 1.05 GB on youtube and 6.4 GB on wikipedia,
# rebuilt every batch.  As SOMEi in submit_somei_times.sh, every vertex run is
# therefore killed after VERTEX_TIMEOUT_SECONDS (wall time, default 3600) and
# recorded as a timeout, which is a result: the speedup over it is then a
# lower bound,
#
#   (VERTEX_TIMEOUT_SECONDS - outside) / time(to)
#
# where `outside` is the largest wall-minus-search time of a finished run in
# the same cell (reading, batching and graph updates, which every arm does
# alike; a warm-up only inserts edges, so it is alike too).  The other arms
# run under RUN_TIMEOUT_SECONDS (default 7200, a cap no earlier job of these
# arms reached), and a timeout of theirs is treated the same way.
#
# BATCH ORDER.  The vertex baseline gets faster as batches grow (a root's graph
# is built once per batch however many of its batch edges there are; locally
# 0.44-0.55 of the b1000 time at b20000), so batch sizes run from the largest
# down, and with SKIP_SMALLER_AFTER_TIMEOUT=1 (default) an arm that timed out
# on a dataset is not run again there at that or any smaller batch size:
# those runs are written as skipped_after_timeout and enter no statistic.
#
# OUT OF MEMORY.  A run that exhausts memory (std::bad_alloc, SLURM's oom-kill,
# or a kill signal before its limit) is recorded as out_of_memory: a result
# like a timeout, not a failure of the job, and no reason to skip the smaller
# batches.  The 64 GB of the job hold the largest vertex bitmap above.
#
# OUTPUTS (in OUTPUT_DIR):
#   table_time.tsv   per cell: each method's median time (ms) -- ">=" a lower
#                    bound when it timed out -- and the speedups
#   table_tree.tsv   per cell: each method's branches and the branch ratios
#   times.tsv        per cell and method: the status of every run, every
#                    run's time, the median, min, max and spread
#   tree.tsv         per cell and method: the search statistics above, the
#                    branch ratio to METHOD_ARM, and whether the runs agreed
#   results.tsv      per cell and comparison: status, times, speedup (or its
#                    lower bound), branches and their ratio
#   geomeans.tsv     per comparison, per batch size within the small and the
#                    large datasets, per group and overall: the geometric mean
#                    of the speedups over the finished cells, a lower bound over
#                    the finished and the timed-out cells together, and the
#                    geometric mean of the branch ratios
#   crossover.tsv    per dataset, batch sizes ascending: the batches run, the
#                    median times of edge and vertex, vertex_vs_edge (as in
#                    results.tsv: time(vertex)/time(edge)) and the faster
#                    baseline; with a third arm run as METHOD_ARM, its time,
#                    edge_vs_ and vertex_vs_ it, and the fastest of the three
#   crossover_summary.tsv  per dataset: vertex_vs_edge at its smallest and its
#                    largest batch size, the last batch size at which edge was
#                    faster and the first at which vertex was, the crossing
#                    estimated between those two (where vertex_vs_edge
#                    reaches 1, on log scales), how many sizes each baseline
#                    won, and whether they crossed inside the grid (yes;
#                    yes_not_monotone if edge was faster again above the
#                    first vertex win; no; vertex_faster_throughout); with a
#                    third arm, in how many cells it was the fastest and its
#                    smallest edge_vs_ and vertex_vs_ speedups.  One run per
#                    cell: a ratio within a few percent of 1 is a tie
#   noise.tsv        per cell and repeat: graph-update-ms over the arms, the
#                    same work in every run; a spread well above 1 marks a
#                    disturbed run
#   timeouts.tsv, summary.tsv (one row per run), manifest.txt, datasets.txt,
#   source_sha256.txt, build.log, probe.log
#
# TIE ORDER.  --tie-order=file: edges go in by increasing timestamp, those
# sharing a timestamp in the dataset's line order -- the only order the code
# has.  Every run checks the order its statistics report.
#
# SETTINGS.  Since the 2026-10-09 cleanup the code has no other switch for
# these searches: the star partition is the vertex cover, and factored root
# preparation, merged and deep re-indexing are always on; since 2026-10-10 the
# edge baseline reads the edge ids its ownership test needs from the global
# adjacency, the dense table and --local-edge-index removed.  Every run checks
# its framework, branch scope, tie order and the count-only path from its
# statistics.
#
# CORRECTNESS.  No certification pass.  Every arm enumerates the same new
# maximal cliques, so inside one cell and one repeat the arms that finish must
# report the same count; a later one that disagrees fails the run.
#
# BATCH SIZES, one list per dataset, placed around the crossing that 18678150
# measured (bryan8), denser near it and running past it.  vertex_vs_edge =
# time(vertex)/time(edge), below 1 once vertex is faster.  The values marked
# "~" are predicted from 18678150 by grid_model.py, in
# experiments/edge_vertex_crossover_20261021 (edge's time interpolated on log
# scales, vertex's as its batches times its per-batch cost); the rest were
# measured there:
#
#   dataset (window)       crossing  sizes and vertex_vs_edge
#   dblp-cite (49,636)     ~9k       3k ~1.45, 4k ~1.36, 5k 1.28, 6k ~1.20,
#                                    8k ~1.06, 10k 0.95, 12.5k ~0.87,
#                                    15k ~0.84, 20k 0.72, 25k ~0.66
#   facebook (335,708)     ~170k     75k ~1.62, 100k 1.32, 125k ~1.24,
#                                    150k 1.11, 175k ~0.95, 200k 0.75,
#                                    250k ~0.74, 300k ~0.57
#   coauthor (8,344,615)   ~320k     50k 1.66, 100k 1.36, 200k 1.19,
#                                    250k ~1.10, 300k ~1.03, 400k ~0.92,
#                                    500k 0.85, 750k ~0.73, 1M 0.65,
#                                    2M ~0.52, 4M ~0.41
#   flickr (first 1M)      ~175k     75k ~1.56, 100k 1.64, 125k ~1.31,
#                                    150k ~1.15, 175k ~0.98, 200k 0.81,
#                                    250k 0.52, 300k ~0.62, 400k ~0.51,
#                                    500k 0.37
#   sx-stackoverflow       ~300k     150k ~1.66, 200k 1.37, 250k 1.16,
#     (first 1M)                     300k ~1.01, 350k ~0.81, 400k ~0.76,
#                                    500k 0.50, 600k ~0.48, 750k ~0.45
#   youtube (9,375,374)    ~1.75M    750k ~1.44, 1M 1.31, 1.25M ~1.19,
#                                    1.5M ~1.09, 2M 0.94, 2.5M ~0.85,
#                                    3M ~0.77, 3.5M ~0.72, 4.5M ~0.62
#   wikipedia (first 12M)  ~4.5M     2M 1.74, 2.5M ~1.51, 3M ~1.31,
#                                    3.5M ~1.23, 4M 1.08, 4.5M ~0.99,
#                                    5M ~0.91, 6M 0.69, 8M ~0.60, 10M ~0.48
#
# On log scales each crossing sits at 42-59% of its dataset's range, with at
# least three sizes on either side.  The smaller sizes of 18678150, where
# vertex took 1.5-7 times edge's time (wikipedia b50000 alone 49 minutes),
# are left to the usual grid of the other scripts.  The lists are
# DBLP_CITE_, FACEBOOK_, COAUTHOR_, FLICKR_, WIKIPEDIA_, YOUTUBE_ and
# STACKOVERFLOW_BATCH_SIZES; DEFAULT_BATCH_SIZES serves any other dataset
# (DATASET_FILTER=all), and BATCH_SIZES, if set, replaces them all.
# LARGE_DATASETS only picks the group a dataset's speedups are averaged in
# (geomeans.tsv).
#
# WINDOWS.  flickr and sx-stackoverflow run their first 1M edges, wikipedia
# its first 12M, the rest whole: the windows of 18678150, whose measurements
# place the crossings.  A larger window moves the crossing itself (locally
# coauthor's first 4M edges crossed near 900k, the whole stream near 320k on
# bryan8; flickr crossed earlier on its first 600k edges than on 1M).  Every
# batch size times the whole window: a dataset with an edge budget stops at
# its target (--max-edges), and a size that does not divide the window
# leaves a partial last batch, as the whole streams always had.  A budgeted
# size must stay below its window, so the window runs in two batches at
# least; checked before anything runs.
#
# WHERE THE GRID STOPS.  At its largest size every dataset runs two batches,
# youtube and coauthor three.  The crossings fall at about 2 batches on
# facebook, 3 on wikipedia and sx-stackoverflow, 5-6 on dblp-cite, flickr
# and youtube, and 26 on coauthor: on the smaller windows the sizes past the
# crossing come close to the static problem.
#
# WARM-UP (off by default): <DATASET>_WARMUP_EDGES > 0 inserts that many edges
# first without searching them (--warmup-search=off) and times only the window
# after it up to the target, the last batch partial as above.
#
# ARM ORDER.  Rotated one position from each cell to the next and from each
# repeat to the next.  REPEATS=1 by default (one run of every arm per cell);
# with more, the repeats are the outer loop.
#
# NODE.  Give the job its node to itself and do not run another timing job on
# the same node meanwhile: 18143576 and 18145683 overlapped on bryan1 for 3.6
# hours, and the first one's youtube cells came out about 7% slower for it.
# OzSTAR refuses --exclusive; ask instead for all of the node's cores and its
# memory spread over them (scontrol show node <node> gives CPUTot, and
# RealMemory in MB):
#   sbatch --nodelist=<node> --cpus-per-task=<CPUTot> \
#          --mem-per-cpu=<RealMemory/CPUTot, rounded down>M <this script>
# The searches are single-threaded either way; the extra cores only build the
# binary faster.  Report one node's times only: across the bryan nodes the
# ratios shifted by dataset (17899068).
#
# BUDGET.  From 18678150 on bryan8, through the model above: edge and vertex
# about 9.6 hours (wikipedia 4.6, youtube 2.8, flickr 1.4, coauthor and
# sx-stackoverflow 0.4 each, dblp-cite and facebook under a minute; reading
# the input is 10-50 s of every run on the large graphs), star about 1.5 more
# (its 50 cells there took 1.2 hours with their reading): about 11 hours.
# The longest vertex run, wikipedia b2000000, took about 18 minutes there,
# far from its hour.  The same binary took 1.28 times as long on bryan5 as on
# bryan1 (17859145 against 17899068), so allow 14.  DEADLINE_SECONDS stops
# cleanly inside the 24-hour allocation.
#
# BUILD.  Unless EXECUTABLE is given, the job builds its own binary from
# ${PROJECT_DIR}/src with the project Makefile, into a new directory under its
# output directory, and refuses a src/ without the vertex baseline (added
# 2026-10-09) or from before that day's cleanup.  Upload src/, Makefile and
# this script; load a C++20 compiler before sbatch or pass MODULES="gcc/13.3.0".
# Create the directory of the --output/--error paths above first: SLURM does
# not.

set -uo pipefail

SCRIPT_DIR="${PROJECT_DIR:-/home/zzhong/zzy_project/DynamicClique/IMMCE_SubsumedClique}"
PROJECT_DIR="${PROJECT_DIR:-$SCRIPT_DIR}"
DATA_PATH="${INPUT:-/fred/oz128/myDataset/DynamicGraph/Group1}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/edge_vertex_crossover_${SLURM_JOB_ID:-local}}"
EXECUTABLE_EXPLICIT=0
[[ -n "${EXECUTABLE:-}" ]] && EXECUTABLE_EXPLICIT=1
EXECUTABLE="${EXECUTABLE:-}"
BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-1}}"
MODULES="${MODULES:-}"

REPEATS="${REPEATS:-1}"
# The batch sizes of each dataset (BATCH SIZES above), in any order; they run
# largest first.  DEFAULT is for a dataset none of the others names.
DBLP_CITE_BATCH_SIZES="${DBLP_CITE_BATCH_SIZES:-3000 4000 5000 6000 8000 10000 12500 15000 20000 25000}"
FACEBOOK_BATCH_SIZES="${FACEBOOK_BATCH_SIZES:-75000 100000 125000 150000 175000 200000 250000 300000}"
COAUTHOR_BATCH_SIZES="${COAUTHOR_BATCH_SIZES:-50000 100000 200000 250000 300000 400000 500000 750000 1000000 2000000 4000000}"
FLICKR_BATCH_SIZES="${FLICKR_BATCH_SIZES:-75000 100000 125000 150000 175000 200000 250000 300000 400000 500000}"
WIKIPEDIA_BATCH_SIZES="${WIKIPEDIA_BATCH_SIZES:-2000000 2500000 3000000 3500000 4000000 4500000 5000000 6000000 8000000 10000000}"
YOUTUBE_BATCH_SIZES="${YOUTUBE_BATCH_SIZES:-750000 1000000 1250000 1500000 2000000 2500000 3000000 3500000 4500000}"
STACKOVERFLOW_BATCH_SIZES="${STACKOVERFLOW_BATCH_SIZES:-150000 200000 250000 300000 350000 400000 500000 600000 750000}"
DEFAULT_BATCH_SIZES="${DEFAULT_BATCH_SIZES:-1000 5000 10000 20000}"
# Set, it replaces every list above for every dataset.
BATCH_SIZES="${BATCH_SIZES:-}"
BATCH_LIST_NAMES="DBLP_CITE_BATCH_SIZES FACEBOOK_BATCH_SIZES COAUTHOR_BATCH_SIZES FLICKR_BATCH_SIZES WIKIPEDIA_BATCH_SIZES YOUTUBE_BATCH_SIZES STACKOVERFLOW_BATCH_SIZES DEFAULT_BATCH_SIZES"
# Only the group a dataset's speedups are averaged in (geomeans.tsv).
LARGE_DATASETS="${LARGE_DATASETS:-processed_dblp_coauthor.txt processed_flickr-growth.txt processed_wikipedia-growth.txt processed_youtube-u-growth.txt processed_sx-stackoverflow.txt}"
ARMS="${ARMS:-edge vertex star}"
# The arm every other one is compared with; it must be in ARMS.
METHOD_ARM="${METHOD_ARM:-star}"
# Space-separated name=from:to entries, each reported as time(from)/time(to).
# Empty means the comparisons described above, built from the arms given.
EXPERIMENTS="${EXPERIMENTS:-}"
TIE_ORDER="${TIE_ORDER:-file}"
# Dataset file names, run in this order; "all" means every processed_*.txt
# found.  A name that matches no file stops the job before anything runs.
DATASET_FILTER="${DATASET_FILTER:-processed_dblp-cite.txt processed_facebook-wosn-links.txt processed_dblp_coauthor.txt processed_flickr-growth.txt processed_wikipedia-growth.txt processed_youtube-u-growth.txt processed_sx-stackoverflow.txt}"
# Wall-time limits of one run: the vertex arm's, and every other arm's (0 =
# none for those; the vertex arm always has one).
VERTEX_TIMEOUT_SECONDS="${VERTEX_TIMEOUT_SECONDS:-3600}"
RUN_TIMEOUT_SECONDS="${RUN_TIMEOUT_SECONDS:-7200}"
SKIP_SMALLER_AFTER_TIMEOUT="${SKIP_SMALLER_AFTER_TIMEOUT:-1}"
# Leaves time for the epilogue (23.5 hours of the 24-hour allocation).  A run
# already started is not cut short by it, only by its own time limit.
DEADLINE_SECONDS="${DEADLINE_SECONDS:-84600}"

# Per-dataset edge budgets, the same as submit_star_branch_tree_times.sh.  0
# means the whole stream.  Every batch size times the same window (target
# minus warm-up); one that does not divide it leaves a partial last batch.
FLICKR_TARGET_EDGES="${FLICKR_TARGET_EDGES:-1000000}"
WIKIPEDIA_TARGET_EDGES="${WIKIPEDIA_TARGET_EDGES:-12000000}"
YOUTUBE_TARGET_EDGES="${YOUTUBE_TARGET_EDGES:-0}"
STACKOVERFLOW_TARGET_EDGES="${STACKOVERFLOW_TARGET_EDGES:-1000000}"
DEFAULT_TARGET_EDGES="${DEFAULT_TARGET_EDGES:-0}"
FLICKR_WARMUP_EDGES="${FLICKR_WARMUP_EDGES:-0}"
WIKIPEDIA_WARMUP_EDGES="${WIKIPEDIA_WARMUP_EDGES:-0}"
YOUTUBE_WARMUP_EDGES="${YOUTUBE_WARMUP_EDGES:-0}"
STACKOVERFLOW_WARMUP_EDGES="${STACKOVERFLOW_WARMUP_EDGES:-0}"
DEFAULT_WARMUP_EDGES="${DEFAULT_WARMUP_EDGES:-0}"

for list_name in DATASET_FILTER LARGE_DATASETS $BATCH_LIST_NAMES \
                 BATCH_SIZES ARMS EXPERIMENTS; do
    printf -v "$list_name" '%s' "${!list_name//,/ }"
done

# Under SLURM, messages from the checks below would otherwise land only in the
# .err file and leave the .out file empty, so a refusal prints on both.
note() {
    echo "$*" >&2
    echo "$*"
}
trap 'status=$?; (( status != 0 )) && echo "exited with status $status"' EXIT

echo "=== immce new-clique search: where the vertex baseline overtakes the edge one, time and branches ==="
echo "Job ID:         ${SLURM_JOB_ID:-local}"
echo "Host:           $(hostname)"
echo "Started:        $(date -Is)"
echo "Project dir:    $PROJECT_DIR"
echo "Data path:      $DATA_PATH"
echo "Dataset filter: $DATASET_FILTER"
echo "Output dir:     $OUTPUT_DIR"
echo "Executable:     ${EXECUTABLE:-(built by this job)}"
echo "Arms:           $ARMS"
echo "Method arm:     $METHOD_ARM"
echo "Tie order:      $TIE_ORDER"
if [[ -n "$BATCH_SIZES" ]]; then
    echo "Batch sizes:    $BATCH_SIZES (every dataset), largest first"
else
    echo "Batch sizes:    per dataset, largest first:"
    for list_name in $BATCH_LIST_NAMES; do
        printf '                %s=%s\n' "$list_name" "${!list_name}"
    done
fi
echo "Large datasets: $LARGE_DATASETS (their group in geomeans.tsv)"
echo "Repeats:        $REPEATS"
echo "Time limits:    vertex ${VERTEX_TIMEOUT_SECONDS}s, other arms ${RUN_TIMEOUT_SECONDS}s (skip smaller batches after a timeout: $SKIP_SMALLER_AFTER_TIMEOUT)"
echo "Deadline:       ${DEADLINE_SECONDS}s after start"
echo

target_edges_for() {
    case "$1" in
        processed_flickr-growth.txt) echo "$FLICKR_TARGET_EDGES" ;;
        processed_wikipedia-growth.txt) echo "$WIKIPEDIA_TARGET_EDGES" ;;
        processed_youtube-u-growth.txt) echo "$YOUTUBE_TARGET_EDGES" ;;
        processed_sx-stackoverflow.txt) echo "$STACKOVERFLOW_TARGET_EDGES" ;;
        *) echo "$DEFAULT_TARGET_EDGES" ;;
    esac
}

warmup_edges_for() {
    case "$1" in
        processed_flickr-growth.txt) echo "$FLICKR_WARMUP_EDGES" ;;
        processed_wikipedia-growth.txt) echo "$WIKIPEDIA_WARMUP_EDGES" ;;
        processed_youtube-u-growth.txt) echo "$YOUTUBE_WARMUP_EDGES" ;;
        processed_sx-stackoverflow.txt) echo "$STACKOVERFLOW_WARMUP_EDGES" ;;
        *) echo "$DEFAULT_WARMUP_EDGES" ;;
    esac
}

# large or small: whether a dataset (file name) is in LARGE_DATASETS.
dataset_class_for() {
    local large
    for large in $LARGE_DATASETS; do
        if [[ "$1" == "$large" ]]; then
            echo large
            return
        fi
    done
    echo small
}

# The batch sizes one dataset (file name) runs, largest first.
batch_sizes_for() {
    local sizes
    if [[ -n "$BATCH_SIZES" ]]; then
        sizes="$BATCH_SIZES"
    else
        case "$1" in
            processed_dblp-cite.txt) sizes="$DBLP_CITE_BATCH_SIZES" ;;
            processed_facebook-wosn-links.txt) sizes="$FACEBOOK_BATCH_SIZES" ;;
            processed_dblp_coauthor.txt) sizes="$COAUTHOR_BATCH_SIZES" ;;
            processed_flickr-growth.txt) sizes="$FLICKR_BATCH_SIZES" ;;
            processed_wikipedia-growth.txt) sizes="$WIKIPEDIA_BATCH_SIZES" ;;
            processed_youtube-u-growth.txt) sizes="$YOUTUBE_BATCH_SIZES" ;;
            processed_sx-stackoverflow.txt) sizes="$STACKOVERFLOW_BATCH_SIZES" ;;
            *) sizes="$DEFAULT_BATCH_SIZES" ;;
        esac
    fi
    tr ' ' '\n' <<< "$sizes" | grep -v '^$' | sort -rn -u | tr '\n' ' '
}

# Parses an arm name into ARM_FRAMEWORK (edge / vertex / star), ARM_SCOPE (the
# branch scope of a star arm, - otherwise) and ARM_LIMIT (its time limit in
# seconds, 0 for none).
set_arm() {
    ARM_SCOPE=-
    ARM_LIMIT="$RUN_TIMEOUT_SECONDS"
    case "$1" in
        edge) ARM_FRAMEWORK=edge ;;
        vertex) ARM_FRAMEWORK=vertex; ARM_LIMIT="$VERTEX_TIMEOUT_SECONDS" ;;
        star) ARM_FRAMEWORK=star; ARM_SCOPE=unbound ;;
        star-root) ARM_FRAMEWORK=star; ARM_SCOPE=root ;;
        *) return 1 ;;
    esac
}
ARM_FORM="edge, vertex, star or star-root"

if [[ ! "$REPEATS" =~ ^[1-9][0-9]*$ ]]; then
    note "Error: REPEATS must be a positive integer: $REPEATS"
    exit 2
fi
if [[ ! "$VERTEX_TIMEOUT_SECONDS" =~ ^[1-9][0-9]*$ ]]; then
    note "Error: VERTEX_TIMEOUT_SECONDS must be a positive integer (the vertex arm always runs under a limit): $VERTEX_TIMEOUT_SECONDS"
    exit 2
fi
if [[ "$SKIP_SMALLER_AFTER_TIMEOUT" != 0 && "$SKIP_SMALLER_AFTER_TIMEOUT" != 1 ]]; then
    note "Error: SKIP_SMALLER_AFTER_TIMEOUT must be 0 or 1: $SKIP_SMALLER_AFTER_TIMEOUT"
    exit 2
fi
for value in "$DEADLINE_SECONDS" "$RUN_TIMEOUT_SECONDS" \
             "$FLICKR_TARGET_EDGES" "$WIKIPEDIA_TARGET_EDGES" "$YOUTUBE_TARGET_EDGES" \
             "$STACKOVERFLOW_TARGET_EDGES" "$DEFAULT_TARGET_EDGES" \
             "$FLICKR_WARMUP_EDGES" "$WIKIPEDIA_WARMUP_EDGES" "$YOUTUBE_WARMUP_EDGES" \
             "$STACKOVERFLOW_WARMUP_EDGES" "$DEFAULT_WARMUP_EDGES"; do
    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        note "Error: deadlines, timeouts, edge targets and warm-ups must be nonnegative integers: $value"
        exit 2
    fi
done
# Edges sharing a timestamp keep the input file's order: the only order the
# reader has since the vertex and shuffle orders were removed (2026-10-08).
if [[ "$TIE_ORDER" != file ]]; then
    note "Error: TIE_ORDER must be file (the vertex and shuffle orders were removed): $TIE_ORDER"
    exit 2
fi
for list_name in $BATCH_LIST_NAMES; do
    if [[ -z "${!list_name// /}" ]]; then
        note "Error: $list_name must contain at least one batch size"
        exit 2
    fi
    for batch_size in ${!list_name}; do
        if [[ ! "$batch_size" =~ ^[1-9][0-9]*$ ]]; then
            note "Error: invalid batch size in $list_name: $batch_size"
            exit 2
        fi
    done
done
for batch_size in $BATCH_SIZES; do
    if [[ ! "$batch_size" =~ ^[1-9][0-9]*$ ]]; then
        note "Error: invalid batch size in BATCH_SIZES: $batch_size"
        exit 2
    fi
done

read -r -a arms <<< "$ARMS"
if (( ${#arms[@]} < 2 )); then
    note "Error: ARMS must name at least two arms to compare: $ARMS"
    exit 2
fi
declare -A arm_seen=()
for arm in "${arms[@]}"; do
    if ! set_arm "$arm"; then
        note "Error: unknown arm: $arm (arms are $ARM_FORM)"
        exit 2
    fi
    if [[ -n "${arm_seen[$arm]:-}" ]]; then
        note "Error: arm listed twice: $arm"
        exit 2
    fi
    arm_seen[$arm]=1
done
if [[ -z "${arm_seen[$METHOD_ARM]:-}" ]]; then
    note "Error: METHOD_ARM=$METHOD_ARM is not one of ARMS: $ARMS"
    exit 2
fi
if (( REPEATS > 1 && REPEATS % ${#arms[@]} != 0 )); then
    echo "Note: REPEATS=$REPEATS is not a multiple of the ${#arms[@]} arms, so inside one cell an arm does not take every position equally often."
fi

known_arm() {
    [[ -n "$1" && -n "${arm_seen[$1]:-}" ]]
}
# Adds name=from:to when both arms are run.
add_experiment() {
    if known_arm "$2" && known_arm "$3"; then
        EXPERIMENTS+=" $1=$2:$3"
    fi
}
# An arm name as it reads inside a comparison name.
arm_tag() {
    printf '%s' "${1//-/_}"
}
if [[ -z "$EXPERIMENTS" ]]; then
    for arm in star star-root; do
        known_arm "$arm" || continue
        add_experiment "edge_vs_$(arm_tag "$arm")" edge "$arm"
        add_experiment "vertex_vs_$(arm_tag "$arm")" vertex "$arm"
    done
    add_experiment vertex_vs_edge vertex edge
    add_experiment star_root_vs_star star-root star
    # The method against any arm the lists above did not reach.
    for arm in "${arms[@]}"; do
        [[ "$arm" == "$METHOD_ARM" ]] && continue
        case " $EXPERIMENTS " in
            *"=$arm:$METHOD_ARM "*) continue ;;
        esac
        add_experiment "$(arm_tag "$arm")_vs_$(arm_tag "$METHOD_ARM")" "$arm" "$METHOD_ARM"
    done
    EXPERIMENTS="${EXPERIMENTS# }"
fi
if [[ -z "$EXPERIMENTS" ]]; then
    note "Error: the arms give no comparison: $ARMS"
    exit 2
fi
for experiment in $EXPERIMENTS; do
    name="${experiment%%=*}"
    IFS=: read -r from to extra <<< "${experiment#*=}"
    if [[ "$experiment" != *=* || -z "$name" || -z "$from" || -z "$to" || -n "$extra" ]]; then
        note "Error: experiment '$experiment' is not name=from:to"
        exit 2
    fi
    if ! known_arm "$from" || ! known_arm "$to"; then
        note "Error: experiment '$experiment' names an arm that is not in ARMS"
        exit 2
    fi
done
echo "Experiments:    $EXPERIMENTS"
echo

mkdir -p "$OUTPUT_DIR" || { note "Error: cannot create $OUTPUT_DIR"; exit 2; }
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

if [[ -n "$MODULES" ]]; then
    if type module > /dev/null 2>&1; then
        # shellcheck disable=SC2086
        if ! module load $MODULES; then
            note "Error: module load $MODULES failed"
            exit 2
        fi
    else
        note "Warning: MODULES is set but this shell has no module command; using the inherited environment."
    fi
fi

if ! command -v timeout > /dev/null 2>&1; then
    note "Error: no timeout on PATH (coreutils); the vertex arm needs it for its time limit."
    exit 2
fi

# ------------------------------------------------------------------ build

BUILD_CXXFLAGS="-std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow"
BUILD_DIR=""
if (( EXECUTABLE_EXPLICIT == 0 )); then
    if [[ ! -f "${PROJECT_DIR}/Makefile" || ! -f "${PROJECT_DIR}/src/star_clique_search.cpp" ]]; then
        note "Error: ${PROJECT_DIR} has no Makefile or no src/; upload both, or set PROJECT_DIR."
        exit 2
    fi
    # Each marker is a line of the current src/ that an older upload lacks.
    for marker in \
        "src/main.cpp|tie-order=file|its reader has no --tie-order=file (2026-10-05)" \
        "src/vertex_clique_search.hpp|class VertexCliqueEnumerator|it has no vertex baseline (added 2026-10-09)" \
        "src/main.cpp|FrameworkOption::Vertex|its driver cannot run the vertex baseline (added 2026-10-09)" \
        "src/main.cpp|\"edge: time-ms=\"|it is from before the 2026-10-09 cleanup (the edge line was still edge-fast:)" \
        "src/star_clique_search.hpp|StarBranchScope branchScope = StarBranchScope::Unbound;|it is from before the 2026-10-09 cleanup (star still had its root shapes and partitions as switches)" \
        "src/star_clique_search.hpp|std::deque<LocalVertexSet> deepMasks_;|its deep re-indexing keeps nested partner masks in a vector, which loses cliques" \
        "src/star_clique_search.hpp|branchUnitBindings|its unbound scope walks the star's whole partner list and scans a pivot for a lone partner (before 2026-10-08)" \
        "src/local_search_graph.hpp|moveToExcludedBits|its BK keeps sorted vertex lists on every sibling move" \
        "src/local_bk.hpp|candidateCount - 1|its pivot scan has no early exit on the P side" \
        "src/subsumed_clique_search.cpp|if (readMemberLists(candidate) >= static_cast<std::size_t>(deferDegree)) {|it is older than the tree of the full-grid jobs (2026-10-14), whose binary this job is meant to share"; do
        IFS='|' read -r marker_file marker_text marker_reason <<< "$marker"
        if ! grep -qF -- "$marker_text" "${PROJECT_DIR}/${marker_file}" 2> /dev/null; then
            note "Error: ${PROJECT_DIR}/src is older than this script: $marker_reason."
            note "Upload the current src/ and Makefile."
            exit 2
        fi
    done
    if ! grep -qF 'src/vertex_clique_search.cpp' "${PROJECT_DIR}/Makefile"; then
        note "Error: ${PROJECT_DIR}/Makefile does not build src/vertex_clique_search.cpp; upload the current Makefile."
        exit 2
    fi
    for tool in make g++; do
        if ! command -v "$tool" > /dev/null 2>&1; then
            note "Error: no $tool on PATH, so the binary cannot be built."
            note "Load a compiler first (e.g. module load gcc/13.3.0), pass"
            note "MODULES=gcc/13.3.0, or build elsewhere and pass EXECUTABLE."
            exit 2
        fi
    done
    gcc_major="$(g++ -dumpversion 2> /dev/null | cut -d. -f1)"
    if [[ "$gcc_major" =~ ^[0-9]+$ ]] && (( gcc_major < 11 )); then
        note "Warning: g++ $gcc_major may be too old for this C++20 tree (the local builds use 13.3)."
    fi
    # A new directory every time: make decides what to rebuild by timestamps,
    # and an upload that keeps the local ones can leave newer stale objects in
    # an old build directory.
    BUILD_DIR="${OUTPUT_DIR}/build_$(date +%Y%m%d_%H%M%S)"
    echo "[$(date -Is)] building into $BUILD_DIR with: $BUILD_CXXFLAGS"
    if ! make -C "$PROJECT_DIR" -j"$BUILD_JOBS" BUILD_DIR="$BUILD_DIR" \
            CXXFLAGS="$BUILD_CXXFLAGS" all > "${OUTPUT_DIR}/build.log" 2>&1; then
        note "Error: build failed; see ${OUTPUT_DIR}/build.log"
        tail -n 30 "${OUTPUT_DIR}/build.log" >&2
        exit 2
    fi
    if grep -q 'warning:' "${OUTPUT_DIR}/build.log"; then
        note "Warning: the build reported compiler warnings; see ${OUTPUT_DIR}/build.log"
    fi
    EXECUTABLE="${BUILD_DIR}/immce_star"
    echo "Built $EXECUTABLE"
fi

if [[ ! -x "$EXECUTABLE" ]]; then
    note "Error: executable not found or not executable: $EXECUTABLE"
    exit 2
fi

usage_text="$("$EXECUTABLE" --help 2>&1 || true)"
for required in \
    framework star-branch-scope tie-order warmup-search maintenance \
    max-batches max-edges warmup-edges warmup-batch-size no-output stats \
    'edge|vertex'; do
    if [[ "$usage_text" != *"$required"* ]]; then
        note "Error: executable lacks required option: $required"
        note "It is not a build of the current IMMCE(SubCliqueSearch) tree."
        exit 2
    fi
done

# ---------------------------------------------------------------- datasets

discover_datasets() {
    local data_path="$1"

    if [[ -f "$data_path" ]]; then
        DATASETS=("$data_path")
        return 0
    fi
    if [[ ! -d "$data_path" ]]; then
        note "Error: dataset path is neither a file nor a directory: $data_path"
        return 2
    fi
    mapfile -d '' -t DATASETS < <(
        find "$data_path" -type f -name 'processed_*.txt' -print0 | sort -z
    )
    if (( ${#DATASETS[@]} == 0 )); then
        note "Error: no processed_*.txt datasets found under: $data_path"
        return 2
    fi
}

DATASETS=()
if ! discover_datasets "$DATA_PATH"; then
    exit 2
fi
# Every dataset the lists name has to exist: a misspelt or misjoined name
# would otherwise quietly drop a dataset from an hours-long job, or run a
# large one at the small batch sizes.
declare -A dataset_present=()
for dataset in "${DATASETS[@]}"; do
    dataset_present["$(basename -- "$dataset")"]=1
done
missing=""
if [[ "$DATASET_FILTER" != all ]]; then
    for wanted in $DATASET_FILTER; do
        [[ -n "${dataset_present[$wanted]:-}" ]] || missing="${missing:+$missing }$wanted"
    done
fi
if [[ -d "$DATA_PATH" ]]; then
    for wanted in $LARGE_DATASETS; do
        [[ -n "${dataset_present[$wanted]:-}" ]] || missing="${missing:+$missing }$wanted"
    done
fi
if [[ -n "$missing" ]]; then
    note "Error: DATASET_FILTER or LARGE_DATASETS names files that are not under $DATA_PATH: $missing"
    note "Lists take dataset file names separated by spaces or commas."
    exit 2
fi
if [[ "$DATASET_FILTER" != all && -n "$DATASET_FILTER" ]]; then
    kept=()
    for wanted in $DATASET_FILTER; do
        for dataset in "${DATASETS[@]}"; do
            if [[ "$(basename -- "$dataset")" == "$wanted" ]]; then
                kept+=("$dataset")
                break
            fi
        done
    done
    DATASETS=("${kept[@]}")
fi
# A dataset with an edge budget times its window (target minus warm-up) at
# every batch size, the last batch partial when the size does not divide it.
# A size of the whole window or more would run it as one batch, the static
# problem: refused here rather than run hours in.
window_errors=""
for dataset in "${DATASETS[@]}"; do
    name="$(basename -- "$dataset")"
    target="$(target_edges_for "$name")"
    warmup="$(warmup_edges_for "$name")"
    (( target > 0 && target > warmup )) || continue
    window=$(( target - warmup ))
    for batch_size in $(batch_sizes_for "$name"); do
        if (( batch_size >= window )); then
            window_errors="${window_errors:+$window_errors }$name:$batch_size"
        fi
    done
done
if [[ -n "$window_errors" ]]; then
    note "Error: these batch sizes would run their dataset's window (target minus warm-up) as a single batch: $window_errors"
    exit 2
fi

LAUNCHER=()
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    LAUNCHER=(srun --exclusive)
fi

# The time limit of one arm, inside the launcher, so that the job step itself
# is what the limit stops; empty for an arm without one.
limit_command_for() {
    LIMIT_COMMAND=()
    if (( ARM_LIMIT > 0 )); then
        LIMIT_COMMAND=(timeout --signal=TERM --kill-after=60 "$ARM_LIMIT")
    fi
}

extract_stat() {
    local label="$1"
    local key="$2"
    local stats_file="$3"
    sed -nE "s/^${label}:.*[, ]${key}=([^, ]+).*$/\\1/p" "$stats_file" | tail -n1
}

dataset_tag_of() {
    local dataset="$1"
    local name
    if [[ -d "$DATA_PATH" && "$dataset" == "$DATA_PATH"/* ]]; then
        name="${dataset#"$DATA_PATH"/}"
    else
        name="$(basename -- "$dataset")"
    fi
    name="${name//\//__}"
    printf '%s' "${name// /_}"
}

smallest_dataset() {
    local smallest="" size best=0
    for dataset in "${DATASETS[@]}"; do
        size="$(stat -c %s -- "$dataset" 2>/dev/null || echo 0)"
        if [[ -z "$smallest" ]] || (( size < best )); then
            smallest="$dataset"
            best="$size"
        fi
    done
    printf '%s' "$smallest"
}

# The command line of one arm, into ARM_ARGUMENTS.
arm_arguments() {
    set_arm "$1"
    ARM_ARGUMENTS=(
        "--framework=${ARM_FRAMEWORK}"
        --maintenance=new-only
        "--tie-order=${TIE_ORDER}"
    )
    if [[ "$ARM_FRAMEWORK" == star ]]; then
        ARM_ARGUMENTS+=("--star-branch-scope=${ARM_SCOPE}")
    fi
}

# Killed by its limit (or past it): a timeout.  Out of memory: std::bad_alloc,
# SLURM's oom-kill message, or a SIGKILL (137) before the limit.  Neither is a
# failure of the job.
run_status() {
    local exit_code="$1" wall="$2" stats_file="$3" limit="$4"
    if (( exit_code == 0 )); then
        echo ok
    elif (( limit > 0 )) && { (( exit_code == 124 )) || \
            awk -v w="$wall" -v l="$limit" 'BEGIN {exit !(w >= l)}'; }; then
        echo timeout
    elif (( exit_code == 137 )) || \
         grep -qiE 'bad_alloc|oom-kill|oom_kill|out of memory' "$stats_file" 2> /dev/null; then
        echo out_of_memory
    else
        echo failed
    fi
}

# One batch per arm proves the parser takes every switch and a run completes
# before hours are spent on the sweep.
probe_binary() {
    local probe_log="${OUTPUT_DIR}/probe.log"
    local smallest
    smallest="$(smallest_dataset)"
    : > "$probe_log"
    echo "[$(date -Is)] probing $EXECUTABLE on $(dataset_tag_of "$smallest")"
    for arm in "${arms[@]}"; do
        arm_arguments "$arm"
        timeout --kill-after=60 600 "$EXECUTABLE" "$smallest" "${ARM_ARGUMENTS[@]}" \
            --batch-size=1000 --max-batches=1 --no-output --stats \
            >> "$probe_log" 2>&1 || return 1
        grep -q "^${ARM_FRAMEWORK}: time-ms=" "$probe_log" || return 1
    done
    grep -q "tie-order=${TIE_ORDER}" "$probe_log" || return 1
    return 0
}

if ! probe_binary; then
    note "Error: $EXECUTABLE did not complete a probe run; see ${OUTPUT_DIR}/probe.log"
    exit 2
fi
echo "Probe runs completed."

JOB_START="$(date +%s)"
SUMMARY_FILE="${OUTPUT_DIR}/summary.tsv"
RESULTS_FILE="${OUTPUT_DIR}/results.tsv"
TIMES_FILE="${OUTPUT_DIR}/times.tsv"
TREE_FILE="${OUTPUT_DIR}/tree.tsv"
TABLE_TIME_FILE="${OUTPUT_DIR}/table_time.tsv"
TABLE_TREE_FILE="${OUTPUT_DIR}/table_tree.tsv"
GEOMEAN_FILE="${OUTPUT_DIR}/geomeans.tsv"
NOISE_FILE="${OUTPUT_DIR}/noise.tsv"
TIMEOUT_FILE="${OUTPUT_DIR}/timeouts.tsv"
CROSSOVER_FILE="${OUTPUT_DIR}/crossover.tsv"
CROSSOVER_SUMMARY_FILE="${OUTPUT_DIR}/crossover_summary.tsv"
# The header and the row are written from one list and one array, and
# write_summary_row refuses a row whose length disagrees with it.
SUMMARY_COLUMNS=(
    dataset dataset_class batch_size warmup_edges max_batches target_edges
    arm framework branch_scope tie_order position repeat
    status exit_code wall_seconds time_limit_seconds
    new_clique_ms graph_update_ms input_preprocess_ms warmup_ms
    clique_counts_match switches_check
    processed_batches processed_edges new_maximal_cliques
    branches search_roots branches_below_roots leaves emitted pivot_calls pivot_candidates
    local_builds local_vertices max_local_vertices peak_local_bytes
    edge_seeds edge_ownership_moves edge_snapshot_lookups
    vertex_roots vertex_pruned_roots vertex_new_edge_prunes vertex_batch_edge_checks
    star_stars star_singleton_stars star_per_partner_star_edges star_pivot_root_stars
    star_partner_prunes star_partnerless_nodes star_node_words
    star_deep_reindexed_subtrees star_ownership_probes
    stats_file
)
(IFS=$'\t'; printf '%s\n' "${SUMMARY_COLUMNS[*]}") > "$SUMMARY_FILE"
printf 'dataset\tbatch_size\tarm\trepeat\ttimeout_seconds\tstats_file\n' > "$TIMEOUT_FILE"

{
    echo "experiment=immce-subcliquesearch-edge-vertex-crossover-time-and-branches"
    echo "window_last_batch=partial when the batch size does not divide the window (--max-edges=target)"
    echo "dynamic_semantics=insertion-only"
    echo "project_dir=$PROJECT_DIR"
    echo "data_path=$DATA_PATH"
    echo "dataset_filter=$DATASET_FILTER"
    echo "executable=$EXECUTABLE"
    echo "executable_sha256=$(sha256sum "$EXECUTABLE" | cut -d' ' -f1)"
    echo "built_by_job=$(( EXECUTABLE_EXPLICIT == 0 ? 1 : 0 ))"
    echo "build_dir=${BUILD_DIR:-(prebuilt)}"
    echo "build_cxxflags=$BUILD_CXXFLAGS"
    echo "compiler=$(g++ --version 2> /dev/null | head -n1)"
    echo "timing_scope=new-only count-only: star: / edge: / vertex: time-ms, no clique materialised"
    echo "branches=recursive-calls, roots included; search_roots: edge seeds, vertex roots - pruned-roots, star singleton-stars + per-partner-star-edges + branch-root-pivot-stars"
    echo "tie_order=$TIE_ORDER"
    echo "method_arm=$METHOD_ARM"
    echo "dataset_count=${#DATASETS[@]}"
    echo "arms=$ARMS"
    echo "experiments=$EXPERIMENTS"
    echo "batch_sizes_override=${BATCH_SIZES:-(none)}"
    for list_name in $BATCH_LIST_NAMES; do
        echo "$(tr '[:upper:]' '[:lower:]' <<< "$list_name")=${!list_name}"
    done
    for dataset in "${DATASETS[@]}"; do
        sizes="$(batch_sizes_for "$(basename -- "$dataset")")"
        echo "batch_sizes.$(basename -- "$dataset")=${sizes% }"
    done
    echo "large_datasets=$LARGE_DATASETS (their group in geomeans.tsv)"
    echo "batch_order=largest first per dataset"
    echo "repeats=$REPEATS"
    echo "arm_order=rotated one position per cell and per repeat"
    echo "vertex_timeout_seconds=$VERTEX_TIMEOUT_SECONDS"
    echo "run_timeout_seconds=$RUN_TIMEOUT_SECONDS (every other arm; 0 = none)"
    echo "skip_smaller_after_timeout=$SKIP_SMALLER_AFTER_TIMEOUT"
    echo "timeout_bound=(the run's limit - largest wall-minus-search time of a finished run in the cell) / time(to)"
    echo "flickr_target_edges=$FLICKR_TARGET_EDGES"
    echo "wikipedia_target_edges=$WIKIPEDIA_TARGET_EDGES"
    echo "youtube_target_edges=$YOUTUBE_TARGET_EDGES"
    echo "stackoverflow_target_edges=$STACKOVERFLOW_TARGET_EDGES"
    echo "default_target_edges=$DEFAULT_TARGET_EDGES"
    echo "flickr_warmup_edges=$FLICKR_WARMUP_EDGES"
    echo "wikipedia_warmup_edges=$WIKIPEDIA_WARMUP_EDGES"
    echo "youtube_warmup_edges=$YOUTUBE_WARMUP_EDGES"
    echo "stackoverflow_warmup_edges=$STACKOVERFLOW_WARMUP_EDGES"
    echo "default_warmup_edges=$DEFAULT_WARMUP_EDGES"
    echo "deadline_seconds=$DEADLINE_SECONDS"
    echo "started=$(date -Is)"
    echo "slurm_job_id=${SLURM_JOB_ID:-local}"
    echo "host=$(hostname)"
    for arm in "${arms[@]}"; do
        arm_arguments "$arm"
        echo "arm_spec.$arm=${ARM_ARGUMENTS[*]} (time limit ${ARM_LIMIT}s)"
    done
} > "${OUTPUT_DIR}/manifest.txt"
if (( EXECUTABLE_EXPLICIT == 0 )); then
    (cd "$PROJECT_DIR" && sha256sum Makefile src/*.cpp src/*.hpp) \
        > "${OUTPUT_DIR}/source_sha256.txt"
fi
printf '%s\n' "${DATASETS[@]}" > "${OUTPUT_DIR}/datasets.txt"

deadline_passed() {
    local now
    now="$(date +%s)"
    (( now - JOB_START >= DEADLINE_SECONDS ))
}

failures=0
timeouts=0
out_of_memory_runs=0
run_count=0
skipped=0
# Set by the first arm of a cell that finishes; every later arm must match it.
cell_reference=""
# dataset tag|arm -> the largest batch size at which that arm timed out.
declare -A timed_out_at=()

is_number() {
    [[ "$1" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]
}

# One summary row, built as an array and joined on tabs.
write_summary_row() {
    local -a row=("$@")
    local index
    for index in "${!row[@]}"; do
        [[ -z "${row[$index]}" ]] && row[$index]='-'
    done
    if (( ${#row[@]} != ${#SUMMARY_COLUMNS[@]} )); then
        note "Error: the summary row has ${#row[@]} fields but the header has ${#SUMMARY_COLUMNS[@]}."
        note "Refusing to write a misaligned row."
        exit 2
    fi
    (IFS=$'\t'; printf '%s\n' "${row[*]}") >> "$SUMMARY_FILE"
}

# The roots of one run's search trees, from the line of its own method; empty
# if a field is missing.
search_roots_of() {
    local framework="$1" stats_file="$2"
    local a b c
    case "$framework" in
        edge)
            extract_stat edge seeds "$stats_file"
            ;;
        vertex)
            a="$(extract_stat vertex roots "$stats_file")"
            b="$(extract_stat vertex pruned-roots "$stats_file")"
            if [[ "$a" =~ ^[0-9]+$ && "$b" =~ ^[0-9]+$ ]]; then
                echo $(( a - b ))
            fi
            ;;
        star)
            a="$(extract_stat star singleton-stars "$stats_file")"
            b="$(extract_stat star per-partner-star-edges "$stats_file")"
            c="$(extract_stat star branch-root-pivot-stars "$stats_file")"
            if [[ "$a" =~ ^[0-9]+$ && "$b" =~ ^[0-9]+$ && "$c" =~ ^[0-9]+$ ]]; then
                echo $(( a + b + c ))
            fi
            ;;
    esac
}

# One measured run: one dataset, one batch size, one arm, one repeat.
run_cell() {
    local dataset="$1" batch_size="$2" max_batches="$3" target="$4"
    local arm="$5" repeat="$6" position="$7" warmup="$8"
    local tag stats_file start end wall exit_code status index class
    local extra_arguments=()

    tag="$(dataset_tag_of "$dataset")"
    class="$(dataset_class_for "$(basename -- "$dataset")")"
    arm_arguments "$arm"
    limit_command_for
    local -a identity=(
        "$tag" "$class" "$batch_size" "$warmup" "${max_batches:-0}" "$target"
        "$arm" "$ARM_FRAMEWORK" "$ARM_SCOPE" "$TIE_ORDER"
        "$position" "$repeat"
    )

    # An arm that already timed out here at this or a larger batch size would
    # only time out again (the vertex baseline slows down as batches shrink).
    if (( SKIP_SMALLER_AFTER_TIMEOUT == 1 )) && \
       [[ -n "${timed_out_at[$tag|$arm]:-}" ]] && \
       (( batch_size <= ${timed_out_at[$tag|$arm]} )); then
        echo "[$(date -Is)] skip dataset=$tag bs=$batch_size arm=$arm (timed out at bs=${timed_out_at[$tag|$arm]})"
        local -a blank=()
        local count=$(( ${#SUMMARY_COLUMNS[@]} - ${#identity[@]} - 1 ))
        for ((index = 0; index < count; ++index)); do blank+=('-'); done
        blank[0]=skipped_after_timeout
        blank[3]="$ARM_LIMIT"
        write_summary_row "${identity[@]}" "${blank[@]}" '-'
        ((skipped += 1))
        return
    fi

    stats_file="${OUTPUT_DIR}/${tag}.b${batch_size}.${arm}.r${repeat}.log"
    if (( max_batches > 0 )); then
        # The stream stops at the target, so a partial last batch ends the
        # window exactly.
        extra_arguments+=("--max-batches=${max_batches}" "--max-edges=${target}")
    fi
    if (( warmup > 0 )); then
        # One insertion, no search: in new-only mode the warm-up's searches
        # leave nothing behind but output.
        extra_arguments+=("--warmup-edges=${warmup}" "--warmup-batch-size=${warmup}"
                          --warmup-search=off)
    fi

    echo "[$(date -Is)] run dataset=$tag bs=$batch_size warmup=$warmup arm=$arm repeat=$repeat position=$position"
    start="$(date +%s.%N)"
    "${LAUNCHER[@]}" "${LIMIT_COMMAND[@]}" "$EXECUTABLE" "$dataset" \
        "${ARM_ARGUMENTS[@]}" \
        "--batch-size=${batch_size}" \
        "${extra_arguments[@]}" \
        --no-output \
        --stats \
        > "$stats_file" 2>&1
    exit_code=$?
    end="$(date +%s.%N)"
    wall="$(awk -v s="$start" -v e="$end" 'BEGIN {printf "%.3f", e - s}')"
    status="$(run_status "$exit_code" "$wall" "$stats_file" "$ARM_LIMIT")"

    local label="$ARM_FRAMEWORK"
    local new_ms new_cliques match='-'
    new_ms="$(extract_stat "$label" time-ms "$stats_file")"
    new_cliques="$(extract_stat graph new-maximal-cliques "$stats_file")"
    # Every arm is exact, so a cell's new-clique count must not depend on it.
    if [[ "$status" == ok ]]; then
        if [[ -z "$cell_reference" ]]; then
            cell_reference="${new_cliques:--}"
            match='reference'
        elif [[ "${new_cliques:--}" == "$cell_reference" ]]; then
            match='yes'
        else
            match='no'
        fi
    fi

    # The arm must have run what its name says, on the count-only path.
    local switches='-'
    if [[ "$status" == ok ]]; then
        switches='ok'
        local other
        if [[ "$(extract_stat graph framework "$stats_file")" != "$ARM_FRAMEWORK" ]]; then
            switches='wrong-framework'
        elif [[ "$(extract_stat graph tie-order "$stats_file")" != "$TIE_ORDER" ]]; then
            switches='wrong-tie-order'
        elif [[ "$(extract_stat graph maintenance "$stats_file")" != new-only || \
                "$(extract_stat graph result-mode "$stats_file")" != count-only ]]; then
            switches='not-count-only'
        elif [[ "$ARM_FRAMEWORK" == star && \
                "$(extract_stat star branch-scope "$stats_file")" != "$ARM_SCOPE" ]]; then
            switches='wrong-branch-scope'
        elif (( warmup > 0 )) && \
             [[ "$(extract_stat graph warmup-search "$stats_file")" != off ]]; then
            switches='warmup-searched'
        else
            for other in edge vertex star; do
                [[ "$other" == "$ARM_FRAMEWORK" ]] && continue
                if [[ -n "$(extract_stat "$other" time-ms "$stats_file")" ]]; then
                    switches='two-producers-ran'
                fi
            done
        fi
    fi

    local branches roots below='-'
    branches="$(extract_stat "$label" recursive-calls "$stats_file")"
    roots="$(search_roots_of "$ARM_FRAMEWORK" "$stats_file")"
    if [[ "$branches" =~ ^[0-9]+$ && "$roots" =~ ^[0-9]+$ ]]; then
        below=$(( branches - roots ))
    fi
    # The fields of one method's own line, empty for the other methods.
    vertex_field() { [[ "$ARM_FRAMEWORK" == vertex ]] && extract_stat vertex "$1" "$stats_file"; }
    star_field() { [[ "$ARM_FRAMEWORK" == star ]] && extract_stat star "$1" "$stats_file"; }

    write_summary_row "${identity[@]}" \
        "$status" "$exit_code" "$wall" "$ARM_LIMIT" \
        "${new_ms:-}" \
        "$(extract_stat timing graph-update-ms "$stats_file")" \
        "$(extract_stat timing input-preprocess-ms "$stats_file")" \
        "$(extract_stat graph warmup-ms "$stats_file")" \
        "$match" "$switches" \
        "$(extract_stat graph batches "$stats_file")" \
        "$(extract_stat graph processed-edges "$stats_file")" \
        "${new_cliques:-}" \
        "${branches:-}" "${roots:-}" "$below" \
        "$(extract_stat "$label" leaves "$stats_file")" \
        "$(extract_stat "$label" emitted "$stats_file")" \
        "$(extract_stat "$label" pivot-calls "$stats_file")" \
        "$(extract_stat "$label" pivot-candidates "$stats_file")" \
        "$(extract_stat "$label" local-builds "$stats_file")" \
        "$(extract_stat "$label" local-vertices "$stats_file")" \
        "$(extract_stat "$label" max-local-vertices "$stats_file")" \
        "$(extract_stat "$label" peak-local-bytes "$stats_file")" \
        "$([[ "$ARM_FRAMEWORK" == edge ]] && extract_stat edge seeds "$stats_file")" \
        "$([[ "$ARM_FRAMEWORK" == edge ]] && extract_stat edge ownership-moves "$stats_file")" \
        "$([[ "$ARM_FRAMEWORK" == edge ]] && extract_stat edge snapshot-edge-lookups "$stats_file")" \
        "$(vertex_field roots)" \
        "$(vertex_field pruned-roots)" \
        "$(vertex_field new-edge-prunes)" \
        "$(vertex_field batch-edge-checks)" \
        "$(star_field stars)" \
        "$(star_field singleton-stars)" \
        "$(star_field per-partner-star-edges)" \
        "$(star_field branch-root-pivot-stars)" \
        "$(star_field partner-prunes)" \
        "$(star_field partnerless-nodes)" \
        "$(star_field node-words)" \
        "$(star_field deep-reindexed-subtrees)" \
        "$(star_field ownership-probes)" \
        "$stats_file"

    case "$status" in
        timeout)
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$tag" "$batch_size" "$arm" "$repeat" \
                "$ARM_LIMIT" "$stats_file" >> "$TIMEOUT_FILE"
            note "Note: arm=$arm hit its ${ARM_LIMIT}s limit on $tag bs=$batch_size; recorded in timeouts.tsv"
            ((timeouts += 1))
            if [[ -z "${timed_out_at[$tag|$arm]:-}" ]] || \
               (( batch_size > ${timed_out_at[$tag|$arm]} )); then
                timed_out_at[$tag|$arm]="$batch_size"
            fi
            ;;
        out_of_memory)
            note "Note: arm=$arm ran out of memory on $tag bs=$batch_size (exit code $exit_code); see $stats_file"
            ((out_of_memory_runs += 1))
            ;;
        failed)
            note "Warning: arm=$arm failed with exit code $exit_code on $tag bs=$batch_size; see $stats_file"
            ((failures += 1))
            ;;
        ok)
            if [[ "$match" == no ]]; then
                note "Error: arm=$arm changed the new-clique count on $tag bs=$batch_size:"
                echo "  ${new_cliques:--} vs $cell_reference (see $stats_file)" >&2
                ((failures += 1))
            elif [[ "$switches" != ok ]]; then
                note "Error: arm=$arm did not do what its name says ($switches) on $tag bs=$batch_size; see $stats_file"
                ((failures += 1))
            elif ! is_number "${new_ms:-}"; then
                note "Warning: the new-clique time is missing from $stats_file"
                ((failures += 1))
            elif [[ ! "${branches:-}" =~ ^[0-9]+$ || "$below" == - ]]; then
                note "Warning: the search statistics are missing from $stats_file"
                ((failures += 1))
            fi
            ;;
    esac
    ((run_count += 1))
}

# Repeats are the outer loop so that stopping early costs precision rather than
# whole datasets; inside a dataset the batch sizes run largest first; the arms
# of one cell run back to back, so the pairing inside a cell is as tight as it
# can be, their order rotated one position per cell and per repeat.
arm_count=${#arms[@]}
for ((repeat = 1; repeat <= REPEATS; ++repeat)); do
    cell_number=0
    for dataset in "${DATASETS[@]}"; do
        name="$(basename -- "$dataset")"
        target="$(target_edges_for "$name")"
        warmup="$(warmup_edges_for "$name")"
        if (( target > 0 && warmup >= target )); then
            if (( repeat == 1 )); then
                note "Warning: $(dataset_tag_of "$dataset"): the $warmup-edge warm-up leaves no window below the $target-edge target; skipped"
            fi
            ((skipped += 1))
            continue
        fi
        for batch_size in $(batch_sizes_for "$name"); do
            max_batches=0
            if (( target > 0 )); then
                window=$(( target - warmup ))
                max_batches=$(( (window + batch_size - 1) / batch_size ))
            fi
            shift_by=$(( (repeat - 1 + cell_number) % arm_count ))
            ((cell_number += 1))
            order=()
            for ((index = 0; index < arm_count; ++index)); do
                order+=("${arms[(index + shift_by) % arm_count]}")
            done
            cell_reference=""
            position=0
            for arm in "${order[@]}"; do
                if deadline_passed; then
                    echo "Deadline reached; stopping after $run_count runs." >&2
                    ((skipped += 1))
                    break 4
                fi
                run_cell "$dataset" "$batch_size" "$max_batches" \
                    "$target" "$arm" "$repeat" "$position" "$warmup"
                ((position += 1))
            done
        done
    done
done

# Ratios are formed inside one cell and one repeat: time(from)/time(to) when
# both finished, its lower bound when `from` hit its limit and `to` finished.
# A cell's speedup is the median of its finished repeats, or of its bounds
# when none finished.  Runs that failed, disagreed on the count or did not run
# what their arm says enter nothing.
awk -F'\t' -v OFS='\t' \
    -v arm_list="$ARMS" \
    -v experiment_list="$EXPERIMENTS" \
    -v repeats="$REPEATS" \
    -v method_arm="$METHOD_ARM" \
    -v results_file="$RESULTS_FILE" \
    -v times_file="$TIMES_FILE" \
    -v tree_file="$TREE_FILE" \
    -v table_time_file="$TABLE_TIME_FILE" \
    -v table_tree_file="$TABLE_TREE_FILE" \
    -v geomean_file="$GEOMEAN_FILE" \
    -v noise_file="$NOISE_FILE" \
    -v crossover_file="$CROSSOVER_FILE" \
    -v crossover_summary_file="$CROSSOVER_SUMMARY_FILE" '
function median(list,   values, n, i, j, key) {
    n = split(list, values, " ")
    if (n == 0) { return "" }
    for (i = 2; i <= n; ++i) {
        key = values[i] + 0
        for (j = i - 1; j >= 1 && values[j] + 0 > key; --j) {
            values[j + 1] = values[j]
        }
        values[j + 1] = key
    }
    if (n % 2 == 1) { return values[(n + 1) / 2] + 0 }
    return (values[n / 2] + values[n / 2 + 1]) / 2
}
function fmt_ms(v) { return (v == "") ? "-" : sprintf("%.3f", v) }
function fmt_ratio(v) { return (v == "") ? "-" : sprintf("%.4f", v) }
function number(v) { return v ~ /^[0-9]+([.][0-9]*)?([eE][-+]?[0-9]+)?$/ }
function count(v) { return v ~ /^[0-9]+$/ }
# The median time of one arm in one cell over its finished repeats, "" if none.
function arm_median(d, b, a,   here, n, i, r, list) {
    list = ""
    n = split(repeat_list[d, b], here, " ")
    for (i = 1; i <= n; ++i) {
        r = here[i]
        if (state[d, b, r, a] == "ok") { list = list " " value[d, b, r, a] }
    }
    return median(list)
}
# What one arm did in one cell when it finished nowhere there: its lower
# bound (from a timeout), or the reason it has no time.
function arm_cell_text(d, b, a,   here, n, i, r, s, low, any) {
    n = split(repeat_list[d, b], here, " ")
    low = ""; any = ""
    for (i = 1; i <= n; ++i) {
        r = here[i]
        s = state[d, b, r, a]
        if (s == "") { continue }
        if (s == "timeout") {
            any = "timeout"
            if ((d, b) in outside_max && limit[d, b, r, a] * 1000 - outside_max[d, b] > 0) {
                if (low == "" || limit[d, b, r, a] * 1000 - outside_max[d, b] < low) {
                    low = limit[d, b, r, a] * 1000 - outside_max[d, b]
                }
            }
        } else if (any == "") {
            any = (s == "skipped") ? "skipped" : (s == "oom") ? "out_of_memory" : "failed"
        }
    }
    if (low != "") { return ">=" fmt_ms(low) }
    return (any == "") ? "-" : any
}
# The last compare() as text: its ratio, a lower bound, or why there is none.
function speedup_text() {
    return (X_status == "finished") ? fmt_ratio(X_speedup) : \
        (X_status == "from_timed_out") ? ">=" fmt_ratio(X_speedup) : \
        (X_status == "missing") ? "-" : X_status
}
# The fastest of edge, vertex and the method in one cell by median time.  An
# arm that timed out, or was skipped after a timeout, is slower than any that
# finished; "-" when one failed or ran out of memory, or none finished.
function fastest_arm(d, b,   names, k, a, m, best, best_arm, text) {
    split("edge vertex " method_arm, names, " ")
    best = ""; best_arm = "-"
    for (k = 1; k <= 3; ++k) {
        a = names[k]
        m = arm_median(d, b, a)
        if (m == "") {
            text = arm_cell_text(d, b, a)
            if (text !~ /^>=/ && text != "timeout" && text != "skipped") { return "-" }
            continue
        }
        if (best == "" || m < best) { best = m; best_arm = a }
    }
    return best_arm
}
# One arm of one cell for the time table.
function time_text(d, b, a,   m) {
    m = arm_median(d, b, a)
    return (m != "") ? fmt_ms(m) : arm_cell_text(d, b, a)
}
# Why one arm of one cell has no search statistics: what its runs did.
function status_text(d, b, a,   text) {
    text = arm_cell_text(d, b, a)
    return (text ~ /^>=/) ? "timeout" : text
}
# The branches of one arm in one cell, "" without a finished run.
function branch_count(d, b, a) {
    return ((d, b, a) in branches) ? branches[d, b, a] : ""
}
function branch_ratio(d, b, from, to,   f, t) {
    f = branch_count(d, b, from); t = branch_count(d, b, to)
    if (f == "" || t == "" || t + 0 <= 0) { return "" }
    return f / t
}
# One comparison in one cell, into X_status, X_fin, X_bound, X_from, X_to,
# X_speedup (a median of ratios, or of bounds when none finished) and
# X_faster.
function compare(d, b, from, to,   here, n, i, r, sf, st, ratios, bounds, from_ms, to_ms, bound_to, n_up, n_skip, n_fail, n_open, n_oom_from, n_oom_to, room) {
    X_fin = 0; X_bound = 0; X_faster = 0
    X_from = ""; X_to = ""; X_speedup = ""
    ratios = ""; bounds = ""; from_ms = ""; to_ms = ""; bound_to = ""
    n_up = 0; n_skip = 0; n_fail = 0; n_open = 0; n_oom_from = 0; n_oom_to = 0
    n = split(repeat_list[d, b], here, " ")
    for (i = 1; i <= n; ++i) {
        r = here[i]
        sf = state[d, b, r, from]; st = state[d, b, r, to]
        if (sf == "" || st == "") { continue }
        if (sf == "ok" && st == "ok") {
            ++X_fin
            from_ms = from_ms " " value[d, b, r, from]
            to_ms = to_ms " " value[d, b, r, to]
            if (value[d, b, r, to] > 0) {
                ratios = ratios " " (value[d, b, r, from] / value[d, b, r, to])
            }
            if (value[d, b, r, to] < value[d, b, r, from]) { ++X_faster }
        } else if (sf == "timeout" && st == "ok") {
            room = ((d, b) in outside_max) ? limit[d, b, r, from] * 1000 - outside_max[d, b] : 0
            if (room <= 0 || value[d, b, r, to] <= 0) { ++n_open; continue }
            ++X_bound
            bounds = bounds " " (room / value[d, b, r, to])
            bound_to = bound_to " " value[d, b, r, to]
            if (X_bound_from == "" || room < X_bound_from) { X_bound_from = room }
        } else if (sf == "ok" && st == "timeout") {
            ++n_up
        } else if (sf == "oom") {
            ++n_oom_from
        } else if (st == "oom") {
            ++n_oom_to
        } else if (sf == "skipped" || st == "skipped") {
            ++n_skip
        } else if (sf == "timeout" || st == "timeout") {
            ++n_open
        } else {
            ++n_fail
        }
    }
    if (X_fin > 0) {
        X_status = "finished"
        X_from = median(from_ms); X_to = median(to_ms); X_speedup = median(ratios)
    } else if (X_bound > 0) {
        X_status = "from_timed_out"
        X_to = median(bound_to); X_speedup = median(bounds)
    } else {
        X_status = (n_up > 0) ? "to_timed_out" : \
            (n_oom_from > 0) ? "from_out_of_memory" : \
            (n_oom_to > 0) ? "to_out_of_memory" : \
            (n_open > 0) ? "timed_out_unbounded" : \
            (n_skip > 0) ? "skipped_after_timeout" : \
            (n_fail > 0) ? "failed" : "missing"
    }
}
# The log-speedup of one cell into its batch size, its group and the whole
# grid; finished cells count in both means, timed-out ones in the bound only.
function add(e, g, b, logv, is_finished,   keys, key) {
    note_batch(e, g, b)
    split(e SUBSEP g SUBSEP b " " e SUBSEP g SUBSEP "all" " " e SUBSEP "all" SUBSEP "all", keys, " ")
    for (key = 1; key <= 3; ++key) {
        if (is_finished) { fin_sum[keys[key]] += logv; ++fin_count[keys[key]] }
        else { bound_sum[keys[key]] += logv; ++bound_count[keys[key]] }
    }
}
function add_branches(e, g, b, logv,   keys, key) {
    note_batch(e, g, b)
    split(e SUBSEP g SUBSEP b " " e SUBSEP g SUBSEP "all" " " e SUBSEP "all" SUBSEP "all", keys, " ")
    for (key = 1; key <= 3; ++key) { br_sum[keys[key]] += logv; ++br_count[keys[key]] }
}
function note_batch(e, g, b) {
    if (!((e, g, b) in seen_experiment_batch)) {
        seen_experiment_batch[e, g, b] = 1
        experiment_batches[e, g] = experiment_batches[e, g] " " b
    }
}
# A cell that enters no mean (skipped, out of memory), counted the same way.
function count_cell(counter, e, g, b) {
    ++counter[e, g, b]; ++counter[e, g, "all"]; ++counter[e, "all", "all"]
    note_batch(e, g, b)
}
BEGIN {
    arm_count = split(arm_list, arm_order, " ")
    experiment_count = split(experiment_list, entry, " ")
    for (e = 1; e <= experiment_count; ++e) {
        split(entry[e], name_and_spec, "=")
        split(name_and_spec[2], spec, ":")
        experiment_name[e] = name_and_spec[1]
        experiment_from[e] = spec[1]
        experiment_to[e] = spec[2]
    }
}
NR == 1 {
    for (i = 1; i <= NF; ++i) { column[$i] = i }
    next
}
{
    d = $(column["dataset"]); b = $(column["batch_size"])
    r = $(column["repeat"]); a = $(column["arm"]); s = $(column["status"])
    class_of[d] = $(column["dataset_class"])
    if (!((d, b) in seen_cell)) {
        seen_cell[d, b] = 1
        cell_key[++cells] = d SUBSEP b
    }
    if (!((d, b, r) in seen_repeat)) {
        seen_repeat[d, b, r] = 1
        repeat_list[d, b] = repeat_list[d, b] " " r
    }
    limit[d, b, r, a] = $(column["time_limit_seconds"]) + 0
    if (s == "ok" && $(column["clique_counts_match"]) != "no" && \
        $(column["switches_check"]) == "ok" && number($(column["new_clique_ms"])) && \
        count($(column["branches"]))) {
        state[d, b, r, a] = "ok"
        value[d, b, r, a] = $(column["new_clique_ms"]) + 0
        if (count($(column["processed_batches"]))) { cell_batches[d, b] = $(column["processed_batches"]) }
        wall[d, b, r, a] = $(column["wall_seconds"]) + 0
        outside = wall[d, b, r, a] * 1000 - value[d, b, r, a]
        if (!((d, b) in outside_max) || outside > outside_max[d, b]) { outside_max[d, b] = outside }
        # The search: deterministic, so the first finished run gives it and
        # the later ones must agree (branches_min / branches_max).
        n_branches = $(column["branches"]) + 0
        if (!((d, b, a) in branches)) {
            branches[d, b, a] = n_branches
            branches_min[d, b, a] = n_branches
            branches_max[d, b, a] = n_branches
            tree_roots[d, b, a] = $(column["search_roots"])
            tree_below[d, b, a] = $(column["branches_below_roots"])
            tree_leaves[d, b, a] = $(column["leaves"])
            tree_emitted[d, b, a] = $(column["emitted"])
            tree_pivot_calls[d, b, a] = $(column["pivot_calls"])
            tree_pivot_candidates[d, b, a] = $(column["pivot_candidates"])
            tree_pruned_roots[d, b, a] = $(column["vertex_pruned_roots"])
            tree_new_edge_prunes[d, b, a] = $(column["vertex_new_edge_prunes"])
            tree_partner_prunes[d, b, a] = $(column["star_partner_prunes"])
        }
        if (n_branches < branches_min[d, b, a]) { branches_min[d, b, a] = n_branches }
        if (n_branches > branches_max[d, b, a]) { branches_max[d, b, a] = n_branches }
        # graph-update-ms is the same work in every run of a cell.
        update = $(column["graph_update_ms"])
        if (number(update)) {
            update += 0
            if (!((d, b, r) in update_runs)) {
                update_runs[d, b, r] = 0
                update_min[d, b, r] = update
                update_max[d, b, r] = update
            }
            ++update_runs[d, b, r]
            if (update < update_min[d, b, r]) { update_min[d, b, r] = update }
            if (update > update_max[d, b, r]) { update_max[d, b, r] = update }
        }
    } else if (s == "timeout") {
        state[d, b, r, a] = "timeout"
    } else if (s == "skipped_after_timeout") {
        state[d, b, r, a] = "skipped"
    } else if (s == "out_of_memory") {
        state[d, b, r, a] = "oom"
    } else {
        state[d, b, r, a] = "failed"
    }
}
END {
    # table_time.tsv: per cell the median time of every arm (or what stopped
    # it) and every comparison; table_tree.tsv: the branches of every arm and
    # every comparison of them; tree.tsv: one row per arm with its search
    # statistics.
    printf "dataset\tbatch_size" > table_time_file
    printf "dataset\tbatch_size" > table_tree_file
    for (k = 1; k <= arm_count; ++k) {
        printf "\t%s_ms", arm_order[k] > table_time_file
        printf "\t%s_branches", arm_order[k] > table_tree_file
    }
    for (e = 1; e <= experiment_count; ++e) {
        printf "\t%s", experiment_name[e] > table_time_file
        printf "\t%s", experiment_name[e] > table_tree_file
    }
    printf "\n" > table_time_file
    printf "\n" > table_tree_file
    print "dataset", "batch_size", "arm", "branches", "branches_vs_" method_arm, \
        "search_roots", "branches_below_roots", "leaves", "emitted", "pivot_calls", \
        "pivot_candidates", "pruned_roots", "new_edge_prunes", "partner_prunes", \
        "runs_agree" > tree_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        printf "%s\t%s", d, b > table_time_file
        printf "%s\t%s", d, b > table_tree_file
        for (k = 1; k <= arm_count; ++k) {
            a = arm_order[k]
            printf "\t%s", time_text(d, b, a) > table_time_file
            printf "\t%s", (((d, b, a) in branches) ? branches[d, b, a] : status_text(d, b, a)) > table_tree_file
        }
        for (e = 1; e <= experiment_count; ++e) {
            X_bound_from = ""
            compare(d, b, experiment_from[e], experiment_to[e])
            text = (X_status == "finished") ? fmt_ratio(X_speedup) : \
                (X_status == "from_timed_out") ? ">=" fmt_ratio(X_speedup) : \
                (X_status == "missing") ? "-" : X_status
            printf "\t%s", text > table_time_file
            ratio = branch_ratio(d, b, experiment_from[e], experiment_to[e])
            printf "\t%s", fmt_ratio(ratio) > table_tree_file
        }
        printf "\n" > table_time_file
        printf "\n" > table_tree_file
        for (k = 1; k <= arm_count; ++k) {
            a = arm_order[k]
            if (!((d, b, a) in branches)) { continue }
            print d, b, a, branches[d, b, a], fmt_ratio(branch_ratio(d, b, a, method_arm)), \
                tree_roots[d, b, a], tree_below[d, b, a], \
                tree_leaves[d, b, a], tree_emitted[d, b, a], \
                tree_pivot_calls[d, b, a], tree_pivot_candidates[d, b, a], \
                tree_pruned_roots[d, b, a], tree_new_edge_prunes[d, b, a], \
                tree_partner_prunes[d, b, a], \
                (branches_min[d, b, a] == branches_max[d, b, a] ? "yes" : "no") > tree_file
        }
    }

    # times.tsv: one row per cell and arm: what its runs did, the time of
    # every repeat, then the median, the extremes and their ratio.
    printf "dataset\tbatch_size\tarm\tfinished\ttimed_out\tskipped\tout_of_memory\tfailed" > times_file
    for (r = 1; r <= repeats; ++r) { printf "\tr%d_ms", r > times_file }
    printf "\tmedian_ms\tmin_ms\tmax_ms\tspread\tmedian_wall_seconds\n" > times_file
    print "dataset", "batch_size", "experiment", "from_arm", "to_arm", "status", \
        "finished_pairs", "bounded_pairs", "from_ms", "to_ms", "speedup", \
        "to_faster_in", "from_branches", "to_branches", "branch_ratio" > results_file
    print "dataset", "batch_size", "repeat", "runs", "graph_update_min_ms", \
        "graph_update_max_ms", "spread" > noise_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        for (m = 1; m <= arm_count; ++m) {
            a = arm_order[m]
            n_ok = 0; n_timeout = 0; n_skip = 0; n_oom = 0; n_fail = 0
            time_list = ""; wall_list = ""; per_repeat = ""
            low = ""; high = ""
            for (r = 1; r <= repeats; ++r) {
                s = state[d, b, r, a]
                if (s == "ok") {
                    t = value[d, b, r, a]
                    ++n_ok
                    time_list = time_list " " t
                    wall_list = wall_list " " wall[d, b, r, a]
                    per_repeat = per_repeat "\t" fmt_ms(t)
                    if (low == "" || t < low) { low = t }
                    if (high == "" || t > high) { high = t }
                } else {
                    if (s == "timeout") { ++n_timeout }
                    else if (s == "skipped") { ++n_skip }
                    else if (s == "oom") { ++n_oom }
                    else if (s == "failed") { ++n_fail }
                    per_repeat = per_repeat "\t" ((s == "") ? "-" : s)
                }
            }
            if (n_ok + n_timeout + n_skip + n_oom + n_fail == 0) { continue }
            printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d%s\t%s\t%s\t%s\t%s\t%s\n", d, b, a, \
                n_ok, n_timeout, n_skip, n_oom, n_fail, per_repeat, \
                fmt_ms(median(time_list)), fmt_ms(low), fmt_ms(high), \
                ((low != "" && low > 0) ? sprintf("%.3f", high / low) : "-"), \
                ((wall_list != "") ? sprintf("%.3f", median(wall_list)) : "-") > times_file
        }
        n = split(repeat_list[d, b], here, " ")
        for (i = 1; i <= n; ++i) {
            r = here[i]
            if (!((d, b, r) in update_runs)) { continue }
            print d, b, r, update_runs[d, b, r], fmt_ms(update_min[d, b, r]), \
                fmt_ms(update_max[d, b, r]), \
                (update_min[d, b, r] > 0 ? sprintf("%.3f", update_max[d, b, r] / update_min[d, b, r]) : "-") > noise_file
        }
        for (e = 1; e <= experiment_count; ++e) {
            from = experiment_from[e]; to = experiment_to[e]
            X_bound_from = ""
            compare(d, b, from, to)
            if (X_status == "missing") { continue }
            from_branches = branch_count(d, b, from)
            to_branches = branch_count(d, b, to)
            ratio = branch_ratio(d, b, from, to)
            if (X_status == "finished") {
                print d, b, experiment_name[e], from, to, X_status, X_fin, X_bound, \
                    fmt_ms(X_from), fmt_ms(X_to), fmt_ratio(X_speedup), X_faster "/" X_fin, \
                    (from_branches == "" ? "-" : from_branches), \
                    (to_branches == "" ? "-" : to_branches), fmt_ratio(ratio) > results_file
                if (X_speedup != "" && X_speedup > 0) { add(e, class_of[d], b, log(X_speedup), 1) }
            } else if (X_status == "from_timed_out") {
                print d, b, experiment_name[e], from, to, X_status, 0, X_bound, \
                    ">=" fmt_ms(X_bound_from), fmt_ms(X_to), ">=" fmt_ratio(X_speedup), \
                    X_bound "/" X_bound, "-", \
                    (to_branches == "" ? "-" : to_branches), "-" > results_file
                add(e, class_of[d], b, log(X_speedup), 0)
            } else {
                print d, b, experiment_name[e], from, to, X_status, 0, 0, "-", "-", "-", "-", \
                    (from_branches == "" ? "-" : from_branches), \
                    (to_branches == "" ? "-" : to_branches), fmt_ratio(ratio) > results_file
                if (X_status == "skipped_after_timeout") { count_cell(skipped_cells, e, class_of[d], b) }
                if (X_status ~ /out_of_memory/) { count_cell(oom_cells, e, class_of[d], b) }
            }
            if (ratio != "" && ratio > 0) { add_branches(e, class_of[d], b, log(ratio)) }
        }
    }

    # crossover.tsv and crossover_summary.tsv: per dataset, batch sizes
    # ascending, vertex against edge; with a third arm run as METHOD_ARM, that
    # arm against both as well.
    for (k = 1; k <= arm_count; ++k) { arm_run[arm_order[k]] = 1 }
    crossover_arms = ("edge" in arm_run) && ("vertex" in arm_run)
    with_method = crossover_arms && (method_arm in arm_run) && \
        method_arm != "edge" && method_arm != "vertex"
    line = "dataset" OFS "batch_size" OFS "batches" OFS "edge_ms" OFS "vertex_ms" OFS \
        "vertex_vs_edge" OFS "faster_baseline"
    if (with_method) {
        line = line OFS method_arm "_ms" OFS "edge_vs_" method_arm OFS \
            "vertex_vs_" method_arm OFS "fastest"
    }
    print line > crossover_file
    line = "dataset" OFS "batch_sizes" OFS "smallest_batch" OFS "largest_batch" OFS \
        "vertex_vs_edge_smallest" OFS "vertex_vs_edge_largest" OFS \
        "last_batch_edge_faster" OFS "first_batch_vertex_faster" OFS \
        "crossing_estimate" OFS "sizes_edge_faster" OFS "sizes_vertex_faster" OFS "crossed"
    if (with_method) {
        line = line OFS method_arm "_fastest_cells" OFS "min_edge_vs_" method_arm OFS \
            "min_vertex_vs_" method_arm
    }
    print line > crossover_summary_file
    dataset_count = 0
    for (c = 1; c <= cells && crossover_arms; ++c) {
        split(cell_key[c], parts, SUBSEP)
        if (!(parts[1] in dataset_batches)) { dataset_name[++dataset_count] = parts[1] }
        dataset_batches[parts[1]] = dataset_batches[parts[1]] " " parts[2]
    }
    for (q = 1; q <= dataset_count; ++q) {
        d = dataset_name[q]
        n = split(dataset_batches[d], sizes, " ")
        for (i = 2; i <= n; ++i) {
            key = sizes[i] + 0
            for (j = i - 1; j >= 1 && sizes[j] + 0 > key; --j) { sizes[j + 1] = sizes[j] }
            sizes[j + 1] = key
        }
        last_edge = ""; first_vertex = ""; back = 0
        edge_sizes = 0; vertex_sizes = 0
        fastest_cells = 0; decided_cells = 0
        min_edge = ""; min_vertex = ""; first_text = "-"; last_text = "-"
        split("", ve_ratio)
        for (i = 1; i <= n; ++i) {
            b = sizes[i]
            X_bound_from = ""
            compare(d, b, "vertex", "edge")
            ve_text = speedup_text()
            if (X_status == "finished" && X_speedup > 0) { ve_ratio[b] = X_speedup }
            # Which baseline was faster: vertex 1, edge 0, undecided -1.  A
            # vertex run skipped after a timeout at a larger size counts as
            # slower, as the skip itself assumes.
            side = (X_status == "finished") ? (X_speedup < 1 ? 1 : 0) : \
                (X_status == "from_timed_out" && X_speedup >= 1) ? 0 : \
                (X_status == "skipped_after_timeout" && arm_median(d, b, "edge") != "") ? 0 : -1
            if (i == 1) { first_text = ve_text }
            if (i == n) { last_text = ve_text }
            if (side == 0) {
                ++edge_sizes
                if (first_vertex == "") { last_edge = b } else { back = 1 }
            } else if (side == 1) {
                ++vertex_sizes
                if (first_vertex == "") { first_vertex = b }
            }
            line = d OFS b OFS (((d, b) in cell_batches) ? cell_batches[d, b] : "-") OFS \
                time_text(d, b, "edge") OFS time_text(d, b, "vertex") OFS ve_text OFS \
                ((side == 1) ? "vertex" : (side == 0) ? "edge" : "-")
            if (with_method) {
                X_bound_from = ""
                compare(d, b, "edge", method_arm)
                em_text = speedup_text()
                if (X_status == "finished" && (min_edge == "" || X_speedup < min_edge)) { min_edge = X_speedup }
                X_bound_from = ""
                compare(d, b, "vertex", method_arm)
                vm_text = speedup_text()
                if (X_status == "finished" && (min_vertex == "" || X_speedup < min_vertex)) { min_vertex = X_speedup }
                fastest = fastest_arm(d, b)
                if (fastest != "-") {
                    ++decided_cells
                    if (fastest == method_arm) { ++fastest_cells }
                }
                line = line OFS time_text(d, b, method_arm) OFS em_text OFS vm_text OFS fastest
            }
            print line > crossover_file
        }
        crossed = (first_vertex == "") ? ((last_edge == "") ? "undecided" : "no") : \
            (last_edge == "") ? "vertex_faster_throughout" : \
            (back ? "yes_not_monotone" : "yes")
        # Where vertex_vs_edge reaches 1, on log scales, between the last size
        # edge won and the first vertex won.
        estimate = "-"
        if (last_edge != "" && first_vertex != "" && (last_edge in ve_ratio) && \
            (first_vertex in ve_ratio) && ve_ratio[last_edge] != ve_ratio[first_vertex]) {
            t = log(ve_ratio[last_edge]) / (log(ve_ratio[last_edge]) - log(ve_ratio[first_vertex]))
            estimate = sprintf("%.0f", exp(log(last_edge) + t * (log(first_vertex) - log(last_edge))))
        }
        line = d OFS n OFS sizes[1] OFS sizes[n] OFS first_text OFS last_text OFS \
            (last_edge == "" ? "-" : last_edge) OFS (first_vertex == "" ? "-" : first_vertex) OFS \
            estimate OFS edge_sizes OFS vertex_sizes OFS crossed
        if (with_method) {
            line = line OFS fastest_cells "/" decided_cells OFS fmt_ratio(min_edge) OFS \
                fmt_ratio(min_vertex)
        }
        print line > crossover_summary_file
    }

    print "experiment", "from_arm", "to_arm", "datasets", "batch_size", \
        "finished_cells", "speedup_geomean_finished", "timed_out_cells", \
        "speedup_geomean_lower_bound", "skipped_cells", "out_of_memory_cells", \
        "branch_cells", "branch_ratio_geomean" > geomean_file
    group_count = split("small large", group_name, " ")
    group_name[++group_count] = "all"
    for (e = 1; e <= experiment_count; ++e) {
        for (q = 1; q <= group_count; ++q) {
            g = group_name[q]
            batch_count = 0
            if (g != "all") { batch_count = split(experiment_batches[e, g], batches_here, " ") }
            batches_here[++batch_count] = "all"
            for (z = 1; z <= batch_count; ++z) {
                b = batches_here[z]
                if (fin_count[e, g, b] + bound_count[e, g, b] + skipped_cells[e, g, b] + \
                    oom_cells[e, g, b] + br_count[e, g, b] == 0) { continue }
                geo_fin = (fin_count[e, g, b] > 0) ? \
                    sprintf("%.4f", exp(fin_sum[e, g, b] / fin_count[e, g, b])) : "-"
                both = fin_count[e, g, b] + bound_count[e, g, b]
                geo_low = (both > 0) ? \
                    sprintf("%.4f", exp((fin_sum[e, g, b] + bound_sum[e, g, b]) / both)) : "-"
                geo_br = (br_count[e, g, b] > 0) ? \
                    sprintf("%.4f", exp(br_sum[e, g, b] / br_count[e, g, b])) : "-"
                print experiment_name[e], experiment_from[e], experiment_to[e], g, b, \
                    fin_count[e, g, b] + 0, geo_fin, bound_count[e, g, b] + 0, geo_low, \
                    skipped_cells[e, g, b] + 0, oom_cells[e, g, b] + 0, \
                    br_count[e, g, b] + 0, geo_br > geomean_file
            }
        }
    }
}
' "$SUMMARY_FILE"

{
    echo "finished=$(date -Is)"
    echo "run_count=$run_count"
    echo "failures=$failures"
    echo "timeouts=$timeouts"
    echo "out_of_memory_runs=$out_of_memory_runs"
    echo "skipped=$skipped"
    echo "elapsed_seconds=$(( $(date +%s) - JOB_START ))"
    echo "table_time_file=$TABLE_TIME_FILE"
    echo "table_tree_file=$TABLE_TREE_FILE"
    echo "results_file=$RESULTS_FILE"
    echo "geomean_file=$GEOMEAN_FILE"
    echo "times_file=$TIMES_FILE"
    echo "tree_file=$TREE_FILE"
    echo "noise_file=$NOISE_FILE"
    echo "crossover_file=$CROSSOVER_FILE"
    echo "crossover_summary_file=$CROSSOVER_SUMMARY_FILE"
    echo "timeout_file=$TIMEOUT_FILE"
    echo "summary_file=$SUMMARY_FILE"
} >> "${OUTPUT_DIR}/manifest.txt"

show_table() {
    if command -v column > /dev/null 2>&1; then
        column -t -s $'\t' "$1"
    else
        cat "$1"
    fi
}

echo "Finished: $(date -Is)"
echo "Runs: $run_count; failures: $failures; timeouts: $timeouts; out of memory: $out_of_memory_runs; skipped: $skipped"
echo "Time table: $TABLE_TIME_FILE"
echo "Branch table: $TABLE_TREE_FILE"
echo "Times: $TIMES_FILE"
echo "Search statistics: $TREE_FILE"
echo "Results: $RESULTS_FILE"
echo "Geometric means: $GEOMEAN_FILE"
echo "Noise: $NOISE_FILE"
echo "Crossover: $CROSSOVER_FILE"
echo "Crossover summary: $CROSSOVER_SUMMARY_FILE"
echo "Timeouts: $TIMEOUT_FILE"
echo "Summary: $SUMMARY_FILE"
echo
echo "Median new-clique time (ms) of every method (>= a lower bound from a timeout), and the speedups time(from)/time(to) (above 1: to is faster):"
show_table "$TABLE_TIME_FILE"
echo
echo "Branches of every method, and branches(from)/branches(to) (above 1: to opens fewer):"
show_table "$TABLE_TREE_FILE"
echo
echo "Geometric means of the per-cell speedups and branch ratios:"
show_table "$GEOMEAN_FILE"
echo
echo "Where vertex overtakes edge, per dataset (vertex_vs_edge = time(vertex)/time(edge): below 1, vertex is faster):"
show_table "$CROSSOVER_SUMMARY_FILE"
echo
echo "Per cell, batch sizes ascending:"
show_table "$CROSSOVER_FILE"
echo
if grep -q $'\tno$' "$TREE_FILE"; then
    echo "Warning: the runs of some arm disagreed on its search (runs_agree=no in tree.tsv):"
    grep $'\tno$' "$TREE_FILE"
    echo
fi
echo "graph-update-ms over the arms of each cell (the same work; a spread well above 1 marks a disturbed run):"
show_table "$NOISE_FILE"
if (( $(wc -l < "$TIMEOUT_FILE") > 1 )); then
    echo
    echo "Runs that hit their time limit (their time is at least the limit):"
    show_table "$TIMEOUT_FILE"
fi
if (( failures != 0 )); then
    exit 1
fi
