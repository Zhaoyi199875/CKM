#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --job-name=immce_overall_hybrid_coauthor_long
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/zzhong/zzy_project/DynamicClique/IMMCE_SubsumedClique/output/overall_hybrid_coauthor_long_%j.out
#SBATCH --error=/home/zzhong/zzy_project/DynamicClique/IMMCE_SubsumedClique/output/overall_hybrid_coauthor_long_%j.err
# OzSTAR refuses --exclusive (NODE below).

# LONG-COAUTHOR VERSION (2026-10-30) of submit_overall_hybrid_times.sh (the
# paper version of submit_mixed_full_times.sh), standalone: the hybrid-batch
# cells of coauthor in the overall comparison of the paper's Section VI-B,
# with SOMEi given days rather than hours.  In job 19365025 it ran out of
# twenty hours at batch 200000 after 12 of its 21 batches, while inserting
# and deleting it finished every coauthor cell.  Everything is that
# script's, except:
#
#   * one PART, coauthor-long: edge-opt, star-rk and somei on coauthor in
#     one job, in that order, at one batch size (200000 unless BATCH_SIZES
#     gives another), so that SOMEi's ratios to the other two are formed in
#     the job and SOMEi has the rest of the allocation (BUDGET below);
#   * the allocation is seven days, the deadline half an hour before its end
#     whatever sbatch --time gives (DEADLINE_SECONDS below), and SOMEi's own
#     limit thirty days, so that the deadline stops it;
#   * its own job name, output files and experiment name.
#
# That script's own changes stand: SOMEi builds the starting graph in one
# batch (--somei-build=one-batch), prints a line per measured batch
# (--progress; somei_progress.tsv shows how far a run its limit stopped had
# got), and the source marker of 2026-10-27.  Since 2026-11-06 that build's
# steps find the maximal cliques of G[R] by Bron-Kerbosch
# (--somei-build-update=bk, one more source marker): the very state the
# paper's construction leaves; the measured batches run the paper's steps.
# And SOMEi's limit bounds its measured batches alone (--measured-time-limit),
# with SOMEI_SETUP_SECONDS more for the reading and the build (TIME LIMITS); a
# run set to finish within SOMEI_EXTEND_SECONDS runs on until then
# (--measured-time-extend), which under thirty days and the deadline never
# comes into play.
#
# It calls no other script.
#
# Mixed batches, each deleting and inserting edges at once: every method's
# running time, split into its new-clique and subsumed searches, and the
# branches of the new-clique searches of star, edge and vertex.
#
# THE WORKLOAD (WORKLOAD, default sliding-window) decides which edges a batch
# inserts and which it deletes:
#
#   sliding-window  --update-mode=sliding-window (see README): the graph
#                   keeps the newest W edges.  A run uses the dataset's first
#                   T edges (its edge budget below, --max-edges; the whole
#                   stream when 0) and W = T/2: the first W go in untimed (the
#                   build), and every batch then inserts the next batch-size
#                   edges and deletes as many of the oldest.  The inserted
#                   edges are a contiguous stretch of the stream and so are
#                   the deleted ones, so no event is split, and the graph is
#                   always a contiguous slice of time.  The last batch may be
#                   short; every batch size replays the same edges.
#   random-labels   --update-mode=random-labels (see README),
#                   the fully dynamic workload of Das et al. (section 6.1):
#                   every edge of the window is labelled insert with
#                   probability INSERT_FRACTION (default 0.5) and delete
#                   otherwise, by a generator seeded with LABEL_SEED (default
#                   1), in stream order.  The build, untimed, inserts the
#                   delete-labelled edges; then every batch inserts its
#                   insert-labelled edges and deletes its delete-labelled
#                   ones.  The i-th edge of the window always takes the i-th
#                   draw, so every batch size and every arm sees the same
#                   labels.  Das labelled two RMAT graphs with random
#                   timestamps; on a real stream the labels split every event
#                   in two and thin the graph at random, which takes away
#                   star's sharing (experiments/mixed_lead_20261010/NOTE.md).
#
# Either way a batch takes effect as a whole -- the graph after it is the
# graph before it less its deletions plus its insertions, the semantics of
# both papers -- and only the batches after the build are timed.
#
# THE STRATEGY.  Both give the same changes and differ in how they get them
# (designs/mixed_semantics_plan_20261010.md):
#
#   single-pass   ours: the cliques the batch destroys (D) are found on the
#                 graph before it and its new cliques (N) on the graph after
#                 it; the subsumed ones (S) come from N, and the cliques the
#                 deletions create (A) are recovered from D, each checked for
#                 maximality on the graph after the batch.  No clique that is
#                 in neither graph ever forms.
#   two-phase     Das et al., section 4.5: the deletions first -- D, then A
#                 recovered on the graph without the deleted edges -- and then
#                 the insertions -- N, then S.  The cliques the first phase
#                 creates and the second subsumes (transient) are netted out
#                 at the end of the batch.
#
# Every arm runs two-phase unless it says -singlepass (below): since
# 2026-10-11 it is the default of the code too, star-rk's single pass having
# come out slower.
#
# ONE RUN of a method times, over the batches,
#
#   total         timing: search-and-maintenance-ms
#                 = new clique + subsumed + index update + netting
#   new clique    <star|edge|vertex>: time-ms.  The new-clique search, run
#                 twice per batch: with the deleted edges as its new edges on
#                 the graph before the batch (finds D), and with the inserted
#                 ones on the graph after it (finds N); mixed.tsv has the two
#                 apart
#   subsumed      subsumed: time-ms.  The subsumed search from N (finds S)
#                 and the recovery from D (finds A); mixed.tsv has the two
#                 apart
#   index update  timing: clique-index-update-ms
#   netting       mixed: transient-netting-ms, two-phase runs only
#
# Reading the input, the build and the graph updates (graph-update-ms, the
# same work in every run) are outside the total.  SOMEi updates a batch as a
# whole -- the inserted edges added, its scope read on the graph holding both
# kinds of edge, the deleted edges removed, its headers rolled back and
# redone -- so it has no new-clique and subsumed parts: those columns are "-",
# and its scope, rebuild and change listing times are reported instead.
#
# THE ARMS are <search>-<subsumed>[-singlepass]:
#
#   search    star        the method: vertex-cover stars, the branch-count rule
#                         at every state whose clique holds no new edge yet
#                         (--star-branch-scope=unbound, the code default)
#             star-root   the rule at each star's centre only (root)
#             edge        one search per batch edge over the common
#                         neighbourhood of its endpoints
#             vertex      one search per batch endpoint over its closed
#                         neighbourhood N[s], the classic Tomita BK with the
#                         new-edge condition pruning it (time-limited)
#   subsumed  rk          --subsumed-algorithm=residual-kernel
#             rk-nogate   the same with --subsumed-size-histogram=off: the
#                         size gate inert, as star-rk-nogate in the overall job
#             optbm       --subsumed-algorithm=optimized-bitmap
#             opt         --subsumed-algorithm=optimized
#   -singlepass           --mixed-strategy=single-pass; without it two-phase
#
# and somei[-noprune][-nochg] (Yu et al., TKDE 2024; see
# submit_somei_times.sh), which runs on SOMEI_DATASETS only (below).  The
# default ARMS (2026-10-12, the user's choice) are
#
#   edge-opt star-rk somei
#
# all two-phase.  star-rk is the method (METHOD_ARM), edge-opt the
# edge-seeded search with the optimized subsumed search, SOMEi the published
# baseline; edge-optbm, vertex-optbm and the -singlepass arms can still be
# named in ARMS.  The build is the same in every run but SOMEi's: star with rk,
# whatever the arm (any exact pair leaves the same cliques), the index
# keeping its size histogram throughout when the arm reads one; SOMEi builds
# its own structure.  Lists take entries separated by spaces or commas.
#
# COMPARISONS (results.tsv, table_speedup.tsv, geomeans.tsv), each
# time(from)/time(to) on one part, above 1 meaning `to` is faster:
#
#   <arm>_vs_<method>_total          every other arm -> METHOD_ARM
#   <arm>_vs_<method>_new_clique     search arms, with the branch ratio
#                                    branches(from)/branches(to)
#   <arm>_vs_<method>_subsumed       search arms
#   <arm>_vs_<twin>_total            a search arm whose strategy is not the
#   <arm>_vs_<twin>_subsumed         method's -> the method's searches with
#                                    that strategy, when that arm runs (e.g.
#                                    edge-optbm-singlepass ->
#                                    star-rk-singlepass); the new-clique
#                                    search does not depend on the strategy
#   <arm>_singlepass_vs_<arm>_total  the two strategies of one search and
#   <arm>_singlepass_vs_<arm>_subsumed  subsumed search when neither is the
#                                    method (whose comparisons cover its pair)
#   <vertex>_vs_<edge>_new_clique    the two baselines' searches
#
# With the default arms: edge-opt against star-rk on the total, new-clique
# and subsumed times, the new-clique one with the branch ratio of edge to
# star; and SOMEi against star-rk on the total.
#
# EXPERIMENTS replaces the list: name=from:to:part entries, part total,
# new_clique, subsumed or index_update (the last three need two search arms).
# A cell's speedup is the median of its paired per-repeat ratios.
#
# BRANCHES, as in submit_edge_vertex_star_times.sh: the recursive calls of the
# new-clique search over the batches, both of its runs per batch, every state
# (R, P, X) it entered, roots and nodes cut on entry included (search_roots:
# edge one per batch edge; vertex one per root that passed the new-edge test
# before its graph was built; star one per single-edge star, one per partner
# of a star whose centre took the new edges, one per star whose centre took
# the pivot).  tree.tsv also gives them by side: on the deleted edges
# (finding D) and on the inserted ones (finding N).  Neither the strategy nor
# the subsumed search changes them (locally every arm with the same search
# reported the same counts), and the search is deterministic, so the runs of
# one arm must agree on it (tree.tsv).
#
# TIME LIMITS, wall time of a whole run, the build included -- except for
# SOMEi since 2026-11-06: its SOMEI_TIMEOUT_SECONDS bound the measured batches
# alone (--measured-time-limit, counted from where they begin, the run then
# ending as a timeout), and the run as a whole gets SOMEI_SETUP_SECONDS
# (default 3600) more for its reading and untimed build, so that a long build
# cannot eat into the measured time; and a run whose measured batches look set
# to finish within SOMEI_EXTEND_SECONDS (default SOMEI_TIMEOUT_SECONDS, off:
# the time so far plus the last batch's time for every batch left, judged
# after each batch) runs on until then (--measured-time-extend), its wall
# limit growing to match:
#
#   vertex arms   VERTEX_TIMEOUT_SECONDS (default 3600): the vertex search
#                 builds one bitmap over N[s] per root, |N[s]|^2/8 bytes, and a
#                 high-degree vertex is an endpoint in almost every batch (up
#                 to 1.05 GB on youtube and 6.4 GB on wikipedia inserting the
#                 whole stream, rebuilt every batch; here the graph holds about
#                 half the edges)
#   somei arms    SOMEI_TIMEOUT_SECONDS (default 2592000, thirty days: the
#                 deadline stops SOMEi first, with what is left of the
#                 allocation after edge-opt and star-rk, some 6.8 days of
#                 seven)
#   the others    RUN_TIMEOUT_SECONDS (default 21600, six hours: edge-opt ran
#                 out of two on coauthor and flickr in job 18429098; 0 = none)
#
# and no run is given more than the time left before DEADLINE_SECONDS, so the
# job ends, and writes its tables, inside its allocation.  A run killed by its
# limit is a timeout, which is a result.  The total of a timed-out star, edge or
# vertex run is at least
#
#   limit - outside
#
# where `outside` is the largest wall-minus-total of a finished star, edge or
# vertex run in the same cell: reading, the build and the graph updates, the
# same work in all of them.  SOMEi builds its own structure, so a SOMEi
# timeout has no bound; and only the total of a timeout is bounded, never its
# parts.
#
# SOMEI_DATASETS: the datasets SOMEi arms run on, set by PART (BUDGET below);
# "all" for every dataset.  On the others its arms are not run (excluded; no
# statistic).  Every mixed run first builds its window, half of its edges,
# untimed but inside its limit, and SOMEi does it in one batch
# (--somei-build=one-batch): one construction.  In the run's own batches it
# would be an insertion run of half the window, spent before a single mixed
# batch.  The ICTree is the same either way.  A timeout skips the smaller
# batches, so a large graph costs SOMEi three hours if it cannot finish the
# largest batch.
#
# BATCH ORDER.  Batch sizes run from the largest down.  The vertex search and
# SOMEi slow down as batches shrink, so with SKIP_SMALLER_AFTER_TIMEOUT=1
# (default) a vertex or SOMEi arm that timed out on a dataset is not run there
# again at that or any smaller batch size (skipped_after_timeout; no
# statistic).  The optimized split gets slower as batches grow, so the other
# arms are never skipped.
#
# OUT OF MEMORY.  A run that exhausts memory (std::bad_alloc, SLURM's
# oom-kill, or a kill signal before its limit) is recorded as out_of_memory: a
# result like a timeout, not a failure of the job, and no reason to skip
# anything.
#
# CORRECTNESS.  No certification pass: the mixed batches were checked against
# brute force, and between every search, subsumed search, strategy and SOMEi,
# when they were added (experiments/mixed_20261010; the sliding window's batch
# cutting against a change stream worked out by hand in make check).  Inside
# one cell and repeat every arm must report the same numbers of inserted and
# deleted edges and of remaining cliques; every arm that lists changes the
# same numbers of new and removed cliques; every search arm the same
# destroyed, recovered, new and subsumed cliques; every two-phase arm the same
# transient cliques.  And every run must have run what its arm says: search,
# subsumed search, size histogram, branch scope, strategy (a single-pass run
# nets nothing), the workload with its labels or its window (half of the edges
# the run used), SOMEi switches, batched full maintenance and tie order.
#
# GRID AND SCALE, the full grid (2026-10-15: the end-to-end jobs' 33 cells,
# shared by the full-grid insertion, deletion and mixed jobs):
#
#   dblp-cite, facebook          SMALL_BATCH_SIZES  1000 5000 10000 20000
#   coauthor, sx-stackoverflow,  LARGE_BATCH_SIZES  10000 20000 50000 100000
#   youtube, wikipedia, flickr                      200000
#
# HEAVY_DATASETS (youtube, wikipedia, flickr), which the one- and two-size
# jobs ran at 20000 only, run HEAVY_BATCH_SIZES, here every large size too.
# BATCH_SIZES, if set, replaces all three lists.  flickr and sx-stackoverflow
# run their first 1M edges, wikipedia its first 12M, the rest whole: a
# sliding-window run is cut there by --max-edges, its last batch possibly
# short, and a batch size larger than its window (half the budget) is skipped
# (on a whole stream the driver refuses such a run, and it fails; the
# smallest window, dblp-cite's 24,818 edges, is above every small size); a
# random-labels window is cut into whole batches by --max-batches, so it must
# be divisible by every batch size (a size that does not divide it is
# skipped).  The datasets run in DATASET_FILTER order, which PART sets
# (BUDGET below).  A dataset name that matches no file stops the job before
# anything runs.
#
# WARM-UP (random labels only, off by default): <DATASET>_WARMUP_EDGES > 0
# inserts that many edges first, in WARMUP_BATCH_SIZE batches, as part of the
# build; they take no label and are never deleted.  Only the window after them
# (target minus warm-up) is labelled and timed.  A sliding window builds its
# first window itself, so with WORKLOAD=sliding-window every warm-up must be 0.
#
# TIE ORDER.  --tie-order=file: edges go in by increasing timestamp, those
# sharing a timestamp in the dataset's line order -- the only order the code
# has.  Every run checks the order its statistics report.
#
# ARM ORDER.  Rotated one position from each cell to the next and from each
# repeat to the next.  REPEATS=1 by default; with more, the repeats are the
# outer loop.
#
# NODE.  Give the job its node to itself and do not run another timing job on
# the same node meanwhile: 18143576 and 18145683 overlapped on bryan1 for 3.6
# hours, and the first one's youtube cells came out about 7% slower for it;
# the mixed job 18429098 and the deletion job 18429099 shared bryan4 for all
# of their 16.7 hours.  OzSTAR refuses --exclusive, and other jobs share the
# node, so the three arms of the cell run back to back in one job and their
# ratios are formed there; report one node's times, since across the bryan
# nodes the ratios shifted by dataset (17899068).
# noise.tsv shows the build time over the star, edge and vertex runs of each
# cell, the same insertion in each.
#
# BUDGET, from the VI-B hybrid jobs on bryan8 (2026-10-28/29).  On coauthor
# edge-opt takes 2.2-3.1 hours a run (3.1 at batch 200000) and star-rk 5-8
# minutes, the three jobs within 5% of each other.  SOMEi took 18.2 hours for
# the first 12 of its 21 batches at 200000 (seconds each for the first six,
# then 5 minutes, 0.9, 0.8, 0.8, 7.2 and 8.4 hours), and its twenty hours ran
# out 1.7 hours into the 13th (19365025).  If the other nine take what the
# last two did, it needs some 3.5-4 days there, more if they keep growing.
# Inserting and deleting, its time at 100000 was 1.5 times that at 200000,
# at 50000 2.2-2.3 times and at 10000 5.2-5.4 times, so a week may hold
# 100000 but hardly the smaller batches.
#
# PART=coauthor-long, the only one: edge-opt, star-rk and somei on coauthor
# at one batch size (BATCH_SIZES, default 200000), in that order, so that the
# two short arms run first and SOMEi has the rest of the allocation.  One job
# a batch size, side by side:
#
#   sbatch --export=ALL,BATCH_SIZES=200000 --nodelist=<node> \
#          submit_overall_hybrid_coauthor_long_times.sh
#   sbatch --export=ALL,BATCH_SIZES=100000 --nodelist=<node> \
#          submit_overall_hybrid_coauthor_long_times.sh
#
# (A variable set in front of sbatch does not reach the job on OzSTAR.)  For
# an allocation other than seven days add --time=<days>-<hh>:<mm>:<ss>,
# longer than 24 hours.
#
# Random labels take less.
#
# BUILD.  Unless EXECUTABLE is given, the job builds its own binary from
# ${PROJECT_DIR}/src with the project Makefile, into a new directory under its
# output directory, and refuses a src/ without mixed batches (added
# 2026-10-10) or from before the 2026-10-09 cleanup.  Upload src/, Makefile and
# this script; load a C++20 compiler before sbatch or pass MODULES="gcc/13.3.0".
# Create the directory of the --output/--error paths above first: SLURM does
# not.
#
# OUTPUTS (in OUTPUT_DIR):
#   table_time.tsv     per cell: every arm's total, new-clique and subsumed
#                      times (median over repeats; ">=" a lower bound from a
#                      timeout)
#   table_speedup.tsv  per cell: every comparison's speedup (">=" a bound)
#   table_tree.tsv     per cell: the branches of every search and the branch
#                      ratios of the new_clique comparisons
#   mixed.tsv          per cell and arm: the inserted and deleted edges; all
#                      the new and removed cliques and the remaining ones; the
#                      destroyed (D), recovered (A), new (N), subsumed (S) and
#                      transient cliques (all_new = new + recovered,
#                      all_removed = destroyed + subsumed), the transient
#                      share of the net change; the median times by side --
#                      the search on the deleted and on the inserted edges, the
#                      subsumed search, the recovery, the index update, the
#                      netting -- and the total
#   times.tsv          per cell and arm: what its runs did, every run's total,
#                      the median of every part, the build and graph-update
#                      times, SOMEi's phases and the wall time
#   tree.tsv           per cell and search arm: the search statistics above,
#                      the branch ratio to METHOD_ARM, the branches by side,
#                      whether the runs agreed
#   recovery.tsv       per cell and search arm: what the recovery did --
#                      candidates, those under two vertices, those already
#                      accepted, maximality checks, rejected, accepted (in a
#                      two-phase run the transient cliques too, netted out
#                      later), the reads of the checks, and for rk its
#                      grouped checks (groups on bit rows, the candidates
#                      checked on them, the entries the rows read; 0 for
#                      the other subsumed searches) -- and its time
#   results.tsv        per cell and comparison: status, both times, speedup (or
#                      its bound), and for new_clique the branches and ratio
#   geomeans.tsv       per comparison, per batch size within the small and the
#                      large datasets, per group and overall: the geometric
#                      mean of the speedups over the finished cells, a lower
#                      bound with the timed-out cells added, and of the branch
#                      ratios
#   noise.tsv          per cell and repeat: the build and the graph-update
#                      times over the star, edge and vertex runs, the same work
#                      in each; a spread well above 1 marks a disturbed run
#   timeouts.tsv, summary.tsv (one row per run), manifest.txt, datasets.txt,
#   source_sha256.txt, build.log, probe.log

set -uo pipefail

SCRIPT_DIR="${PROJECT_DIR:-/home/zzhong/zzy_project/DynamicClique/IMMCE_SubsumedClique}"
PROJECT_DIR="${PROJECT_DIR:-$SCRIPT_DIR}"
DATA_PATH="${INPUT:-/fred/oz128/myDataset/DynamicGraph/Group1}"
# PART picks the arms, the datasets, SOMEi's datasets and limit, and the
# large graphs' batch sizes when they are not given (BUDGET above).
PART="${PART:-coauthor-long}"
case "$PART" in
    coauthor-long)
        PART_ARMS="edge-opt star-rk somei"
        PART_FILTER="processed_dblp_coauthor.txt"
        PART_SOMEI_DATASETS="all"
        PART_SOMEI_TIMEOUT="2592000"
        PART_LARGE_BATCH_SIZES="200000"
        ;;
    *)
        echo "Error: PART must be coauthor-long: $PART" >&2
        echo "Error: PART must be coauthor-long: $PART"
        exit 2
        ;;
esac
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DIR}/output/overall_hybrid_${PART}_${SLURM_JOB_ID:-local}}"
EXECUTABLE_EXPLICIT=0
[[ -n "${EXECUTABLE:-}" ]] && EXECUTABLE_EXPLICIT=1
EXECUTABLE="${EXECUTABLE:-}"
BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-1}}"
MODULES="${MODULES:-}"

REPEATS="${REPEATS:-1}"
SMALL_BATCH_SIZES="${SMALL_BATCH_SIZES:-1000 5000 10000 20000}"
LARGE_BATCH_SIZES="${LARGE_BATCH_SIZES:-$PART_LARGE_BATCH_SIZES}"
LARGE_DATASETS="${LARGE_DATASETS:-processed_dblp_coauthor.txt processed_flickr-growth.txt processed_wikipedia-growth.txt processed_youtube-u-growth.txt processed_sx-stackoverflow.txt}"
# Large datasets that run HEAVY_BATCH_SIZES instead of LARGE_BATCH_SIZES; in
# the full grid they are the same sizes.
HEAVY_DATASETS="${HEAVY_DATASETS:-processed_youtube-u-growth.txt processed_wikipedia-growth.txt processed_flickr-growth.txt}"
HEAVY_BATCH_SIZES="${HEAVY_BATCH_SIZES:-$LARGE_BATCH_SIZES}"
# Set, it replaces both lists above for every dataset.
BATCH_SIZES="${BATCH_SIZES:-}"
ARMS="${ARMS:-$PART_ARMS}"
# The arm every other one is compared with; it must be in ARMS.
METHOD_ARM="${METHOD_ARM:-star-rk}"
# Space-separated name=from:to:part entries.  Empty means the comparisons
# described above, built from the arms given.
EXPERIMENTS="${EXPERIMENTS:-}"
# sliding-window or random-labels (see THE WORKLOAD above).
WORKLOAD="${WORKLOAD:-sliding-window}"
# Random labels: the chance that an edge of the window is inserted rather
# than deleted, and the seed of the generator that draws them.
INSERT_FRACTION="${INSERT_FRACTION:-0.5}"
LABEL_SEED="${LABEL_SEED:-1}"
TIE_ORDER="${TIE_ORDER:-file}"
# Dataset file names, run in this order; "all" means every processed_*.txt
# found.  Unset, PART picks them, coauthor last: here it is the most
# expensive by far.
DATASET_FILTER="${DATASET_FILTER:-$PART_FILTER}"
# The datasets SOMEi arms run on; "all" means every one.
SOMEI_DATASETS="${SOMEI_DATASETS:-$PART_SOMEI_DATASETS}"
# Wall-time limits of one run, the build included: the vertex arms', the
# SOMEi arms', and every other arm's (0 = none for those).
VERTEX_TIMEOUT_SECONDS="${VERTEX_TIMEOUT_SECONDS:-3600}"
SOMEI_TIMEOUT_SECONDS="${SOMEI_TIMEOUT_SECONDS:-$PART_SOMEI_TIMEOUT}"
# Wall seconds a SOMEi run gets on top of SOMEI_TIMEOUT_SECONDS for its
# reading and untimed build (2026-11-06): SOMEI_TIMEOUT_SECONDS bounds its
# measured batches alone (--measured-time-limit), and this the rest.
SOMEI_SETUP_SECONDS="${SOMEI_SETUP_SECONDS:-3600}"
# A SOMEi run whose measured batches look set to finish within
# SOMEI_EXTEND_SECONDS -- the time so far plus the last batch's time for every
# batch left, judged after each batch -- runs on until then instead of
# stopping at SOMEI_TIMEOUT_SECONDS (2026-11-06, --measured-time-extend), so
# that a cell a little over the limit still gets its numbers.  Default
# SOMEI_TIMEOUT_SECONDS, that is off: since 2026-11-06 the limit is a flat
# five hours in every workload.
SOMEI_EXTEND_SECONDS="${SOMEI_EXTEND_SECONDS:-$SOMEI_TIMEOUT_SECONDS}"
RUN_TIMEOUT_SECONDS="${RUN_TIMEOUT_SECONDS:-21600}"
SKIP_SMALLER_AFTER_TIMEOUT="${SKIP_SMALLER_AFTER_TIMEOUT:-1}"
# Half an hour before the end of the job's allocation, which squeue gives
# (seven days unless sbatch --time gave another).  It has to be longer than
# 24 hours: a shorter one would repeat the twenty-hour run.  Given, as it has
# to be outside SLURM, DEADLINE_SECONDS stands.  No run gets more than the
# time left before it, and none starts with less than MIN_RUN_SECONDS left.
if [[ -z "${DEADLINE_SECONDS:-}" ]]; then
    allocation="$(squeue -h -j "${SLURM_JOB_ID:-none}" -o %l 2> /dev/null)"
    if [[ ! "$allocation" =~ ^(([0-9]+)-)?([0-9]+):([0-9]+):([0-9]+)$ ]]; then
        echo "Error: the deadline comes from the job's allocation, which squeue did not give ('$allocation'); give DEADLINE_SECONDS" >&2
        echo "Error: the deadline comes from the job's allocation, which squeue did not give ('$allocation'); give DEADLINE_SECONDS"
        exit 2
    fi
    days=$(( 10#${BASH_REMATCH[2]:-0} ))
    allocation=$(( ((days * 24 + 10#${BASH_REMATCH[3]}) * 60 + 10#${BASH_REMATCH[4]}) * 60 + 10#${BASH_REMATCH[5]} ))
    if (( allocation <= 86400 )); then
        echo "Error: the allocation is ${allocation}s, not longer than 24 hours: leave the script's seven days or give --time=<days>-<hh>:<mm>:<ss>" >&2
        echo "Error: the allocation is ${allocation}s, not longer than 24 hours: leave the script's seven days or give --time=<days>-<hh>:<mm>:<ss>"
        exit 2
    fi
    DEADLINE_SECONDS=$(( allocation - 1800 ))
fi
MIN_RUN_SECONDS="${MIN_RUN_SECONDS:-300}"

# Per-dataset edge budgets, the end-to-end job's.  0 means the whole stream.
# With random labels the window (target minus warm-up) must be divisible by
# every batch size, so every batch size labels and replays exactly the same
# edges; a sliding window cuts the run at the budget itself.
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
WARMUP_BATCH_SIZE="${WARMUP_BATCH_SIZE:-20000}"

for list_name in DATASET_FILTER LARGE_DATASETS HEAVY_DATASETS SMALL_BATCH_SIZES \
                 LARGE_BATCH_SIZES HEAVY_BATCH_SIZES BATCH_SIZES ARMS EXPERIMENTS \
                 SOMEI_DATASETS; do
    printf -v "$list_name" '%s' "${!list_name//,/ }"
done

# Under SLURM, messages from the checks below would otherwise land only in the
# .err file and leave the .out file empty, so a refusal prints on both.
note() {
    echo "$*" >&2
    echo "$*"
}
trap 'status=$?; (( status != 0 )) && echo "exited with status $status"' EXIT

echo "=== immce overall comparison (paper Section VI-B), hybrid batches, part $PART: total, new-clique and subsumed times ==="
echo "Job ID:         ${SLURM_JOB_ID:-local}"
echo "Host:           $(hostname)"
echo "Started:        $(date -Is)"
echo "Project dir:    $PROJECT_DIR"
echo "Data path:      $DATA_PATH"
echo "Dataset filter: $DATASET_FILTER"
echo "SOMEi datasets: $SOMEI_DATASETS"
echo "Output dir:     $OUTPUT_DIR"
echo "Executable:     ${EXECUTABLE:-(built by this job)}"
echo "Arms:           $ARMS"
echo "Method arm:     $METHOD_ARM"
if [[ "$WORKLOAD" == random-labels ]]; then
    echo "Workload:       random-labels, insert with probability $INSERT_FRACTION, seed $LABEL_SEED"
else
    echo "Workload:       $WORKLOAD, the window half of each run's edges"
fi
echo "Tie order:      $TIE_ORDER"
if [[ -n "$BATCH_SIZES" ]]; then
    echo "Batch sizes:    $BATCH_SIZES (every dataset), largest first"
else
    echo "Batch sizes:    $SMALL_BATCH_SIZES (small datasets), $LARGE_BATCH_SIZES (large), $HEAVY_BATCH_SIZES (heavy), largest first"
    echo "Large datasets: $LARGE_DATASETS"
    echo "Heavy datasets: $HEAVY_DATASETS"
fi
echo "Repeats:        $REPEATS"
echo "Time limits:    vertex ${VERTEX_TIMEOUT_SECONDS}s, somei ${SOMEI_TIMEOUT_SECONDS}s of measured batches, up to ${SOMEI_EXTEND_SECONDS}s for one set to finish by then (${SOMEI_SETUP_SECONDS}s more for its reading and build), other arms ${RUN_TIMEOUT_SECONDS}s (skip smaller batches after a vertex or SOMEi timeout: $SKIP_SMALLER_AFTER_TIMEOUT)"
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

# Whether a dataset (file name) is in HEAVY_DATASETS.
dataset_is_heavy() {
    local heavy
    for heavy in $HEAVY_DATASETS; do
        [[ "$1" == "$heavy" ]] && return 0
    done
    return 1
}

# The batch sizes one dataset (file name) runs, largest first.
batch_sizes_for() {
    local sizes
    if [[ -n "$BATCH_SIZES" ]]; then
        sizes="$BATCH_SIZES"
    elif dataset_is_heavy "$1"; then
        sizes="$HEAVY_BATCH_SIZES"
    elif [[ "$(dataset_class_for "$1")" == large ]]; then
        sizes="$LARGE_BATCH_SIZES"
    else
        sizes="$SMALL_BATCH_SIZES"
    fi
    tr ' ' '\n' <<< "$sizes" | grep -v '^$' | sort -rn -u | tr '\n' ' '
}

# Whether SOMEi arms run on a dataset (file name).
somei_runs_on() {
    local name
    [[ "$SOMEI_DATASETS" == all ]] && return 0
    for name in $SOMEI_DATASETS; do
        [[ "$1" == "$name" ]] && return 0
    done
    return 1
}

# Parses an arm into ARM_KIND (search or somei) and its settings:
# ARM_FRAMEWORK (star / edge / vertex / somei), ARM_SCOPE (a star arm's branch
# scope, - otherwise), ARM_ALGORITHM (the subsumed search, none for SOMEi),
# ARM_GATE (on when the index keeps the size histogram the size gate reads,
# off when not, - for SOMEi), ARM_STRATEGY (single-pass, two-phase or somei),
# ARM_PRUNING and ARM_CHANGES (SOMEi's switches, - otherwise), ARM_LIMIT (its
# time limit in seconds, 0 for none) and ARM_SKIPS (1 if it skips the smaller
# batch sizes after a timeout).
set_arm() {
    local rest
    ARM_KIND=search
    ARM_SCOPE=-
    ARM_ALGORITHM=none
    ARM_GATE=-
    ARM_STRATEGY=two-phase
    ARM_PRUNING=-
    ARM_CHANGES=-
    ARM_LIMIT="$RUN_TIMEOUT_SECONDS"
    ARM_SKIPS=0
    case "$1" in
        somei|somei-*)
            ARM_KIND=somei
            ARM_FRAMEWORK=somei
            ARM_STRATEGY=somei
            ARM_PRUNING=on
            ARM_CHANGES=on
            ARM_LIMIT=$(( SOMEI_EXTEND_SECONDS + SOMEI_SETUP_SECONDS ))
            ARM_SKIPS=1
            rest="${1#somei}"
            if [[ "$rest" == -noprune* ]]; then
                ARM_PRUNING=off
                rest="${rest#-noprune}"
            fi
            if [[ "$rest" == -nochg ]]; then
                ARM_CHANGES=off
                rest=""
            fi
            [[ -z "$rest" ]]
            return
            ;;
        star-root-*) ARM_FRAMEWORK=star; ARM_SCOPE=root; rest="${1#star-root-}" ;;
        star-*) ARM_FRAMEWORK=star; ARM_SCOPE=unbound; rest="${1#star-}" ;;
        edge-*) ARM_FRAMEWORK=edge; rest="${1#edge-}" ;;
        vertex-*)
            ARM_FRAMEWORK=vertex
            ARM_LIMIT="$VERTEX_TIMEOUT_SECONDS"
            ARM_SKIPS=1
            rest="${1#vertex-}"
            ;;
        *) return 1 ;;
    esac
    if [[ "$rest" == *-singlepass ]]; then
        ARM_STRATEGY=single-pass
        rest="${rest%-singlepass}"
    fi
    case "$rest" in
        rk) ARM_ALGORITHM=residual-kernel; ARM_GATE=on ;;
        rk-nogate) ARM_ALGORITHM=residual-kernel; ARM_GATE=off ;;
        optbm) ARM_ALGORITHM=optimized-bitmap; ARM_GATE=off ;;
        opt) ARM_ALGORITHM=optimized; ARM_GATE=off ;;
        *) return 1 ;;
    esac
}
ARM_FORM="<star|star-root|edge|vertex>-<rk|rk-nogate|optbm|opt>[-singlepass] or somei[-noprune][-nochg]"

for value in "$REPEATS" "$VERTEX_TIMEOUT_SECONDS" "$SOMEI_TIMEOUT_SECONDS" \
             "$WARMUP_BATCH_SIZE" "$MIN_RUN_SECONDS"; do
    if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
        note "Error: REPEATS, VERTEX_TIMEOUT_SECONDS, SOMEI_TIMEOUT_SECONDS, WARMUP_BATCH_SIZE and MIN_RUN_SECONDS must be positive integers: $value"
        exit 2
    fi
done
if [[ ! "$SOMEI_SETUP_SECONDS" =~ ^[0-9]+$ ]]; then
    note "Error: SOMEI_SETUP_SECONDS must be a nonnegative integer: $SOMEI_SETUP_SECONDS"
    exit 2
fi
if [[ ! "$SOMEI_EXTEND_SECONDS" =~ ^[1-9][0-9]*$ ]] || \
        (( SOMEI_EXTEND_SECONDS < SOMEI_TIMEOUT_SECONDS )); then
    note "Error: SOMEI_EXTEND_SECONDS must be an integer no smaller than SOMEI_TIMEOUT_SECONDS: $SOMEI_EXTEND_SECONDS"
    exit 2
fi
if [[ "$SKIP_SMALLER_AFTER_TIMEOUT" != 0 && "$SKIP_SMALLER_AFTER_TIMEOUT" != 1 ]]; then
    note "Error: SKIP_SMALLER_AFTER_TIMEOUT must be 0 or 1: $SKIP_SMALLER_AFTER_TIMEOUT"
    exit 2
fi
# One cell a job, so that the arms run in the order given and SOMEi, last,
# has the rest of the allocation.
read -r -a job_datasets <<< "$DATASET_FILTER"
job_sizes=()
if (( ${#job_datasets[@]} == 1 )) && [[ "${job_datasets[0]}" != all ]]; then
    read -r -a job_sizes <<< "$(batch_sizes_for "${job_datasets[0]}")"
fi
if (( ${#job_sizes[@]} != 1 )); then
    note "Error: a job runs one dataset at one batch size (DATASET_FILTER, BATCH_SIZES): $DATASET_FILTER at ${BATCH_SIZES:-$LARGE_BATCH_SIZES}"
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
# Every run prints its labels back -- the fraction with six significant
# digits, the seed as an integer -- and is checked against these, so both take
# their plain forms.
if [[ ! "$INSERT_FRACTION" =~ ^[01]([.][0-9]+)?$ ]] || \
   ! awk -v p="$INSERT_FRACTION" 'BEGIN {exit !(p + 0 <= 1)}'; then
    note "Error: INSERT_FRACTION must be a decimal from 0 to 1, such as 0.5: $INSERT_FRACTION"
    exit 2
fi
if [[ ! "$LABEL_SEED" =~ ^(0|[1-9][0-9]{0,18})$ ]]; then
    note "Error: LABEL_SEED must be a nonnegative integer of at most 19 digits without leading zeros: $LABEL_SEED"
    exit 2
fi
case "$WORKLOAD" in
    sliding-window)
        for value in "$FLICKR_WARMUP_EDGES" "$WIKIPEDIA_WARMUP_EDGES" "$YOUTUBE_WARMUP_EDGES" \
                     "$STACKOVERFLOW_WARMUP_EDGES" "$DEFAULT_WARMUP_EDGES"; do
            if (( value != 0 )); then
                note "Error: a sliding window builds its first window itself, so every <DATASET>_WARMUP_EDGES must be 0 with WORKLOAD=sliding-window: $value"
                exit 2
            fi
        done
        ;;
    random-labels)
        if awk -v p="$INSERT_FRACTION" 'BEGIN {exit !(p + 0 == 0 || p + 0 == 1)}'; then
            echo "Note: INSERT_FRACTION=$INSERT_FRACTION, so every batch only $( [[ "${INSERT_FRACTION:0:1}" == 1 ]] && echo inserts || echo deletes ) and the two strategies coincide."
        fi
        ;;
    *)
        note "Error: WORKLOAD must be sliding-window or random-labels: $WORKLOAD"
        exit 2
        ;;
esac
# Edges sharing a timestamp keep the input file's order: the only order the
# reader has since the vertex and shuffle orders were removed (2026-10-08).
if [[ "$TIE_ORDER" != file ]]; then
    note "Error: TIE_ORDER must be file (the vertex and shuffle orders were removed): $TIE_ORDER"
    exit 2
fi
for list_name in SMALL_BATCH_SIZES LARGE_BATCH_SIZES HEAVY_BATCH_SIZES; do
    if [[ -z "${!list_name// /}" ]]; then
        note "Error: $list_name must contain at least one batch size"
        exit 2
    fi
done
for batch_size in $SMALL_BATCH_SIZES $LARGE_BATCH_SIZES $HEAVY_BATCH_SIZES $BATCH_SIZES; do
    if [[ ! "$batch_size" =~ ^[1-9][0-9]*$ ]]; then
        note "Error: invalid batch size: $batch_size"
        exit 2
    fi
done

read -r -a arms <<< "$ARMS"
if (( ${#arms[@]} < 2 )); then
    note "Error: ARMS must name at least two arms to compare: $ARMS"
    exit 2
fi
declare -A arm_seen=()
declare -A arm_kind=()
declare -A arm_strategy=()
for arm in "${arms[@]}"; do
    if ! set_arm "$arm"; then
        note "Error: unknown arm: $arm (arms have the form $ARM_FORM)"
        exit 2
    fi
    if [[ -n "${arm_seen[$arm]:-}" ]]; then
        note "Error: arm listed twice: $arm"
        exit 2
    fi
    arm_seen[$arm]=1
    arm_kind[$arm]="$ARM_KIND"
    arm_strategy[$arm]="$ARM_STRATEGY"
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
# An arm name as it reads inside a comparison name.
arm_tag() {
    printf '%s' "${1//-/_}"
}
# Adds name=from:to:part when both arms are run.
add_experiment() {
    if known_arm "$2" && known_arm "$3" && [[ "$2" != "$3" ]]; then
        EXPERIMENTS+=" $1=$2:$3:$4"
    fi
}
# The same search and subsumed search with the other strategy.
twin_of() {
    if [[ "$1" == *-singlepass ]]; then
        printf '%s' "${1%-singlepass}"
    else
        printf '%s' "$1-singlepass"
    fi
}
# The edge arm a vertex arm's search is compared with: the one with the same
# subsumed search and strategy if there is one, else one with the same
# strategy, else any (the new-clique search depends on neither).
edge_reference_for() {
    local same="" other="-singlepass" candidate
    if [[ "$1" == *-singlepass ]]; then
        same="-singlepass"
        other=""
    fi
    for candidate in "edge-$1" edge-{rk,rk-nogate,optbm,opt}"$same" \
                     edge-{rk,rk-nogate,optbm,opt}"$other"; do
        if known_arm "$candidate"; then
            echo "$candidate"
            return
        fi
    done
}
if [[ -z "$EXPERIMENTS" ]]; then
    for arm in "${arms[@]}"; do
        [[ "$arm" == "$METHOD_ARM" ]] && continue
        base="$(arm_tag "$arm")_vs_$(arm_tag "$METHOD_ARM")"
        add_experiment "${base}_total" "$arm" "$METHOD_ARM" total
        if [[ "${arm_kind[$arm]}" == search && "${arm_kind[$METHOD_ARM]}" == search ]]; then
            add_experiment "${base}_new_clique" "$arm" "$METHOD_ARM" new_clique
            add_experiment "${base}_subsumed" "$arm" "$METHOD_ARM" subsumed
        fi
    done
    # A search arm whose strategy is not the method's, against the method's
    # searches with that strategy, when that arm runs: the searches alone.
    method_twin="$(twin_of "$METHOD_ARM")"
    if [[ "${arm_kind[$METHOD_ARM]}" == search ]] && known_arm "$method_twin"; then
        for arm in "${arms[@]}"; do
            [[ "${arm_kind[$arm]}" == search && "$arm" != "$METHOD_ARM" && \
               "$arm" != "$method_twin" && \
               "${arm_strategy[$arm]}" == "${arm_strategy[$method_twin]}" ]] || continue
            base="$(arm_tag "$arm")_vs_$(arm_tag "$method_twin")"
            add_experiment "${base}_total" "$arm" "$method_twin" total
            add_experiment "${base}_subsumed" "$arm" "$method_twin" subsumed
        done
    fi
    # The two strategies of one search and subsumed search, unless one of them
    # is the method, whose comparisons cover that pair.
    for arm in "${arms[@]}"; do
        [[ "$arm" == *-singlepass ]] || continue
        twin="$(twin_of "$arm")"
        known_arm "$twin" || continue
        [[ "$arm" == "$METHOD_ARM" || "$twin" == "$METHOD_ARM" ]] && continue
        base="$(arm_tag "$arm")_vs_$(arm_tag "$twin")"
        add_experiment "${base}_total" "$arm" "$twin" total
        add_experiment "${base}_subsumed" "$arm" "$twin" subsumed
    done
    for arm in "${arms[@]}"; do
        [[ "$arm" == vertex-* ]] || continue
        reference="$(edge_reference_for "${arm#vertex-}")"
        [[ -z "$reference" || "$reference" == "$METHOD_ARM" || "$arm" == "$METHOD_ARM" ]] && continue
        add_experiment "$(arm_tag "$arm")_vs_$(arm_tag "$reference")_new_clique" \
            "$arm" "$reference" new_clique
    done
    EXPERIMENTS="${EXPERIMENTS# }"
fi
if [[ -z "$EXPERIMENTS" ]]; then
    note "Error: the arms give no comparison: $ARMS"
    exit 2
fi
for experiment in $EXPERIMENTS; do
    name="${experiment%%=*}"
    IFS=: read -r from to part extra <<< "${experiment#*=}"
    if [[ "$experiment" != *=* || -z "$name" || -z "$from" || -z "$to" || \
          -z "$part" || -n "$extra" ]]; then
        note "Error: experiment '$experiment' is not name=from:to:part"
        exit 2
    fi
    if ! known_arm "$from" || ! known_arm "$to"; then
        note "Error: experiment '$experiment' names an arm that is not in ARMS"
        exit 2
    fi
    case "$part" in
        total) ;;
        new_clique|subsumed|index_update)
            if [[ "${arm_kind[$from]}" != search || "${arm_kind[$to]}" != search ]]; then
                note "Error: experiment '$experiment': only star, edge and vertex arms have a $part part"
                exit 2
            fi
            ;;
        *)
            note "Error: experiment '$experiment': the part must be total, new_clique, subsumed or index_update"
            exit 2
            ;;
    esac
done
arm_kind_list=""
arm_strategy_list=""
for arm in "${arms[@]}"; do
    arm_kind_list+=" $arm:${arm_kind[$arm]}"
    arm_strategy_list+=" $arm:${arm_strategy[$arm]}"
done
arm_kind_list="${arm_kind_list# }"
arm_strategy_list="${arm_strategy_list# }"
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
    note "Error: no timeout on PATH (coreutils); every run needs it for its time limit."
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
        "src/maintenance_engine.hpp|enum class MixedStrategy|it has no mixed batches (added 2026-10-10)" \
        "src/main.cpp|UpdateModeOption::SlidingWindow|its driver has no --update-mode=sliding-window (added 2026-10-10)" \
        "src/main.cpp|\"--max-edges=\"|its driver cannot cut the stream by edge count (added 2026-10-10)" \
        "src/main.cpp|UpdateModeOption::RandomLabels|its driver has no --update-mode=random-labels (added 2026-10-10)" \
        "src/main.cpp|\"mixed: strategy=\"|its driver reports no mixed statistics (added 2026-10-10)" \
        "src/main.cpp|\"recovery: destroyed=\"|its driver reports no recovery statistics (added 2026-10-10)" \
        "src/incremental_graph.hpp|void removeBatch|its graph cannot delete edges (added 2026-10-10)" \
        "src/subsumed_clique_search.hpp|recoverBatch|its subsumed search cannot recover the cliques a deletion creates (added 2026-10-10)" \
        "src/subsumed_clique_search.hpp|memberOrder_|its recovery finds every vertex extending a candidate before deciding (the early exit came on 2026-10-10)" \
        "src/subsumed_clique_search.cpp|minRowGroup|its residual-kernel recovery checks every candidate one by one (grouping by shared member came on 2026-10-11)" \
        "src/subsumed_clique_search.hpp|acceptedInBatch|its residual-kernel recovery probes the clique index for every candidate (the faster checks came on 2026-10-12)" \
        "src/subsumed_clique_search.cpp|gatherLimit|its residual-kernel recovery can gather more vertices than 32-bit offsets reach, which corrupted the clique index on coauthor batch 200000 (fixed 2026-10-12)" \
        "src/subsumed_clique_search.cpp|collectExtensionMasks|its residual-kernel recovery checks every candidate of a general kernel (the extension masks came on 2026-10-12)" \
        "src/subsumed_clique_search.cpp|maskFactor|its residual-kernel recovery puts every general kernel under extension masks, 6-15 times slower on flickr (fixed 2026-10-13)" \
        "src/subsumed_clique_search.hpp|radixCounts_|its residual-kernel recovery sorts the deferred candidates by comparison and ANDs whole rows (the faster pass came on 2026-10-13)" \
        "src/subsumed_clique_search.hpp|wideDegree_|its residual-kernel search sends cliques wider than 64 vertices with complete kernels through the complement search (the closed form came on 2026-10-13)" \
        "src/subsumed_clique_search.hpp|centreBits_|its residual-kernel recovery builds group rows by merging lists and numbers each centre in arrays over every vertex (the bit lookup came on 2026-10-13)" \
        "src/subsumed_clique_search.cpp|if (readMemberLists(candidate) >= static_cast<std::size_t>(deferDegree)) {|its opt and optbm recovery checks every candidate on the spot (since 2026-10-14 the baselines stop at the first extender and check the high-degree candidates grouped at the end of the batch, as residual-kernel does)" \
        "src/main.cpp|someiOneBatchBuild|its SOMEi builds the starting graph batch by batch and prints no progress lines (--somei-build and --progress came on 2026-10-27)" \
        "src/main.cpp|immceStopAtMeasuredTimeLimit|its driver has no --measured-time-limit, so a slow SOMEi build eats into its measured time (added 2026-11-06)" \
        "src/main.cpp|measuredTimeExtendSeconds|its driver cannot let a SOMEi run that is set to finish soon run past its limit (--measured-time-extend, added 2026-11-06)" \
        "src/somei.hpp|enumerateBuild|its SOMEi builds the starting graph by the paper's own construction, hours on youtube, wikipedia and flickr (the same state by Bron-Kerbosch on G[R], --somei-build-update, came on 2026-11-06)" \
        "src/somei.hpp|beginBatch|its SOMEi cannot delete edges (added 2026-10-10)" \
        "src/main.cpp|tie-order=file|its reader has no --tie-order=file (2026-10-05)" \
        "src/vertex_clique_search.hpp|class VertexCliqueEnumerator|it has no vertex baseline (added 2026-10-09)" \
        "src/main.cpp|FrameworkOption::Vertex|its driver cannot run the vertex baseline (added 2026-10-09)" \
        "src/main.cpp|\"edge: time-ms=\"|it is from before the 2026-10-09 cleanup (the edge line was still edge-fast:)" \
        "src/star_clique_search.hpp|StarBranchScope branchScope = StarBranchScope::Unbound;|it is from before the 2026-10-09 cleanup (star still had its root shapes and partitions as switches)" \
        "src/star_clique_search.hpp|std::deque<LocalVertexSet> deepMasks_;|its deep re-indexing keeps nested partner masks in a vector, which loses cliques" \
        "src/star_clique_search.hpp|branchUnitBindings|its unbound scope walks the star's whole partner list and scans a pivot for a lone partner (before 2026-10-08)" \
        "src/local_search_graph.hpp|moveToExcludedBits|its BK keeps sorted vertex lists on every sibling move" \
        "src/local_bk.hpp|candidateCount - 1|its pivot scan has no early exit on the P side" \
        "src/wide_residual.hpp|class DynamicSearch|its complement search splits kernels past 256 endpoints (before 2026-10-07)"; do
        IFS='|' read -r marker_file marker_text marker_reason <<< "$marker"
        if ! grep -qF -- "$marker_text" "${PROJECT_DIR}/${marker_file}" 2> /dev/null; then
            note "Error: ${PROJECT_DIR}/src is older than this script: $marker_reason."
            note "Upload the current src/ and Makefile."
            exit 2
        fi
    done
    for required_source in src/maintenance_engine.cpp src/vertex_clique_search.cpp src/somei.cpp; do
        if ! grep -qF "$required_source" "${PROJECT_DIR}/Makefile"; then
            note "Error: ${PROJECT_DIR}/Makefile does not build $required_source; upload the current Makefile."
            exit 2
        fi
    done
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
    somei-build somei-build-update progress measured-time-limit measured-time-extend \
    framework update-mode random-labels insert-fraction label-seed \
    sliding-window window-edges max-edges \
    mixed-strategy two-phase maintenance maintenance-pipeline \
    subsumed-algorithm optimized-bitmap residual-kernel subsumed-size-histogram \
    star-branch-scope somei-pruning somei-changes tie-order max-batches \
    warmup-edges warmup-batch-size no-output stats 'vertex|both|all|somei'; do
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
# would otherwise quietly drop a dataset from an hours-long job, run a large
# one at the small batch sizes, or let SOMEi run where it was not meant to.
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
    for wanted in $LARGE_DATASETS $HEAVY_DATASETS; do
        [[ -n "${dataset_present[$wanted]:-}" ]] || missing="${missing:+$missing }$wanted"
    done
    if [[ "$SOMEI_DATASETS" != all ]]; then
        for wanted in $SOMEI_DATASETS; do
            [[ -n "${dataset_present[$wanted]:-}" ]] || missing="${missing:+$missing }$wanted"
        done
    fi
fi
if [[ -n "$missing" ]]; then
    note "Error: DATASET_FILTER, LARGE_DATASETS, HEAVY_DATASETS or SOMEI_DATASETS names files that are not under $DATA_PATH: $missing"
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

LAUNCHER=()
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    LAUNCHER=(srun --exclusive)
fi

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
        --maintenance=full
        "--update-mode=${WORKLOAD}"
        --maintenance-pipeline=batched
        "--tie-order=${TIE_ORDER}"
    )
    # A sliding window keeps half of the run's edges, the driver's default.
    if [[ "$WORKLOAD" == random-labels ]]; then
        ARM_ARGUMENTS+=("--insert-fraction=${INSERT_FRACTION}" "--label-seed=${LABEL_SEED}")
    fi
    if [[ "$ARM_KIND" == somei ]]; then
        # The untimed build in one batch, its steps by Bron-Kerbosch on G[R]
        # (the paper's state, 2026-11-06), so that the limit falls on the
        # measured batches; a line per measured batch, so that a run the
        # limit stops shows how far it got (somei_progress.tsv); and the
        # limit on the measured batches alone (2026-11-06).
        ARM_ARGUMENTS+=("--somei-pruning=${ARM_PRUNING}" "--somei-changes=${ARM_CHANGES}"
                        --somei-build=one-batch --somei-build-update=bk --progress
                        "--measured-time-limit=${SOMEI_TIMEOUT_SECONDS}"
                        "--measured-time-extend=${SOMEI_EXTEND_SECONDS}")
    else
        ARM_ARGUMENTS+=("--mixed-strategy=${ARM_STRATEGY}" "--subsumed-algorithm=${ARM_ALGORITHM}")
        if [[ "$ARM_ALGORITHM" == residual-kernel && "$ARM_GATE" == off ]]; then
            ARM_ARGUMENTS+=(--subsumed-size-histogram=off)
        fi
        if [[ "$ARM_FRAMEWORK" == star ]]; then
            ARM_ARGUMENTS+=("--star-branch-scope=${ARM_SCOPE}")
        fi
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

# One batch per arm, on the first 20000 edges (a sliding window starts at
# half the stream, 2026-11-07), proves the parser takes every switch and a mixed batch
# completes before hours are spent on the sweep.
probe_binary() {
    local probe_log="${OUTPUT_DIR}/probe.log"
    local smallest arm_log="${OUTPUT_DIR}/probe.arm.log"
    smallest="$(smallest_dataset)"
    : > "$probe_log"
    echo "[$(date -Is)] probing $EXECUTABLE on $(dataset_tag_of "$smallest")"
    for arm in "${arms[@]}"; do
        arm_arguments "$arm"
        echo "== $arm: ${ARM_ARGUMENTS[*]}" >> "$probe_log"
        if ! timeout --kill-after=60 600 "$EXECUTABLE" "$smallest" "${ARM_ARGUMENTS[@]}" \
                --batch-size=1000 --max-batches=1 --max-edges=20000 --no-output --stats \
                > "$arm_log" 2>&1; then
            cat "$arm_log" >> "$probe_log"
            return 1
        fi
        cat "$arm_log" >> "$probe_log"
        grep -q "update-mode=${WORKLOAD}" "$arm_log" || return 1
        grep -q "tie-order=${TIE_ORDER}" "$arm_log" || return 1
        grep -q "^mixed: strategy=${ARM_STRATEGY}," "$arm_log" || return 1
        if [[ "$ARM_KIND" == somei ]]; then
            grep -q '^somei: time-ms=' "$arm_log" || return 1
        else
            grep -q "^${ARM_FRAMEWORK}: time-ms=" "$arm_log" || return 1
            grep -q '^recovery: destroyed=' "$arm_log" || return 1
        fi
    done
    rm -f "$arm_log"
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
MIXED_FILE="${OUTPUT_DIR}/mixed.tsv"
TREE_FILE="${OUTPUT_DIR}/tree.tsv"
RECOVERY_FILE="${OUTPUT_DIR}/recovery.tsv"
TABLE_TIME_FILE="${OUTPUT_DIR}/table_time.tsv"
TABLE_SPEEDUP_FILE="${OUTPUT_DIR}/table_speedup.tsv"
TABLE_TREE_FILE="${OUTPUT_DIR}/table_tree.tsv"
GEOMEAN_FILE="${OUTPUT_DIR}/geomeans.tsv"
NOISE_FILE="${OUTPUT_DIR}/noise.tsv"
TIMEOUT_FILE="${OUTPUT_DIR}/timeouts.tsv"
# The header and the row are written from one list and one array, and
# write_summary_row refuses a row whose length disagrees with it.
SUMMARY_COLUMNS=(
    dataset dataset_class batch_size warmup_edges max_batches target_edges
    arm kind framework branch_scope subsumed_algorithm size_histogram strategy
    somei_pruning somei_changes workload insert_fraction label_seed tie_order
    position repeat
    status exit_code wall_seconds time_limit_seconds
    total_ms new_clique_ms subsumed_ms index_update_ms netting_ms
    destroyed_search_ms new_search_ms subsumed_search_ms recovery_ms
    graph_update_ms build_ms input_preprocess_ms
    somei_scope_ms somei_rebuild_ms somei_change_ms
    clique_counts_match switches_check
    batches window_edges inserted_edges deleted_edges
    new_cliques removed_cliques current_maximal_cliques
    mixed_destroyed mixed_recovered mixed_new mixed_subsumed mixed_transient
    branches destroyed_side_branches new_side_branches
    search_roots branches_below_roots leaves emitted pivot_calls pivot_candidates
    local_builds max_local_vertices peak_local_bytes
    edge_seeds vertex_roots vertex_pruned_roots vertex_new_edge_prunes
    star_stars star_singleton_stars star_per_partner_star_edges star_pivot_root_stars
    star_partner_prunes
    recovery_candidates recovery_below_two recovery_duplicates recovery_checks
    recovery_rejected recovery_accepted recovery_neighbor_reads recovery_adjacency_tests
    recovery_groups recovery_grouped_checks recovery_row_reads
    size_gate_skipped_cliques size_gate_skipped_probes fallback_cliques kernel_dynamic
    somei_rebuilt_headers somei_rolled_back_bodies somei_redo_steps somei_body_comparisons
    stats_file
)
(IFS=$'\t'; printf '%s\n' "${SUMMARY_COLUMNS[*]}") > "$SUMMARY_FILE"
printf 'dataset\tbatch_size\tarm\trepeat\ttimeout_seconds\tcut_by_deadline\tstats_file\n' > "$TIMEOUT_FILE"

{
    echo "experiment=immce-subcliquesearch-overall-hybrid-coauthor-long"
    echo "part=$PART"
    echo "workload=$WORKLOAD"
    if [[ "$WORKLOAD" == random-labels ]]; then
        echo "dynamic_semantics=mixed batches: --update-mode=random-labels, every window edge labelled insert with probability $INSERT_FRACTION (seed $LABEL_SEED, stream order); the delete-labelled edges inserted untimed (the build: star + residual-kernel; SOMEi its own), then every batch inserting its insert-labelled and deleting its delete-labelled edges at once (timed)"
        echo "insert_fraction=$INSERT_FRACTION"
        echo "label_seed=$LABEL_SEED"
    else
        echo "dynamic_semantics=mixed batches: --update-mode=sliding-window, the graph keeps the newest W of the run's first T edges (T the edge budget, --max-edges, or the whole stream; W = T/2); the first W inserted untimed (the build: star + residual-kernel; SOMEi its own), then every batch inserting the next edges and deleting as many of the oldest at once (timed)"
    fi
    echo "strategies=two-phase (Das et al. section 4.5; the default, arms without -singlepass) and single-pass (ours; -singlepass arms); SOMEi updates each batch its own way"
    echo "project_dir=$PROJECT_DIR"
    echo "data_path=$DATA_PATH"
    echo "dataset_filter=$DATASET_FILTER"
    echo "somei_datasets=$SOMEI_DATASETS"
    echo "executable=$EXECUTABLE"
    echo "executable_sha256=$(sha256sum "$EXECUTABLE" | cut -d' ' -f1)"
    echo "built_by_job=$(( EXECUTABLE_EXPLICIT == 0 ? 1 : 0 ))"
    echo "build_dir=${BUILD_DIR:-(prebuilt)}"
    echo "build_cxxflags=$BUILD_CXXFLAGS"
    echo "compiler=$(g++ --version 2> /dev/null | head -n1)"
    echo "timing_scope=the mixed batches only; total = timing: search-and-maintenance-ms = new clique (<search>: time-ms, both runs: the deleted edges on the graph before the batch, finding the destroyed cliques, and the inserted ones on the graph after it, finding the new cliques) + subsumed (subsumed: time-ms, the subsumed search and the recovery) + index update (clique-index-update-ms) + netting (mixed: transient-netting-ms, two-phase only); somei: its own update and change listing, no new/subsumed split; graph-update-ms and the build (warmup-ms) outside"
    echo "branches=recursive-calls of the new-clique search over the batches, both runs, roots included (mixed: destroyed-side-recursive-calls and new-side-recursive-calls split them); search_roots: edge seeds, vertex roots - pruned-roots, star singleton-stars + per-partner-star-edges + branch-root-pivot-stars"
    echo "tie_order=$TIE_ORDER"
    echo "method_arm=$METHOD_ARM"
    echo "dataset_count=${#DATASETS[@]}"
    echo "arms=$ARMS"
    echo "experiments=$EXPERIMENTS"
    echo "batch_sizes_override=${BATCH_SIZES:-(none)}"
    echo "small_batch_sizes=$SMALL_BATCH_SIZES"
    echo "large_batch_sizes=$LARGE_BATCH_SIZES"
    echo "large_datasets=$LARGE_DATASETS"
    echo "heavy_batch_sizes=$HEAVY_BATCH_SIZES"
    echo "heavy_datasets=$HEAVY_DATASETS"
    echo "batch_order=largest first per dataset"
    echo "repeats=$REPEATS"
    echo "arm_order=rotated one position per cell and per repeat"
    echo "vertex_timeout_seconds=$VERTEX_TIMEOUT_SECONDS"
    echo "somei_timeout_seconds=$SOMEI_TIMEOUT_SECONDS (the measured batches, --measured-time-limit)"
    echo "somei_setup_seconds=$SOMEI_SETUP_SECONDS (more, for the reading and the untimed build)"
    echo "somei_extend_seconds=$SOMEI_EXTEND_SECONDS (a run set to finish by then runs until then, --measured-time-extend)"
    echo "run_timeout_seconds=$RUN_TIMEOUT_SECONDS (every other arm; 0 = none)"
    echo "skip_smaller_after_timeout=$SKIP_SMALLER_AFTER_TIMEOUT (vertex and somei arms)"
    echo "timeout_bound=total of a star/edge/vertex run >= its limit - the largest wall-minus-total of a finished star/edge/vertex run in the cell; somei timeouts unbounded"
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
    echo "warmup_batch_size=$WARMUP_BATCH_SIZE"
    echo "deadline_seconds=$DEADLINE_SECONDS"
    echo "min_run_seconds=$MIN_RUN_SECONDS"
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

# Seconds left before DEADLINE_SECONDS.
seconds_left() {
    local now
    now="$(date +%s)"
    echo $(( DEADLINE_SECONDS - (now - JOB_START) ))
}

failures=0
timeouts=0
out_of_memory_runs=0
run_count=0
skipped=0
excluded=0
# count name -> what the first finished run of the cell that reports it
# reported; every later run of the cell must match.
declare -A cell_counts=()
# dataset tag|arm -> the largest batch size at which that arm timed out.
declare -A timed_out_at=()

is_number() {
    [[ "$1" =~ ^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]]
}

# Whether a sliding-window run kept half of the edges it used: its first
# window is its warm-up, the batches after it replayed the rest, and the run
# was cut at the dataset's budget (none when 0).
window_is_half() {
    local stats_file="$1" target="$2" window warm processed cap
    window="$(extract_stat mixed window-edges "$stats_file")"
    warm="$(extract_stat graph warmup-edges "$stats_file")"
    processed="$(extract_stat graph processed-edges "$stats_file")"
    cap="$(extract_stat graph max-edges "$stats_file")"
    [[ "$window" =~ ^[0-9]+$ && "$processed" =~ ^[0-9]+$ && "$warm" == "$window" ]] || return 1
    if (( target > 0 )); then
        [[ "$cap" == "$target" ]] || return 1
    else
        [[ "$cap" == all ]] || return 1
    fi
    (( window == (window + processed) / 2 ))
}

# Whether a fraction a run printed is the one given (to the six digits it
# prints).
same_fraction() {
    awk -v a="$1" -v b="$2" 'BEGIN {
        if (a == "") { exit 1 }
        d = a - b
        if (d < 0) { d = -d }
        exit !(d <= 1e-6)
    }'
}

# Holds one count of a finished run to the first finished run of its cell that
# reported the same count, through run_cell's match and mismatch.
match_count() {
    local name="$1" value="$2"
    if [[ -z "$value" ]]; then
        match=no
        mismatch+="${mismatch:+, }$name missing"
    elif [[ -z "${cell_counts[$name]:-}" ]]; then
        cell_counts[$name]="$value"
    elif [[ "${cell_counts[$name]}" != "$value" ]]; then
        match=no
        mismatch+="${mismatch:+, }$name $value vs ${cell_counts[$name]}"
    elif [[ "$match" != no ]]; then
        match=yes
    fi
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

# The roots of one run's search trees, from the line of its own search; empty
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
# Returns 1 when the deadline leaves no room for it.
run_cell() {
    local dataset="$1" batch_size="$2" max_batches="$3" target="$4"
    local arm="$5" repeat="$6" position="$7" warmup="$8"
    local tag class stats_file start end wall exit_code status index
    local extra_arguments=()

    tag="$(dataset_tag_of "$dataset")"
    class="$(dataset_class_for "$(basename -- "$dataset")")"
    arm_arguments "$arm"
    local fraction=- seed=-
    if [[ "$WORKLOAD" == random-labels ]]; then
        fraction="$INSERT_FRACTION"
        seed="$LABEL_SEED"
    fi
    local -a identity=(
        "$tag" "$class" "$batch_size" "$warmup" "${max_batches:-0}" "$target"
        "$arm" "$ARM_KIND" "$ARM_FRAMEWORK" "$ARM_SCOPE" "$ARM_ALGORITHM" "$ARM_GATE"
        "$ARM_STRATEGY" "$ARM_PRUNING" "$ARM_CHANGES" "$WORKLOAD" "$fraction" "$seed"
        "$TIE_ORDER" "$position" "$repeat"
    )

    # Not run: a SOMEi arm outside SOMEI_DATASETS, or an arm that already
    # timed out here at this or a larger batch size and slows down as
    # batches shrink.
    local not_run=""
    if [[ "$ARM_KIND" == somei ]] && ! somei_runs_on "$(basename -- "$dataset")"; then
        not_run=excluded
    elif (( ARM_SKIPS == 1 && SKIP_SMALLER_AFTER_TIMEOUT == 1 )) && \
         [[ -n "${timed_out_at[$tag|$arm]:-}" ]] && \
         (( batch_size <= ${timed_out_at[$tag|$arm]} )); then
        not_run=skipped_after_timeout
    fi
    if [[ -n "$not_run" ]]; then
        if [[ "$not_run" == excluded ]]; then
            echo "[$(date -Is)] exclude dataset=$tag bs=$batch_size arm=$arm (not in SOMEI_DATASETS)"
            ((excluded += 1))
        else
            echo "[$(date -Is)] skip dataset=$tag bs=$batch_size arm=$arm (timed out at bs=${timed_out_at[$tag|$arm]})"
            ((skipped += 1))
        fi
        local -a blank=()
        local count=$(( ${#SUMMARY_COLUMNS[@]} - ${#identity[@]} - 1 ))
        for ((index = 0; index < count; ++index)); do blank+=('-'); done
        blank[0]="$not_run"
        blank[3]="$ARM_LIMIT"
        write_summary_row "${identity[@]}" "${blank[@]}" '-'
        return 0
    fi

    # The arm's own limit, cut to the time left before the deadline.
    local left limit cut_by_deadline=no
    left="$(seconds_left)"
    if (( left < MIN_RUN_SECONDS )); then
        return 1
    fi
    limit="$ARM_LIMIT"
    if (( limit == 0 || left < limit )); then
        limit="$left"
        cut_by_deadline=yes
    fi
    local -a limit_command=(timeout --signal=TERM --kill-after=60 "$limit")

    stats_file="${OUTPUT_DIR}/${tag}.b${batch_size}.${arm}.r${repeat}.log"
    if (( max_batches > 0 )); then
        extra_arguments+=("--max-batches=${max_batches}")
    fi
    if [[ "$WORKLOAD" == sliding-window ]] && (( target > 0 )); then
        extra_arguments+=("--max-edges=${target}")
    fi
    if (( warmup > 0 )); then
        extra_arguments+=("--warmup-edges=${warmup}" "--warmup-batch-size=${WARMUP_BATCH_SIZE}")
    fi

    echo "[$(date -Is)] run dataset=$tag bs=$batch_size warmup=$warmup arm=$arm repeat=$repeat position=$position limit=${limit}s"
    start="$(date +%s.%N)"
    "${LAUNCHER[@]}" "${limit_command[@]}" "$EXECUTABLE" "$dataset" \
        "${ARM_ARGUMENTS[@]}" \
        "--batch-size=${batch_size}" \
        "${extra_arguments[@]}" \
        --no-output \
        --stats \
        > "$stats_file" 2>&1
    exit_code=$?
    end="$(date +%s.%N)"
    wall="$(awk -v s="$start" -v e="$end" 'BEGIN {printf "%.3f", e - s}')"
    status="$(run_status "$exit_code" "$wall" "$stats_file" "$limit")"

    local label="" total_ms new_ms subsumed_ms index_ms netting_ms
    total_ms="$(extract_stat timing search-and-maintenance-ms "$stats_file")"
    if [[ "$ARM_KIND" == search ]]; then
        label="$ARM_FRAMEWORK"
        new_ms="$(extract_stat "$label" time-ms "$stats_file")"
        subsumed_ms="$(extract_stat subsumed time-ms "$stats_file")"
        index_ms="$(extract_stat timing clique-index-update-ms "$stats_file")"
        netting_ms="$(extract_stat mixed transient-netting-ms "$stats_file")"
    fi

    # Every arm is exact and draws the same labels, so no count of a cell may
    # depend on the arm.  A SOMEi arm without change listing has the edge and
    # remaining counts only.
    local inserted deleted created removed current key match='-' mismatch=""
    inserted="$(extract_stat mixed inserted-edges "$stats_file")"
    deleted="$(extract_stat mixed deleted-edges "$stats_file")"
    created="$(extract_stat graph new-maximal-cliques "$stats_file")"
    removed="$(extract_stat graph subsumed-maximal-cliques "$stats_file")"
    current="$(extract_stat graph current-maximal-cliques "$stats_file")"
    if [[ "$status" == ok ]]; then
        match=reference
        match_count inserted-edges "$inserted"
        match_count deleted-edges "$deleted"
        match_count remaining "$current"
        if [[ "$ARM_CHANGES" != off ]]; then
            match_count new "$created"
            match_count removed "$removed"
        fi
        if [[ "$ARM_KIND" == search ]]; then
            for key in destroyed recovered new subsumed; do
                match_count "mixed-$key" "$(extract_stat mixed "$key" "$stats_file")"
            done
            if [[ "$ARM_STRATEGY" == two-phase ]]; then
                match_count transient "$(extract_stat mixed transient "$stats_file")"
            fi
        fi
    fi

    # The run must have done what its arm says.
    local switches='-'
    if [[ "$status" == ok ]]; then
        switches=ok
        local other
        if [[ "$(extract_stat graph framework "$stats_file")" != "$ARM_FRAMEWORK" ]]; then
            switches=wrong-framework
        elif [[ "$(extract_stat graph update-mode "$stats_file")" != "$WORKLOAD" ]]; then
            switches=wrong-workload
        elif [[ "$(extract_stat graph maintenance "$stats_file")" != full || \
                "$(extract_stat graph maintenance-pipeline "$stats_file")" != batched ]]; then
            switches=not-batched-full-maintenance
        elif [[ "$(extract_stat graph tie-order "$stats_file")" != "$TIE_ORDER" ]]; then
            switches=wrong-tie-order
        elif [[ "$(extract_stat mixed strategy "$stats_file")" != "$ARM_STRATEGY" ]]; then
            switches=wrong-strategy
        elif [[ "$WORKLOAD" == random-labels ]] && \
             { ! same_fraction "$(extract_stat mixed insert-fraction "$stats_file")" "$INSERT_FRACTION" || \
               [[ "$(extract_stat mixed label-seed "$stats_file")" != "$LABEL_SEED" ]]; }; then
            switches=wrong-labels
        elif [[ "$WORKLOAD" == sliding-window ]] && ! window_is_half "$stats_file" "$target"; then
            switches=wrong-window
        elif [[ "$ARM_KIND" == somei ]]; then
            if [[ "$(extract_stat somei pruning "$stats_file")" != "$ARM_PRUNING" ]]; then
                switches=wrong-somei-pruning
            elif [[ "$(extract_stat somei changes "$stats_file")" != "$ARM_CHANGES" ]]; then
                switches=wrong-somei-changes
            elif [[ "$(extract_stat somei build "$stats_file")" != one-batch ]]; then
                switches=wrong-somei-build
            elif [[ "$(extract_stat somei build-update "$stats_file")" != bk ]]; then
                switches=wrong-somei-build-update
            elif [[ "$(extract_stat graph measured-time-limit "$stats_file")" != "$SOMEI_TIMEOUT_SECONDS" ]]; then
                switches=wrong-measured-time-limit
            elif [[ "$(extract_stat graph measured-time-extend "$stats_file")" != "$SOMEI_EXTEND_SECONDS" ]]; then
                switches=wrong-measured-time-extend
            elif [[ "$(extract_stat graph subsumed-algorithm "$stats_file")" != none || \
                    -n "$(extract_stat subsumed time-ms "$stats_file")" ]]; then
                switches=somei-ran-a-subsumed-search
            fi
        elif [[ "$(extract_stat graph subsumed-algorithm "$stats_file")" != "$ARM_ALGORITHM" || \
                "$(extract_stat subsumed algorithm "$stats_file")" != "$ARM_ALGORITHM" ]]; then
            switches=wrong-subsumed-search
        elif [[ "$ARM_FRAMEWORK" == star && \
                "$(extract_stat star branch-scope "$stats_file")" != "$ARM_SCOPE" ]]; then
            switches=wrong-branch-scope
        elif [[ "$(extract_stat index size-histogram "$stats_file")" != "$ARM_GATE" ]]; then
            switches=wrong-size-histogram
        elif [[ -z "$(extract_stat recovery checks "$stats_file")" ]]; then
            switches=no-recovery-statistics
        elif [[ "$ARM_STRATEGY" == single-pass && \
                "$(extract_stat mixed transient "$stats_file")" != 0 ]]; then
            switches=single-pass-netted-cliques
        else
            for other in edge vertex star; do
                [[ "$other" == "$ARM_FRAMEWORK" ]] && continue
                if [[ -n "$(extract_stat "$other" time-ms "$stats_file")" ]]; then
                    switches=two-searches-ran
                fi
            done
        fi
    fi

    local branches="" roots="" below='-'
    if [[ "$ARM_KIND" == search ]]; then
        branches="$(extract_stat "$label" recursive-calls "$stats_file")"
        roots="$(search_roots_of "$ARM_FRAMEWORK" "$stats_file")"
        if [[ "$branches" =~ ^[0-9]+$ && "$roots" =~ ^[0-9]+$ ]]; then
            below=$(( branches - roots ))
        fi
    fi
    # The fields of one line, empty where the arm has no such line (SOMEi's
    # mixed line holds its edge counts only).
    search_field() { [[ "$ARM_KIND" == search ]] && extract_stat "$label" "$1" "$stats_file"; }
    mixed_field() { [[ "$ARM_KIND" == search ]] && extract_stat mixed "$1" "$stats_file"; }
    edge_field() { [[ "$ARM_FRAMEWORK" == edge ]] && extract_stat edge "$1" "$stats_file"; }
    vertex_field() { [[ "$ARM_FRAMEWORK" == vertex ]] && extract_stat vertex "$1" "$stats_file"; }
    star_field() { [[ "$ARM_FRAMEWORK" == star ]] && extract_stat star "$1" "$stats_file"; }
    recovery_field() { [[ "$ARM_KIND" == search ]] && extract_stat recovery "$1" "$stats_file"; }
    subsumed_field() { [[ "$ARM_KIND" == search ]] && extract_stat subsumed "$1" "$stats_file"; }
    somei_field() { [[ "$ARM_KIND" == somei ]] && extract_stat somei "$1" "$stats_file"; }

    local -a measured=(
        "$status" "$exit_code" "$wall" "$limit"
        "${total_ms:-}" "${new_ms:-}" "${subsumed_ms:-}" "${index_ms:-}" "${netting_ms:-}"
        "$(mixed_field destroyed-search-ms)" "$(mixed_field new-search-ms)"
        "$(mixed_field subsumed-ms)" "$(mixed_field recovery-ms)"
        "$(extract_stat timing graph-update-ms "$stats_file")"
        "$(extract_stat graph warmup-ms "$stats_file")"
        "$(extract_stat timing input-preprocess-ms "$stats_file")"
        "$(somei_field scope-ms)" "$(somei_field rebuild-ms)" "$(somei_field change-ms)"
        "$match" "$switches"
        "$(extract_stat graph batches "$stats_file")"
        "$([[ "$WORKLOAD" == sliding-window ]] && extract_stat mixed window-edges "$stats_file")"
        "${inserted:-}" "${deleted:-}"
        "${created:-}" "${removed:-}" "${current:-}"
        "$(mixed_field destroyed)" "$(mixed_field recovered)" "$(mixed_field new)"
        "$(mixed_field subsumed)" "$(mixed_field transient)"
        "${branches:-}"
        "$(mixed_field destroyed-side-recursive-calls)" "$(mixed_field new-side-recursive-calls)"
        "${roots:-}" "$below"
        "$(search_field leaves)" "$(search_field emitted)"
        "$(search_field pivot-calls)" "$(search_field pivot-candidates)"
        "$(search_field local-builds)" "$(search_field max-local-vertices)"
        "$(search_field peak-local-bytes)"
        "$(edge_field seeds)"
        "$(vertex_field roots)" "$(vertex_field pruned-roots)" "$(vertex_field new-edge-prunes)"
        "$(star_field stars)" "$(star_field singleton-stars)"
        "$(star_field per-partner-star-edges)" "$(star_field branch-root-pivot-stars)"
        "$(star_field partner-prunes)"
        "$(recovery_field candidates)" "$(recovery_field below-two)"
        "$(recovery_field duplicates)" "$(recovery_field checks)"
        "$(recovery_field rejected)" "$(recovery_field accepted)"
        "$(recovery_field neighbor-reads)" "$(recovery_field adjacency-tests)"
        "$(recovery_field groups)" "$(recovery_field grouped-checks)"
        "$(recovery_field row-reads)"
        "$(subsumed_field size-gate-skipped-cliques)" "$(subsumed_field size-gate-skipped-probes)"
        "$(subsumed_field fallback-cliques)"
        "$([[ "$ARM_KIND" == search ]] && extract_stat wideresidual kernel-dynamic "$stats_file")"
        "$(somei_field rebuilt-headers)" "$(somei_field rolled-back-bodies)"
        "$(somei_field redo-steps)" "$(somei_field body-comparisons)"
        "$stats_file"
    )
    write_summary_row "${identity[@]}" "${measured[@]}"

    case "$status" in
        timeout)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tag" "$batch_size" "$arm" "$repeat" \
                "$limit" "$cut_by_deadline" "$stats_file" >> "$TIMEOUT_FILE"
            if [[ "$cut_by_deadline" == yes ]]; then
                note "Note: arm=$arm hit the deadline's ${limit}s limit on $tag bs=$batch_size; recorded in timeouts.tsv"
            else
                note "Note: arm=$arm hit its ${limit}s limit on $tag bs=$batch_size; recorded in timeouts.tsv"
            fi
            ((timeouts += 1))
            if (( ARM_SKIPS == 1 )) && \
               { [[ -z "${timed_out_at[$tag|$arm]:-}" ]] || \
                 (( batch_size > ${timed_out_at[$tag|$arm]} )); }; then
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
                note "Error: arm=$arm changed the counts on $tag bs=$batch_size: $mismatch (see $stats_file)"
                ((failures += 1))
            elif [[ "$switches" != ok ]]; then
                note "Error: arm=$arm did not do what its name says ($switches) on $tag bs=$batch_size; see $stats_file"
                ((failures += 1))
            elif ! is_number "${total_ms:-}" || \
                 { [[ "$ARM_KIND" == search ]] && \
                   { ! is_number "${new_ms:-}" || ! is_number "${subsumed_ms:-}" || \
                     ! is_number "${index_ms:-}" || ! is_number "${netting_ms:-}"; }; }; then
                note "Warning: a time is missing from $stats_file"
                ((failures += 1))
            elif [[ "$ARM_KIND" == search ]] && \
                 [[ ! "${branches:-}" =~ ^[0-9]+$ || "$below" == - ]]; then
                note "Warning: the search statistics are missing from $stats_file"
                ((failures += 1))
            fi
            ;;
    esac
    ((run_count += 1))
    return 0
}

# Repeats are the outer loop so that stopping early costs precision rather than
# whole datasets; inside a dataset the batch sizes run largest first; the arms
# of one cell run back to back, so the pairing inside a cell is as tight as it
# can be, their order rotated one position per cell and per repeat.
arm_count=${#arms[@]}
deadline_reached=0
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
            # A sliding window is cut at the budget by --max-edges (run_cell);
            # a random-labels window into whole batches.
            if [[ "$WORKLOAD" == sliding-window ]] && \
               (( target > 0 && target / 2 < batch_size )); then
                if (( repeat == 1 )); then
                    echo "Skipping dataset=$(dataset_tag_of "$dataset") bs=$batch_size: the $(( target / 2 ))-edge window is smaller than a batch" >&2
                fi
                ((skipped += 1))
                continue
            fi
            max_batches=0
            if (( target > 0 )) && [[ "$WORKLOAD" == random-labels ]]; then
                window=$(( target - warmup ))
                if (( window % batch_size != 0 )); then
                    if (( repeat == 1 )); then
                        echo "Skipping dataset=$(dataset_tag_of "$dataset") bs=$batch_size: the $window-edge window is not divisible by it" >&2
                    fi
                    ((skipped += 1))
                    continue
                fi
                max_batches=$(( window / batch_size ))
            fi
            shift_by=$(( (repeat - 1 + cell_number) % arm_count ))
            ((cell_number += 1))
            order=()
            for ((index = 0; index < arm_count; ++index)); do
                order+=("${arms[(index + shift_by) % arm_count]}")
            done
            cell_counts=()
            position=0
            for arm in "${order[@]}"; do
                if ! run_cell "$dataset" "$batch_size" "$max_batches" \
                        "$target" "$arm" "$repeat" "$position" "$warmup"; then
                    echo "Deadline reached; stopping after $run_count runs." >&2
                    deadline_reached=1
                    break 4
                fi
                ((position += 1))
            done
        done
    done
done

# Ratios are formed inside one cell and one repeat: time(from)/time(to) when
# both finished, a lower bound on the total when `from` hit its limit and `to`
# finished.  A cell's speedup is the median of its finished repeats, or of its
# bounds when none finished.  Runs that failed, disagreed on the counts or did
# not run what their arm says enter nothing.
awk -F'\t' -v OFS='\t' \
    -v arm_list="$ARMS" \
    -v arm_kind_list="$arm_kind_list" \
    -v arm_strategy_list="$arm_strategy_list" \
    -v experiment_list="$EXPERIMENTS" \
    -v repeats="$REPEATS" \
    -v method_arm="$METHOD_ARM" \
    -v results_file="$RESULTS_FILE" \
    -v times_file="$TIMES_FILE" \
    -v mixed_file="$MIXED_FILE" \
    -v tree_file="$TREE_FILE" \
    -v recovery_file="$RECOVERY_FILE" \
    -v table_time_file="$TABLE_TIME_FILE" \
    -v table_speedup_file="$TABLE_SPEEDUP_FILE" \
    -v table_tree_file="$TABLE_TREE_FILE" \
    -v geomean_file="$GEOMEAN_FILE" \
    -v noise_file="$NOISE_FILE" '
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
function dash(v) { return (v == "") ? "-" : v }
# The median of one stored value of one arm in one cell over its finished
# repeats, "" if none.
function arm_median(d, b, a, p,   here, n, i, r, list) {
    list = ""
    n = split(repeat_list[d, b], here, " ")
    for (i = 1; i <= n; ++i) {
        r = here[i]
        if (state[d, b, r, a] == "ok" && ((d, b, r, a, p) in value)) {
            list = list " " value[d, b, r, a, p]
        }
    }
    return median(list)
}
# A lower bound on the total of a search run that hit its limit: the limit
# less the most a finished search run of the cell spent outside its timers
# (reading, the build, the graph updates).  "" for SOMEi or without room.
function bound_ms(d, b, r, a,   room) {
    if (kind[a] != "search" || !((d, b) in outside_max)) { return "" }
    room = limit[d, b, r, a] * 1000 - outside_max[d, b]
    return (room > 0) ? room : ""
}
# What one arm did in one cell when it finished nowhere there: for the total,
# a lower bound from a timeout if there is one; otherwise what stopped it.
function arm_cell_text(d, b, a, p,   here, n, i, r, s, low, bound, any) {
    n = split(repeat_list[d, b], here, " ")
    low = ""; any = ""
    for (i = 1; i <= n; ++i) {
        r = here[i]
        s = state[d, b, r, a]
        if (s == "" || s == "ok") { continue }
        if (s == "timeout") {
            any = "timeout"
            if (p == "total") {
                bound = bound_ms(d, b, r, a)
                if (bound != "" && (low == "" || bound < low)) { low = bound }
            }
        } else if (any == "") {
            any = (s == "skipped") ? "skipped" : (s == "excluded") ? "excluded" : \
                (s == "oom") ? "out_of_memory" : "failed"
        }
    }
    if (low != "") { return ">=" fmt_ms(low) }
    return (any == "") ? "-" : any
}
function time_text(d, b, a, p,   m) {
    m = arm_median(d, b, a, p)
    return (m != "") ? fmt_ms(m) : arm_cell_text(d, b, a, p)
}
function branch_count(d, b, a) {
    return ((d, b, a) in branches) ? branches[d, b, a] : ""
}
function branch_ratio(d, b, from, to,   f, t) {
    f = branch_count(d, b, from); t = branch_count(d, b, to)
    if (f == "" || t == "" || t + 0 <= 0) { return "" }
    return f / t
}
# One comparison on one part in one cell, into X_status, X_fin, X_bound,
# X_from, X_to, X_speedup (a median of ratios, or of bounds when none
# finished), X_faster and X_bound_from.
function compare(d, b, from, to, p,   here, n, i, r, sf, st, ratios, bounds, from_ms, to_ms, bound_to, n_up, n_skip, n_excl, n_fail, n_open, n_oom_from, n_oom_to, room) {
    X_fin = 0; X_bound = 0; X_faster = 0
    X_from = ""; X_to = ""; X_speedup = ""; X_bound_from = ""
    ratios = ""; bounds = ""; from_ms = ""; to_ms = ""; bound_to = ""
    n_up = 0; n_skip = 0; n_excl = 0; n_fail = 0; n_open = 0; n_oom_from = 0; n_oom_to = 0
    n = split(repeat_list[d, b], here, " ")
    for (i = 1; i <= n; ++i) {
        r = here[i]
        sf = state[d, b, r, from]; st = state[d, b, r, to]
        if (sf == "" || st == "") { continue }
        if (sf == "ok" && st == "ok") {
            if (!((d, b, r, from, p) in value) || !((d, b, r, to, p) in value)) {
                ++n_fail
                continue
            }
            ++X_fin
            from_ms = from_ms " " value[d, b, r, from, p]
            to_ms = to_ms " " value[d, b, r, to, p]
            if (value[d, b, r, to, p] > 0) {
                ratios = ratios " " (value[d, b, r, from, p] / value[d, b, r, to, p])
            }
            if (value[d, b, r, to, p] < value[d, b, r, from, p]) { ++X_faster }
        } else if (sf == "timeout" && st == "ok") {
            room = (p == "total") ? bound_ms(d, b, r, from) : ""
            if (room == "" || value[d, b, r, to, "total"] <= 0) { ++n_open; continue }
            ++X_bound
            bounds = bounds " " (room / value[d, b, r, to, "total"])
            bound_to = bound_to " " value[d, b, r, to, "total"]
            if (X_bound_from == "" || room < X_bound_from) { X_bound_from = room }
        } else if (sf == "ok" && st == "timeout") {
            ++n_up
        } else if (sf == "oom") {
            ++n_oom_from
        } else if (st == "oom") {
            ++n_oom_to
        } else if (sf == "excluded" || st == "excluded") {
            ++n_excl
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
            (n_excl > 0) ? "excluded" : \
            (n_skip > 0) ? "skipped_after_timeout" : \
            (n_fail > 0) ? "failed" : "missing"
    }
}
# The log-speedup of one cell into its batch size, its group and the whole
# grid; finished cells count in both means, timed-out ones in the bound only.
function add(e, g, b, logv, is_finished,   key) {
    note_batch(e, g, b)
    split(e SUBSEP g SUBSEP b " " e SUBSEP g SUBSEP "all" " " e SUBSEP "all" SUBSEP "all", keys, " ")
    for (key = 1; key <= 3; ++key) {
        if (is_finished) { fin_sum[keys[key]] += logv; ++fin_count[keys[key]] }
        else { bound_sum[keys[key]] += logv; ++bound_count[keys[key]] }
    }
}
function add_branches(e, g, b, logv,   key) {
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
# A cell that enters no mean (skipped, excluded, out of memory), counted the
# same way.
function count_cell(counter, e, g, b) {
    ++counter[e, g, b]; ++counter[e, g, "all"]; ++counter[e, "all", "all"]
    note_batch(e, g, b)
}
# Tracks the spread of one value over the runs of one cell and repeat.
function track(name, d, b, r, v) {
    if (!number(v)) { return }
    v += 0
    if (!((name, d, b, r) in runs_of)) {
        runs_of[name, d, b, r] = 0
        low_of[name, d, b, r] = v
        high_of[name, d, b, r] = v
    }
    ++runs_of[name, d, b, r]
    if (v < low_of[name, d, b, r]) { low_of[name, d, b, r] = v }
    if (v > high_of[name, d, b, r]) { high_of[name, d, b, r] = v }
}
function spread_fields(name, d, b, r) {
    if (!((name, d, b, r) in runs_of)) { return "0" OFS "-" OFS "-" OFS "-" }
    return runs_of[name, d, b, r] OFS fmt_ms(low_of[name, d, b, r]) OFS \
        fmt_ms(high_of[name, d, b, r]) OFS \
        (low_of[name, d, b, r] > 0 ? sprintf("%.3f", high_of[name, d, b, r] / low_of[name, d, b, r]) : "-")
}
BEGIN {
    arm_count = split(arm_list, arm_order, " ")
    kind_count = split(arm_kind_list, kind_entry, " ")
    for (k = 1; k <= kind_count; ++k) {
        split(kind_entry[k], kind_pair, ":")
        kind[kind_pair[1]] = kind_pair[2]
    }
    strategy_count = split(arm_strategy_list, strategy_entry, " ")
    for (k = 1; k <= strategy_count; ++k) {
        split(strategy_entry[k], strategy_pair, ":")
        strategy[strategy_pair[1]] = strategy_pair[2]
    }
    experiment_count = split(experiment_list, entry, " ")
    for (e = 1; e <= experiment_count; ++e) {
        split(entry[e], name_and_spec, "=")
        split(name_and_spec[2], spec, ":")
        experiment_name[e] = name_and_spec[1]
        experiment_from[e] = spec[1]
        experiment_to[e] = spec[2]
        experiment_part[e] = spec[3]
    }
    # What each finished run stores, and the summary column it comes from.
    stored_count = split("total new_clique subsumed index_update netting destroyed_search new_search subsumed_search recovery graph_update build somei_scope somei_rebuild somei_change wall", stored, " ")
    for (k = 1; k <= stored_count; ++k) { stored_column[stored[k]] = stored[k] "_ms" }
    stored_column["wall"] = "wall_seconds"
    tree_field_count = split("search_roots branches_below_roots leaves emitted pivot_calls pivot_candidates vertex_pruned_roots vertex_new_edge_prunes star_partner_prunes destroyed_side_branches new_side_branches", tree_field, " ")
    recovery_field_count = split("mixed_destroyed mixed_recovered recovery_candidates recovery_below_two recovery_duplicates recovery_checks recovery_rejected recovery_accepted recovery_neighbor_reads recovery_adjacency_tests recovery_groups recovery_grouped_checks recovery_row_reads", recovery_field, " ")
    mixed_field_count = split("inserted_edges deleted_edges new_cliques removed_cliques current_maximal_cliques mixed_destroyed mixed_recovered mixed_new mixed_subsumed mixed_transient", mixed_field, " ")
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
    searched = (kind[a] == "search")
    if (s == "ok" && $(column["clique_counts_match"]) != "no" && \
        $(column["switches_check"]) == "ok" && number($(column["total_ms"])) && \
        (!searched || (number($(column["new_clique_ms"])) && \
                       number($(column["subsumed_ms"])) && count($(column["branches"]))))) {
        state[d, b, r, a] = "ok"
        for (k = 1; k <= stored_count; ++k) {
            v = $(column[stored_column[stored[k]]])
            if (number(v)) { value[d, b, r, a, stored[k]] = v + 0 }
        }
        # The counts are the same in every run of a cell (checked), so the
        # first finished run gives them.
        if (!((d, b, a) in has_mix)) {
            has_mix[d, b, a] = 1
            for (k = 1; k <= mixed_field_count; ++k) {
                mix[d, b, a, mixed_field[k]] = $(column[mixed_field[k]])
            }
        }
        if (searched) {
            outside = value[d, b, r, a, "wall"] * 1000 - value[d, b, r, a, "total"]
            if (!((d, b) in outside_max) || outside > outside_max[d, b]) { outside_max[d, b] = outside }
            # Every search arm builds with the same insertion and changes the
            # graph by the same edges.  SOMEi does too, but its own structures
            # keep the graph out of the caches (deleting, its removals took
            # 1.4-1.8 times as long in the local smoke), so it is left out of
            # both.
            track("build", d, b, r, $(column["build_ms"]))
            track("graph_update", d, b, r, $(column["graph_update_ms"]))
            # The search and the recovery are deterministic: the first
            # finished run gives them and the later ones must agree.
            n_branches = $(column["branches"]) + 0
            if (!((d, b, a) in branches)) {
                branches[d, b, a] = n_branches
                branches_min[d, b, a] = n_branches
                branches_max[d, b, a] = n_branches
                for (k = 1; k <= tree_field_count; ++k) {
                    tree[d, b, a, tree_field[k]] = $(column[tree_field[k]])
                }
                for (k = 1; k <= recovery_field_count; ++k) {
                    recovery[d, b, a, recovery_field[k]] = $(column[recovery_field[k]])
                }
            }
            if (n_branches < branches_min[d, b, a]) { branches_min[d, b, a] = n_branches }
            if (n_branches > branches_max[d, b, a]) { branches_max[d, b, a] = n_branches }
        }
    } else if (s == "timeout") {
        state[d, b, r, a] = "timeout"
    } else if (s == "skipped_after_timeout") {
        state[d, b, r, a] = "skipped"
    } else if (s == "excluded") {
        state[d, b, r, a] = "excluded"
    } else if (s == "out_of_memory") {
        state[d, b, r, a] = "oom"
    } else {
        state[d, b, r, a] = "failed"
    }
}
END {
    # table_time.tsv: the total of every arm, and the new-clique and subsumed
    # times of every search arm.
    header = "dataset" OFS "batch_size"
    for (m = 1; m <= arm_count; ++m) {
        a = arm_order[m]
        header = header OFS a "_total_ms"
        if (kind[a] == "search") { header = header OFS a "_new_clique_ms" OFS a "_subsumed_ms" }
    }
    print header > table_time_file
    # table_speedup.tsv: every comparison.
    header = "dataset" OFS "batch_size"
    for (e = 1; e <= experiment_count; ++e) { header = header OFS experiment_name[e] }
    print header > table_speedup_file
    # table_tree.tsv: the branches of every search arm and the branch ratios
    # of the new_clique comparisons.
    header = "dataset" OFS "batch_size"
    for (m = 1; m <= arm_count; ++m) {
        a = arm_order[m]
        if (kind[a] == "search") { header = header OFS a "_branches" }
    }
    for (e = 1; e <= experiment_count; ++e) {
        if (experiment_part[e] == "new_clique") { header = header OFS experiment_name[e] "_branches" }
    }
    print header > table_tree_file
    print "dataset", "batch_size", "arm", "branches", "branches_vs_" method_arm, \
        "destroyed_side", "new_side", "search_roots", "branches_below_roots", "leaves", \
        "emitted", "pivot_calls", "pivot_candidates", "pruned_roots", "new_edge_prunes", \
        "partner_prunes", "runs_agree" > tree_file
    print "dataset", "batch_size", "arm", "strategy", "destroyed", "recovered", \
        "candidates", "below_two", "already_accepted", "checks", "rejected", \
        "accepted", "neighbor_reads", "adjacency_tests", "groups", "grouped_checks", \
        "row_reads", "recovery_ms" > recovery_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        line_time = d OFS b
        line_tree = d OFS b
        for (m = 1; m <= arm_count; ++m) {
            a = arm_order[m]
            line_time = line_time OFS time_text(d, b, a, "total")
            if (kind[a] == "search") {
                line_time = line_time OFS time_text(d, b, a, "new_clique") OFS time_text(d, b, a, "subsumed")
                text = arm_cell_text(d, b, a, "branches")
                line_tree = line_tree OFS (((d, b, a) in branches) ? branches[d, b, a] : text)
            }
        }
        print line_time > table_time_file
        line_speedup = d OFS b
        for (e = 1; e <= experiment_count; ++e) {
            compare(d, b, experiment_from[e], experiment_to[e], experiment_part[e])
            text = (X_status == "finished") ? fmt_ratio(X_speedup) : \
                (X_status == "from_timed_out") ? ">=" fmt_ratio(X_speedup) : \
                (X_status == "missing") ? "-" : X_status
            line_speedup = line_speedup OFS text
            if (experiment_part[e] == "new_clique") {
                line_tree = line_tree OFS fmt_ratio(branch_ratio(d, b, experiment_from[e], experiment_to[e]))
            }
        }
        print line_speedup > table_speedup_file
        print line_tree > table_tree_file
        for (k = 1; k <= arm_count; ++k) {
            a = arm_order[k]
            if (!((d, b, a) in branches)) { continue }
            print d, b, a, branches[d, b, a], \
                (kind[method_arm] == "search" ? fmt_ratio(branch_ratio(d, b, a, method_arm)) : "-"), \
                dash(tree[d, b, a, "destroyed_side_branches"]), dash(tree[d, b, a, "new_side_branches"]), \
                dash(tree[d, b, a, "search_roots"]), dash(tree[d, b, a, "branches_below_roots"]), \
                dash(tree[d, b, a, "leaves"]), dash(tree[d, b, a, "emitted"]), \
                dash(tree[d, b, a, "pivot_calls"]), dash(tree[d, b, a, "pivot_candidates"]), \
                dash(tree[d, b, a, "vertex_pruned_roots"]), dash(tree[d, b, a, "vertex_new_edge_prunes"]), \
                dash(tree[d, b, a, "star_partner_prunes"]), \
                (branches_min[d, b, a] == branches_max[d, b, a] ? "yes" : "no") > tree_file
            print d, b, a, strategy[a], \
                dash(recovery[d, b, a, "mixed_destroyed"]), dash(recovery[d, b, a, "mixed_recovered"]), \
                dash(recovery[d, b, a, "recovery_candidates"]), dash(recovery[d, b, a, "recovery_below_two"]), \
                dash(recovery[d, b, a, "recovery_duplicates"]), dash(recovery[d, b, a, "recovery_checks"]), \
                dash(recovery[d, b, a, "recovery_rejected"]), dash(recovery[d, b, a, "recovery_accepted"]), \
                dash(recovery[d, b, a, "recovery_neighbor_reads"]), dash(recovery[d, b, a, "recovery_adjacency_tests"]), \
                dash(recovery[d, b, a, "recovery_groups"]), dash(recovery[d, b, a, "recovery_grouped_checks"]), \
                dash(recovery[d, b, a, "recovery_row_reads"]), \
                fmt_ms(arm_median(d, b, a, "recovery")) > recovery_file
        }
    }

    # mixed.tsv: what the batches changed, and the parts by side.  The share
    # of transient cliques is over the net change, destroyed + recovered +
    # new + subsumed, and only two-phase runs net any.
    print "dataset", "batch_size", "arm", "strategy", "inserted_edges", "deleted_edges", \
        "all_new", "all_removed", "remaining", "destroyed", "recovered", "new", \
        "subsumed", "transient", "transient_share", "destroyed_search_ms", \
        "new_search_ms", "subsumed_search_ms", "recovery_ms", "index_update_ms", \
        "netting_ms", "total_ms" > mixed_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        for (m = 1; m <= arm_count; ++m) {
            a = arm_order[m]
            if (!((d, b, a) in has_mix)) { continue }
            share = "-"
            if (strategy[a] == "two-phase" && count(mix[d, b, a, "mixed_transient"])) {
                net = mix[d, b, a, "mixed_destroyed"] + mix[d, b, a, "mixed_recovered"] + \
                    mix[d, b, a, "mixed_new"] + mix[d, b, a, "mixed_subsumed"]
                if (net > 0) { share = sprintf("%.4f", mix[d, b, a, "mixed_transient"] / net) }
            }
            print d, b, a, strategy[a], \
                dash(mix[d, b, a, "inserted_edges"]), dash(mix[d, b, a, "deleted_edges"]), \
                dash(mix[d, b, a, "new_cliques"]), dash(mix[d, b, a, "removed_cliques"]), \
                dash(mix[d, b, a, "current_maximal_cliques"]), \
                dash(mix[d, b, a, "mixed_destroyed"]), dash(mix[d, b, a, "mixed_recovered"]), \
                dash(mix[d, b, a, "mixed_new"]), dash(mix[d, b, a, "mixed_subsumed"]), \
                dash(mix[d, b, a, "mixed_transient"]), share, \
                fmt_ms(arm_median(d, b, a, "destroyed_search")), \
                fmt_ms(arm_median(d, b, a, "new_search")), \
                fmt_ms(arm_median(d, b, a, "subsumed_search")), \
                fmt_ms(arm_median(d, b, a, "recovery")), \
                fmt_ms(arm_median(d, b, a, "index_update")), \
                fmt_ms(arm_median(d, b, a, "netting")), \
                fmt_ms(arm_median(d, b, a, "total")) > mixed_file
        }
    }

    # times.tsv: one row per cell and arm.
    printf "dataset\tbatch_size\tarm\tfinished\ttimed_out\tskipped\texcluded\tout_of_memory\tfailed" > times_file
    for (r = 1; r <= repeats; ++r) { printf "\tr%d_total_ms", r > times_file }
    printf "\ttotal_ms\tnew_clique_ms\tsubsumed_ms\tindex_update_ms\tnetting_ms\tgraph_update_ms\tbuild_ms" > times_file
    printf "\tsomei_scope_ms\tsomei_rebuild_ms\tsomei_change_ms\ttotal_min_ms\ttotal_max_ms\tspread\twall_seconds\n" > times_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        for (m = 1; m <= arm_count; ++m) {
            a = arm_order[m]
            n_ok = 0; n_timeout = 0; n_skip = 0; n_excl = 0; n_oom = 0; n_fail = 0
            per_repeat = ""; low = ""; high = ""
            for (r = 1; r <= repeats; ++r) {
                s = state[d, b, r, a]
                if (s == "ok") {
                    t = value[d, b, r, a, "total"]
                    ++n_ok
                    per_repeat = per_repeat "\t" fmt_ms(t)
                    if (low == "" || t < low) { low = t }
                    if (high == "" || t > high) { high = t }
                } else {
                    if (s == "timeout") { ++n_timeout }
                    else if (s == "skipped") { ++n_skip }
                    else if (s == "excluded") { ++n_excl }
                    else if (s == "oom") { ++n_oom }
                    else if (s == "failed") { ++n_fail }
                    per_repeat = per_repeat "\t" ((s == "") ? "-" : s)
                }
            }
            if (n_ok + n_timeout + n_skip + n_excl + n_oom + n_fail == 0) { continue }
            printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d%s", d, b, a, n_ok, n_timeout, \
                n_skip, n_excl, n_oom, n_fail, per_repeat > times_file
            printf "\t%s\t%s\t%s\t%s\t%s\t%s\t%s", time_text(d, b, a, "total"), \
                fmt_ms(arm_median(d, b, a, "new_clique")), fmt_ms(arm_median(d, b, a, "subsumed")), \
                fmt_ms(arm_median(d, b, a, "index_update")), fmt_ms(arm_median(d, b, a, "netting")), \
                fmt_ms(arm_median(d, b, a, "graph_update")), \
                fmt_ms(arm_median(d, b, a, "build")) > times_file
            printf "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", fmt_ms(arm_median(d, b, a, "somei_scope")), \
                fmt_ms(arm_median(d, b, a, "somei_rebuild")), fmt_ms(arm_median(d, b, a, "somei_change")), \
                fmt_ms(low), fmt_ms(high), \
                ((low != "" && low > 0) ? sprintf("%.3f", high / low) : "-"), \
                fmt_ms(arm_median(d, b, a, "wall")) > times_file
        }
    }

    # noise.tsv: the same work in several runs of one cell and repeat.
    print "dataset", "batch_size", "repeat", "build_runs", "build_min_ms", "build_max_ms", \
        "build_spread", "graph_update_runs", "graph_update_min_ms", "graph_update_max_ms", \
        "graph_update_spread" > noise_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        n = split(repeat_list[d, b], here, " ")
        for (i = 1; i <= n; ++i) {
            r = here[i]
            if (!(("build", d, b, r) in runs_of) && !(("graph_update", d, b, r) in runs_of)) { continue }
            print d, b, r, spread_fields("build", d, b, r), spread_fields("graph_update", d, b, r) > noise_file
        }
    }

    # results.tsv and the means.
    print "dataset", "batch_size", "experiment", "from_arm", "to_arm", "part", "status", \
        "finished_pairs", "bounded_pairs", "from_ms", "to_ms", "speedup", \
        "to_faster_in", "from_branches", "to_branches", "branch_ratio" > results_file
    for (c = 1; c <= cells; ++c) {
        split(cell_key[c], parts, SUBSEP)
        d = parts[1]; b = parts[2]
        for (e = 1; e <= experiment_count; ++e) {
            from = experiment_from[e]; to = experiment_to[e]; p = experiment_part[e]
            compare(d, b, from, to, p)
            if (X_status == "missing") { continue }
            from_branches = "-"; to_branches = "-"; ratio = ""
            if (p == "new_clique") {
                from_branches = dash(branch_count(d, b, from))
                to_branches = dash(branch_count(d, b, to))
                ratio = branch_ratio(d, b, from, to)
            }
            if (X_status == "finished") {
                print d, b, experiment_name[e], from, to, p, X_status, X_fin, X_bound, \
                    fmt_ms(X_from), fmt_ms(X_to), fmt_ratio(X_speedup), X_faster "/" X_fin, \
                    from_branches, to_branches, fmt_ratio(ratio) > results_file
                if (X_speedup != "" && X_speedup > 0) { add(e, class_of[d], b, log(X_speedup), 1) }
            } else if (X_status == "from_timed_out") {
                print d, b, experiment_name[e], from, to, p, X_status, 0, X_bound, \
                    ">=" fmt_ms(X_bound_from), fmt_ms(X_to), ">=" fmt_ratio(X_speedup), \
                    X_bound "/" X_bound, "-", to_branches, "-" > results_file
                add(e, class_of[d], b, log(X_speedup), 0)
            } else {
                print d, b, experiment_name[e], from, to, p, X_status, 0, 0, "-", "-", "-", "-", \
                    from_branches, to_branches, fmt_ratio(ratio) > results_file
                if (X_status == "skipped_after_timeout") { count_cell(skipped_cells, e, class_of[d], b) }
                if (X_status == "excluded") { count_cell(excluded_cells, e, class_of[d], b) }
                if (X_status ~ /out_of_memory/) { count_cell(oom_cells, e, class_of[d], b) }
            }
            if (ratio != "" && ratio > 0) { add_branches(e, class_of[d], b, log(ratio)) }
        }
    }

    print "experiment", "from_arm", "to_arm", "part", "datasets", "batch_size", \
        "finished_cells", "speedup_geomean_finished", "timed_out_cells", \
        "speedup_geomean_lower_bound", "skipped_cells", "excluded_cells", \
        "out_of_memory_cells", "branch_cells", "branch_ratio_geomean" > geomean_file
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
                    excluded_cells[e, g, b] + oom_cells[e, g, b] + br_count[e, g, b] == 0) { continue }
                geo_fin = (fin_count[e, g, b] > 0) ? \
                    sprintf("%.4f", exp(fin_sum[e, g, b] / fin_count[e, g, b])) : "-"
                both = fin_count[e, g, b] + bound_count[e, g, b]
                geo_low = (both > 0) ? \
                    sprintf("%.4f", exp((fin_sum[e, g, b] + bound_sum[e, g, b]) / both)) : "-"
                geo_br = (br_count[e, g, b] > 0) ? \
                    sprintf("%.4f", exp(br_sum[e, g, b] / br_count[e, g, b])) : "-"
                print experiment_name[e], experiment_from[e], experiment_to[e], \
                    experiment_part[e], g, b, fin_count[e, g, b] + 0, geo_fin, \
                    bound_count[e, g, b] + 0, geo_low, skipped_cells[e, g, b] + 0, \
                    excluded_cells[e, g, b] + 0, oom_cells[e, g, b] + 0, \
                    br_count[e, g, b] + 0, geo_br > geomean_file
            }
        }
    }
}
' "$SUMMARY_FILE"

# somei_progress.tsv: how far every SOMEi run got, read off its --progress
# lines.  measured_begin_seconds is when its measured batches began, after
# the reading and the untimed build ("-": the limit came first); batches_done
# of batches, with their total so far, is what a run its limit stopped had
# finished.  linear_estimate_ms scales that total to all the batches, a rough
# figure only: batches differ in cost.  measured_limit_seconds is the limit on
# its measured batches and measured_extend_seconds the one a run set to finish
# by then could have (2026-11-06); past_limit says whether its measured
# batches ran past measured_limit_seconds (by more than the 2 s an alarm
# may overshoot), as only that extension lets them;
# stopped_by is what ended a timeout -- measured-time-limit (the
# run's own limit; stopped_measured_seconds is the measured time it then had)
# or wall-limit (the wall limit or the deadline came first: the warning below).
SOMEI_PROGRESS_FILE="${OUTPUT_DIR}/somei_progress.tsv"
printf 'dataset\tbatch_size\trepeat\tstatus\tlimit_seconds\tmeasured_begin_seconds\tbatches_done\tbatches\ttotal_ms_done\tlinear_estimate_ms\tstats_file\tmeasured_limit_seconds\tstopped_by\tmeasured_extend_seconds\tpast_limit\tstopped_measured_seconds\n' \
    > "$SOMEI_PROGRESS_FILE"
awk -F'\t' -v OFS='\t' '
NR == 1 {
    for (i = 1; i <= NF; ++i) { column[$i] = i }
    next
}
$(column["kind"]) == "somei" && $(column["stats_file"]) != "-" {
    print $(column["dataset"]), $(column["batch_size"]), $(column["repeat"]), \
        $(column["status"]), $(column["time_limit_seconds"]), $(column["stats_file"])
}' "$SUMMARY_FILE" |
while IFS=$'\t' read -r p_dataset p_batch p_repeat p_status p_limit p_log; do
    [[ -f "$p_log" ]] || continue
    awk -F', ' -v OFS='\t' -v d="$p_dataset" -v b="$p_batch" -v r="$p_repeat" \
        -v s="$p_status" -v limit="$p_limit" -v logfile="$p_log" \
        -v measured="$SOMEI_TIMEOUT_SECONDS" -v extend="$SOMEI_EXTEND_SECONDS" '
    /^progress: measured-begin/ {
        for (i = 2; i <= NF; ++i) {
            if ($i ~ /^measured-batches=/) {
                split(substr($i, 18), counts, "/")
                total = counts[2]
            }
            if ($i ~ /^wall-ms=/) { begin = substr($i, 9) / 1000 }
        }
    }
    /^progress: measured-batches=/ {
        split(substr($1, 28), counts, "/")
        done = counts[1]
        total = counts[2]
        for (i = 2; i <= NF; ++i) {
            if ($i ~ /^search-and-maintenance-ms=/) { ms = substr($i, 27) }
            if ($i ~ /^wall-ms=/) { last_wall = substr($i, 9) / 1000 }
        }
    }
    /^progress: measured-time-limit/ {
        stopped = "measured-time-limit"
        for (i = 2; i <= NF; ++i) {
            if ($i ~ /^measured-wall-ms=/) { stopped_ms = substr($i, 18) }
        }
    }
    END {
        estimate = (done > 0 && total > 0 && ms != "") ? sprintf("%.0f", ms * total / done) : "-"
        if (stopped == "") { stopped = (s == "timeout" ? "wall-limit" : "-") }
        if (stopped_ms != "") { ran = stopped_ms / 1000 }
        else if (last_wall != "" && begin != "") { ran = last_wall - begin }
        else { ran = 0 }
        print d, b, r, s, limit, (begin == "" ? "-" : sprintf("%.1f", begin)), \
            (done == "" ? 0 : done), (total == "" ? "-" : total), \
            (ms == "" ? "-" : ms), estimate, logfile, measured, stopped, extend,
            (ran > measured + 2 ? "yes" : "no"),
            (stopped_ms == "" ? "-" : sprintf("%.1f", stopped_ms / 1000))
    }' "$p_log"
done >> "$SOMEI_PROGRESS_FILE"

{
    echo "finished=$(date -Is)"
    echo "run_count=$run_count"
    echo "failures=$failures"
    echo "timeouts=$timeouts"
    echo "out_of_memory_runs=$out_of_memory_runs"
    echo "skipped=$skipped"
    echo "excluded=$excluded"
    echo "deadline_reached=$deadline_reached"
    echo "elapsed_seconds=$(( $(date +%s) - JOB_START ))"
    echo "table_time_file=$TABLE_TIME_FILE"
    echo "table_speedup_file=$TABLE_SPEEDUP_FILE"
    echo "table_tree_file=$TABLE_TREE_FILE"
    echo "mixed_file=$MIXED_FILE"
    echo "results_file=$RESULTS_FILE"
    echo "geomean_file=$GEOMEAN_FILE"
    echo "times_file=$TIMES_FILE"
    echo "tree_file=$TREE_FILE"
    echo "recovery_file=$RECOVERY_FILE"
    echo "noise_file=$NOISE_FILE"
    echo "timeout_file=$TIMEOUT_FILE"
    echo "summary_file=$SUMMARY_FILE"
    echo "somei_progress_file=$SOMEI_PROGRESS_FILE"
} >> "${OUTPUT_DIR}/manifest.txt"

show_table() {
    if command -v column > /dev/null 2>&1; then
        column -t -s $'\t' "$1"
    else
        cat "$1"
    fi
}

echo "Finished: $(date -Is)"
echo "Runs: $run_count; failures: $failures; timeouts: $timeouts; out of memory: $out_of_memory_runs; skipped: $skipped; excluded: $excluded"
if (( deadline_reached == 1 )); then
    echo "The deadline stopped the sweep before every cell had run."
fi
echo "Time table: $TABLE_TIME_FILE"
echo "Speedups: $TABLE_SPEEDUP_FILE"
echo "Branch table: $TABLE_TREE_FILE"
echo "Changes and parts by side: $MIXED_FILE"
echo "Times: $TIMES_FILE"
echo "Search statistics: $TREE_FILE"
echo "Recovery statistics: $RECOVERY_FILE"
echo "Results: $RESULTS_FILE"
echo "Geometric means: $GEOMEAN_FILE"
echo "Noise: $NOISE_FILE"
echo "Timeouts: $TIMEOUT_FILE"
echo "Summary: $SUMMARY_FILE"
echo "SOMEi progress: $SOMEI_PROGRESS_FILE"
echo
echo "How far every SOMEi run got (from its progress lines; a run its limit stopped shows the batches it finished):"
show_table "$SOMEI_PROGRESS_FILE"
short_somei="$(awk -F'\t' 'NR > 1 && $4 == "timeout" && $13 != "measured-time-limit" {
    print "  " $1 " b" $2 " repeat " $3 ": measured from " $6 " s of a " $5 " s run"
}' "$SOMEI_PROGRESS_FILE")"
if [[ -n "$short_somei" ]]; then
    echo "Warning: SOMEi runs that the wall limit or the deadline stopped before ${SOMEI_TIMEOUT_SECONDS}s of measured batches (a larger SOMEI_SETUP_SECONDS gives them their time):"
    echo "$short_somei"
fi
echo
echo "Median mixed-batch times (ms) of every arm -- total, then new clique and subsumed for the searches; >= a lower bound from a timeout:"
show_table "$TABLE_TIME_FILE"
echo
echo "Speedups time(from)/time(to) (above 1: to is faster):"
show_table "$TABLE_SPEEDUP_FILE"
echo
echo "Branches of every search, and branches(from)/branches(to) of the new_clique comparisons (above 1: to opens fewer):"
show_table "$TABLE_TREE_FILE"
echo
echo "The changes and the median parts by side (transient_share: the netted cliques over the net change):"
show_table "$MIXED_FILE"
echo
echo "Geometric means of the per-cell speedups and branch ratios:"
show_table "$GEOMEAN_FILE"
echo
