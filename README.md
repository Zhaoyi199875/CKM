# CKM: Maximal Clique Maintenance in Dynamic Graphs

This repository contains the code for the paper
*Beyond One Edge at a Time: Maximal Clique Maintenance in Dynamic Graphs*.

CKM (Cover and Kernel Maintenance) maintains the maximal cliques of a dynamic
graph under batches of edge insertions, deletions or both. After each
insertion batch, CoverSeed finds the new maximal cliques and KernelEnum finds
the subsumed cliques. Deletion and hybrid batches are reduced to insertions.
One binary, `immce_star`, implements CKM and every method compared in the
paper.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | C++20 sources: the maintenance engine, CKM (CoverSeed and KernelEnum), IMCE (EdgeSeed and Edge-Based Splitting), VertexSeed and SOMEi |
| `tests/` | Tests run by `make check` |
| `preprocessing/` | Converter from KONECT edge lists to the input format |
| `scripts/` | SLURM scripts used for the experiments in Section VI |

## Build

Requirements: Linux, GCC 13 or newer and GNU make.

```bash
make            # builds build/immce_star
make check      # optional: runs the tests
```

## Datasets

The experiments use six timestamped graphs from KONECT (http://konect.cc/).

| Name in the paper | KONECT network | Edges used |
|---|---|---|
| cit-dblp | `dblp-cite` | all |
| youtube | `youtube-u-growth` | all |
| dblp-coauthor | `dblp_coauthor` | all |
| flickr-growth | `flickr-growth` | first 1M |
| wikipedia-growth | `wikipedia-growth` | first 12M |
| stack-overflow | `sx-stackoverflow` | first 1M |

The scripts keep these prefixes by passing `--max-batches`, the prefix length
divided by the batch size.

Convert a downloaded `out.<name>` file with

```bash
cd preprocessing
make
./preprocess_dynamic_graph out.flickr-growth    # writes processed_flickr-growth.txt
```

The converter keeps the two endpoint columns and the timestamp column, drops
records without a timestamp, and maps the timestamps to consecutive integers
in increasing order. `immce_star` then treats the graph as undirected, drops
self-loops, keeps the earliest timestamp of each duplicated edge, and
processes the edges in timestamp order (ties in file order).

## Running

```bash
./build/immce_star processed_dblp-cite.txt --batch-size=1000 \
    --framework=star --subsumed-algorithm=residual-kernel \
    --maintenance=full --update-mode=insert --stats --no-output
```

`--stats` prints the counts and the timing breakdown. `./build/immce_star --help`
lists all options.

Methods in the paper and the corresponding options:

| Paper | Options |
|---|---|
| CKM | `--framework=star --subsumed-algorithm=residual-kernel --maintenance=full` |
| IMCE | `--framework=edge --subsumed-algorithm=optimized --maintenance=full` |
| SOMEi | `--framework=somei --maintenance=full` |
| EdgeSeed, VertexSeed, CoverSeed (new maximal cliques only) | `--framework=edge`, `--framework=vertex`, `--framework=star`, each with `--maintenance=new-only` |
| Edge-Based Splitting, KernelEnum | `--subsumed-algorithm=optimized`, `--subsumed-algorithm=residual-kernel` |
| Pivot, Anchor, Anchor-or-Pivot branching | `--star-branching=pivot`, `--star-branching=anchor`, `--star-branching=rule` (default) |
| +Cover (Section VI-E) | `--framework=vertex --vertex-seeds=cover` |
| +Assignment | `--framework=star --star-root-filter=off --star-branching=pivot` |
| +Domain | `--framework=star --star-root-filter=on --star-branching=pivot` |

Workloads:

| Paper | Options |
|---|---|
| Insertion batches | `--update-mode=insert`: the batches are inserted into the empty graph in timestamp order |
| Deletion batches | `--update-mode=reverse-delete`: starting from the graph with all edges, the batches are deleted in reverse order |
| Hybrid batches | `--update-mode=sliding-window`: the window first holds the first half of the edges, and each batch inserts the next `b` edges and deletes the `b` oldest edges |

`--batch-size=b` sets the batch size and `--max-edges=N` keeps the first `N`
edges. `--measured-time-limit=SECONDS` stops a run after the given measured
maintenance time (the experiments use 18000 seconds).

## Experiment scripts

The scripts in `scripts/` were run on a SLURM cluster. Their headers document
the compared arms, datasets, batch sizes and output files. Set the `#SBATCH`
output paths and `PROJECT_DIR` for your environment before submitting them.

| Section | Script |
|---|---|
| VI-B (insertions, deletions, hybrid batches) | `submit_overall_insert_times.sh`, `submit_overall_delete_times.sh`, `submit_overall_hybrid_times.sh`, `submit_overall_hybrid_coauthor_long_times.sh` |
| VI-C (EdgeSeed, VertexSeed and CoverSeed) | `submit_edge_vertex_crossover_times.sh` |
| VI-D (Edge-Based Splitting and KernelEnum) | `submit_overall_insert_times.sh`, with the subsumed-search arms selected through `ARMS` |
| VI-E and VI-F (ablation and branching) | `submit_star_ablation_times.sh` |
