# CKM: Maximal Clique Maintenance in Dynamic Graphs

This repository contains the code for the paper
*Beyond One Edge at a Time: Maximal Clique Maintenance in Dynamic Graphs*.

CKM (Cover and Kernel Maintenance) maintains the maximal cliques of a dynamic
graph under batches of edge insertions, deletions or both. After each
insertion batch, CoverSeed (Section III) finds the new maximal cliques and
KernelEnum (Section IV) finds the subsumed cliques. Deletion and hybrid
batches are handled as described in Section II-D. The binary `ckm`
implements CKM, and VertexSeed (Section III-A) can take the place of
CoverSeed.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | C++20 sources: CoverSeed (`cover_seed`), VertexSeed (`vertex_seed`), KernelEnum (`subsumed_search`, `kernel_enum_bk_pivot.hpp`), the maximal-clique index (`maximal_clique_index`) and the maintenance engine (`maintenance_engine`) |
| `tests/` | Tests run by `make check` |
| `preprocessing/` | Converter from KONECT edge lists to the input format |
| `scripts/` | Runs of CKM for Section V |

## Build

Requirements: Linux, GCC 13 or newer and GNU make.

```bash
make            # builds build/ckm
make check      # optional: runs the tests
```

## Datasets

The experiments use six timestamped graphs from KONECT (http://konect.cc/).

| Name in the paper | KONECT network | Edges used |
|---|---|---|
| cit-dblp | `dblp-cite` | all |
| dblp-coauthor | `dblp_coauthor` | all |
| flickr-growth | `flickr-growth` | first 1M |
| wikipedia-growth | `wikipedia-growth` | first 12M |
| youtube | `youtube-u-growth` | all |
| stack-overflow | `sx-stackoverflow` | first 1M |

Convert a downloaded `out.<name>` file with

```bash
cd preprocessing
make
./preprocess_dynamic_graph out.flickr-growth    # writes processed_flickr-growth.txt
```

The converter keeps the two endpoint columns and the timestamp column, drops
records without a timestamp, and maps the timestamps to consecutive integers
in increasing order. `ckm` then treats the graph as undirected, drops
self-loops, keeps the earliest timestamp of each duplicated edge, and
processes the edges in timestamp order (ties in file order). `--max-edges=N`
keeps the first `N` edges of that order.

## Running

```bash
./build/ckm processed_dblp-cite.txt --batch-size=1000 --stats
```

maintains the maximal cliques of `cit-dblp` with CKM in batches of 1000
insertions and writes the changes of every batch to standard output, one line
per changed maximal clique: the batch's timestamp, `new` or `del`, and its
vertices. `--no-output` drops these lines, and `--stats` prints the counts
and the time of every step to standard error. `./build/ckm --help` lists all
options.

| Option | Meaning |
|---|---|
| `--update-mode=insert\|reverse-delete\|sliding-window` | insertion, deletion or hybrid batches (default `insert`) |
| `--batch-size=N` | `N` edges per batch (default 1000) |
| `--max-edges=N`, `--max-batches=N` | use only the first edges or batches |
| `--maintenance=full\|new-only` | maintain every maximal clique, or only find the new maximal cliques of insertion batches (default `full`) |
| `--new-clique-search=coverseed\|vertexseed` | the search for new maximal cliques (default `coverseed`) |
| `--change-output=new\|subsumed\|both` | which changes are written |
| `--measured-time-limit=SECONDS` | stop after the measured batches have run this long (exit status 124) |

### Workloads

| Paper | Option |
|---|---|
| Insertion batches | `--update-mode=insert`: the batches are inserted into the empty graph in timestamp order |
| Deletion batches | `--update-mode=reverse-delete`: starting from the graph with all edges, the batches are deleted in reverse order |
| Hybrid batches | `--update-mode=sliding-window`: the window first holds the first half of the edges, and each batch inserts the next `b` edges and deletes the `b` oldest edges |

CoverSeed always uses Anchor-or-Pivot branching. Anchor branching and pivot
branching are both implemented in full in `src/cover_seed.cpp`, and the
choice between them is made in two places: `CoverSeed::searchSeed` at the
root of a seed's search and `CoverSeed::branchUnanchored` below it. Fixing
that choice in both places runs either branching alone.

## Experiments

Each script in `scripts/` runs CKM for one part of Section V on the processed
files in `data/` and writes a table of one row per run, with the log of every
run next to it. The other methods compared in the paper are not included.

| Section | Script | Runs | Output |
|---|---|---|---|
| V-B | `maintenance.sh` | CKM under insertion, deletion and hybrid batches | `results/maintenance/maintenance.tsv` |
| V-C | `new_cliques.sh` | CoverSeed and VertexSeed | `results/new_cliques/new_cliques.tsv` |
| V-D | `subsumed.sh` | KernelEnum | `results/subsumed/subsumed.tsv` |

```bash
bash scripts/maintenance.sh
DATA_DIR=/path/to/data ONLY="cit-dblp" REPEATS=5 bash scripts/new_cliques.sh
```

`DATA_DIR`, `OUT`, `CKM`, `REPEATS` and `ONLY` (a subset of the datasets)
override the defaults. The paper runs every setting five times
(`REPEATS=5`) and reports the average time. Runs are stopped after 5 hours
of measured batches in V-B and V-D, and after 1 hour (VertexSeed) or 2 hours
(CoverSeed) in V-C. On a cluster, submit each script as a job of its own.
