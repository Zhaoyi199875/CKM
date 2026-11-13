# CKM: Maximal Clique Maintenance in Dynamic Graphs

This repository contains the code for the paper
*Beyond One Edge at a Time: Maximal Clique Maintenance in Dynamic Graphs*.

CKM (Cover and Kernel Maintenance) maintains the maximal cliques of a dynamic
graph under batches of edge insertions, deletions or both. After each
insertion batch, CoverSeed (Section III) finds the new maximal cliques and
KernelEnum (Section IV) finds the subsumed cliques. Deletion and hybrid
batches are handled as described in Section II-D. The binary `ckm`
implements CKM.

## Build

Requirements: Linux, GCC 13 or newer and GNU make.

```bash
make            # builds build/ckm
```

`ckm` reads edge lists with timestamps. Convert a KONECT edge list
(http://konect.cc/), downloaded as `out.<name>`, with

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

maintains the maximal cliques of the graph in `processed_dblp-cite.txt` with
CKM in batches of 1000 insertions and writes the changes of every batch to
standard output, one line per changed maximal clique: the batch's timestamp,
`new` or `del`, and its vertices. `--no-output` drops these lines, and
`--stats` prints the counts and the time of every step to standard error.
`./build/ckm --help` lists all options.

| Option | Meaning |
|---|---|
| `--update-mode=insert\|reverse-delete\|sliding-window` | insertion, deletion or hybrid batches (default `insert`) |
| `--batch-size=N` | `N` edges per batch (default 1000) |
| `--max-edges=N`, `--max-batches=N` | use only the first edges or batches |
| `--change-output=new\|subsumed\|both` | which changes are written |
| `--measured-time-limit=SECONDS` | stop after the measured batches have run this long (exit status 124) |

### Workloads

| Paper | Option |
|---|---|
| Insertion batches | `--update-mode=insert`: the batches are inserted into the empty graph in timestamp order |
| Deletion batches | `--update-mode=reverse-delete`: starting from the graph with all edges, the batches are deleted in reverse order |
| Hybrid batches | `--update-mode=sliding-window`: the window first holds the first half of the edges, and each batch inserts the next `b` edges and deletes the `b` oldest edges |
