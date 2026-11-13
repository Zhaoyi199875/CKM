# CKM: Maximal Clique Maintenance in Dynamic Graphs

This repository contains the code for the paper
*Beyond One Edge at a Time: Maximal Clique Maintenance in Dynamic Graphs*.

CKM (Cover and Kernel Maintenance) maintains the maximal cliques of a dynamic
graph under batches of edge insertions, deletions or both. After each
insertion batch, CoverSeed (Section III) finds the new maximal cliques and
KernelEnum (Section IV) finds the subsumed cliques. Deletion and hybrid
batches are handled as described in Section II-D.

## Build

Requirements: Linux, GCC 13 or newer and GNU make.

To build CKM, run `make` in the repository root. This creates the binary
`build/ckm`, and `make clean` removes it again.

```bash
make
```

`ckm` reads a text file with one edge per line: two integer vertex ids and an
integer timestamp, separated by spaces or tabs. Text after `#` is ignored.
The graph is undirected. Self-loops are dropped, a repeated edge keeps its
earliest timestamp, and the edges are processed in timestamp order, ties in
file order.

To convert an edge list whose lines hold the two endpoints first and the
timestamp last, build the converter in `preprocessing/` and run it on the
file:

```bash
cd preprocessing
make
./preprocess_dynamic_graph example_dataset
```

This writes `processed_example_dataset` next to the input file. The converter
keeps the first two columns and the last one, skips lines that start with `%`
or `#`, drops edges whose timestamp is 0, and maps the timestamps to
consecutive integers in increasing order.

## Running

```bash
./build/ckm processed_example_dataset --batch-size=1000 --stats
```

inserts the edges of `processed_example_dataset` in batches of 1000 and
maintains the maximal cliques after every batch. It prints the changed
maximal cliques of every batch and, because of `--stats`, the counts and the
running time of every step at the end. A changed maximal clique is printed
as, for example,

```
5	new	1 2 3
5	del	1 2
```

which says that the batch with timestamp 5 created the maximal clique
{1, 2, 3} and removed {1, 2}. A vertex without edges is a maximal clique of
one vertex, and each vertex is in the graph from just before the batch of its
first edge.

| Option | Meaning |
|---|---|
| `--update-mode=insertion\|deletion\|hybrid` | the kind of batches, see below (default `insertion`) |
| `--batch-size=N` | `N` edges per batch (default 1000) |
| `--batch-by-timestamp` | one batch per timestamp instead |
| `--max-edges=N` | use only the first `N` edges |
| `--max-batches=N` | use only the first `N` batches |
| `--change-output=new\|removed\|both` | write the new maximal cliques, the removed ones, or both (default `both`) |
| `--output=FILE` | write the changes to `FILE` instead of standard output |
| `--no-output` | write no changes |
| `--stats` | print counts and times |

The three update modes are the workloads of the paper:

| Mode | Batches |
|---|---|
| `insertion` | the batches are inserted into the empty graph in timestamp order |
| `deletion` | starting from the graph with all edges, the batches are deleted in reverse order |
| `hybrid` | a window first holds the first half of the edges, and each batch then inserts the next `b` edges and deletes the `b` oldest ones |

In the `deletion` and `hybrid` modes, the graph the batches start from is
built first and is not counted in the times or the output.
