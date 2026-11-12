"""Writes tests/somei_build_fixture.txt (2026-11-06): a small edge stream for
check-somei's comparison of SOMEi's two untimed builds (--somei-build-update
=bk|paper).  Cliques of 3-10 vertices planted around a few low-numbered hubs,
plus random edges, in a shuffled order with tied timestamps -- so that in
vertex-id order the hubs head many bodies and the build's steps have R of
every size.  Deterministic; run from the project directory:

    python3 tests/make_somei_build_fixture.py > tests/somei_build_fixture.txt
"""
import random
import sys

generator = random.Random(20261106)
vertices = list(range(1, 91))
hubs = vertices[:8]
others = vertices[8:]
edges = set()
for _ in range(30):
    size = generator.randint(3, 10)
    members = generator.sample(hubs, generator.randint(1, 2))
    members += generator.sample(others, size - len(members))
    for i, u in enumerate(members):
        for v in members[i + 1:]:
            edges.add((min(u, v), max(u, v)))
while len(edges) < 900:
    u, v = generator.sample(vertices, 2)
    edges.add((min(u, v), max(u, v)))
ordered = sorted(edges)
generator.shuffle(ordered)
out = sys.stdout
out.write("# SOMEi build fixture: planted cliques around low-numbered hubs plus\n")
out.write("# random edges (tests/make_somei_build_fixture.py).\n")
for position, (u, v) in enumerate(ordered):
    if generator.random() < 0.5:
        u, v = v, u
    out.write(f"{u} {v} {position // 3}\n")
