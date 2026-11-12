CXX ?= g++
BUILD_DIR ?= build
CPPFLAGS ?= -Isrc
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDFLAGS ?=

CORE_SOURCES := src/insertion_stream.cpp src/incremental_graph.cpp \
	src/maximal_clique_index.cpp src/new_clique_search.cpp \
	src/star_clique_search.cpp src/subsumed_clique_search.cpp \
	src/vertex_clique_search.cpp src/somei.cpp src/maintenance_engine.cpp
APP_SOURCES := $(CORE_SOURCES) src/main.cpp
APP_OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(APP_SOURCES))
TEST_OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(CORE_SOURCES)) \
	$(BUILD_DIR)/test_search.o
SOMEI_TEST_OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(CORE_SOURCES)) \
	$(BUILD_DIR)/test_somei.o

.PHONY: all check check-somei sanitize clean

all: $(BUILD_DIR)/immce_star

$(BUILD_DIR)/immce_star: $(APP_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/test_search: $(TEST_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/test_somei: $(SOMEI_TEST_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/test_search.o: tests/test_search.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/test_somei.o: tests/test_somei.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

# The SOMEi baseline (src/somei.hpp): the paper's worked example, every
# four- and five-vertex old/new edge state, and random multi-batch streams,
# each against brute force and against CM-Constructor run from scratch; then
# the binary on the semantics fixture, under both update scopes, with its
# change stream held to the star/residual-kernel one.  The untimed build by
# enumeration (--somei-build-update=bk) is held to the paper's construction:
# the same state in test_somei, the same change streams after it on
# tests/somei_build_fixture.txt.
check-somei: $(BUILD_DIR)/test_somei $(BUILD_DIR)/immce_star
	$(BUILD_DIR)/test_somei
	for pruning in on off; do \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-by-timestamp --framework=somei --maintenance=full \
			--somei-pruning=$$pruning --change-output=both --no-output \
			--stats 2>&1 | \
			grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=6, current-maximal-cliques=2' || \
			{ echo "somei pruning=$$pruning failed"; exit 1; }; \
	done
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=somei --maintenance=full --somei-changes=off \
		--somei-profile=on --no-output --stats 2>&1 | \
		grep -q 'subsumed-maximal-cliques=0, current-maximal-cliques=2'
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=somei --maintenance=full --change-output=both \
		--output=$(BUILD_DIR)/somei_changes.txt
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=star --maintenance=full --maintenance-pipeline=batched \
		--subsumed-algorithm=residual-kernel --change-output=both \
		--output=$(BUILD_DIR)/star_changes.txt
	cmp -s $(BUILD_DIR)/somei_changes.txt $(BUILD_DIR)/star_changes.txt || \
		{ echo "somei change stream differs from star/residual-kernel"; exit 1; }
	rm -f $(BUILD_DIR)/somei_changes.txt $(BUILD_DIR)/star_changes.txt
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=somei \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=somei \
		--maintenance=full --maintenance-pipeline=streaming \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=somei \
		--maintenance=full --somei-changes=off > /dev/null 2>&1
	# The untimed build by Bron-Kerbosch on G[R] (2026-11-06) leaves the
	# state the paper's construction leaves: the deletions and the mixed
	# batches after it list the same changes.
	for mode in reverse-delete random-labels; do \
		for update in bk paper; do \
			$(BUILD_DIR)/immce_star tests/somei_build_fixture.txt \
				--batch-size=50 --framework=somei --maintenance=full \
				--update-mode=$$mode --somei-build=one-batch \
				--somei-build-update=$$update --change-output=both \
				--output=$(BUILD_DIR)/somei_build_$$update.txt --stats \
				2> $(BUILD_DIR)/somei_build_$$update.err || exit 1; \
			grep -q "build=one-batch, build-update=$$update" \
				$(BUILD_DIR)/somei_build_$$update.err || \
				{ echo "somei $$mode build-update=$$update not reported"; exit 1; }; \
		done; \
		cmp -s $(BUILD_DIR)/somei_build_bk.txt $(BUILD_DIR)/somei_build_paper.txt || \
			{ echo "somei $$mode: the enumerated build changes the change stream"; exit 1; }; \
	done
	rm -f $(BUILD_DIR)/somei_build_bk.* $(BUILD_DIR)/somei_build_paper.*
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=star \
		--somei-build-update=bk --no-output > /dev/null 2>&1

# The semantics fixture has four new maximal cliques when cut by timestamp
# (two in one batch) and, under full maintenance, six subsumed ones over two
# current cliques: two cliques of two, and the single-vertex cliques {10},
# {20}, {30}, {40}, each subsumed by its vertex's first edge (a vertex that
# appears is an isolated vertex of the graph before its batch, see
# MaintenanceEngine).  Every check below has to reproduce those numbers: both
# star methods, the edge baseline, the vertex baseline as producer under both
# pipelines, the cross-framework comparisons, both maintenance pipelines, and
# each subsumed framework in turn (the residual-kernel size gate on and off)
# -- the loop below is where a new one gets added.  Reverse delete runs the
# same batches backwards, where every count swaps; random labels mixes
# insertions and deletions in one batch, where every configuration must
# write the same change stream; a sliding window keeps the newest edges,
# where every configuration must write the change stream worked out by hand.
check: check-somei $(BUILD_DIR)/test_search $(BUILD_DIR)/immce_star
	$(BUILD_DIR)/test_search
	for scope in unbound root; do for filter in on off; do \
	for branching in rule anchor pivot; do \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-by-timestamp --framework=star --star-branch-scope=$$scope \
			--star-branching=$$branching \
			--star-root-filter=$$filter --no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
			{ echo "star branch-scope=$$scope branching=$$branching root-filter=$$filter by timestamp failed"; exit 1; }; \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--framework=star --star-branch-scope=$$scope \
			--star-branching=$$branching \
			--star-root-filter=$$filter --no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=2, subsumed-maximal-cliques=0' || \
			{ echo "star branch-scope=$$scope branching=$$branching root-filter=$$filter in one batch failed"; exit 1; }; \
	done; done; done
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--batch-by-timestamp --framework=edge --no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
		{ echo "edge by timestamp failed"; exit 1; }
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--batch-by-timestamp --framework=vertex --no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
		{ echo "vertex by timestamp failed"; exit 1; }
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--framework=vertex --no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=2, subsumed-maximal-cliques=0' || \
		{ echo "vertex in one batch failed"; exit 1; }
	# The vertex branching ablation: every branching finds the same cliques,
	# and the pivot is the default.
	for branching in pivot new-edges rule-root rule-unbound; do \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-by-timestamp --framework=vertex --vertex-branching=$$branching \
			--no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
			{ echo "vertex branching=$$branching by timestamp failed"; exit 1; }; \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--framework=vertex --vertex-branching=$$branching \
			--no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=2, subsumed-maximal-cliques=0' || \
			{ echo "vertex branching=$$branching in one batch failed"; exit 1; }; \
	done
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=vertex \
		--no-output --stats 2>&1 | grep -q 'branching=pivot' || \
		{ echo "vertex branching does not default to pivot"; exit 1; }
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --vertex-branching=maybe --no-output > /dev/null 2>&1
	# The seed-cover ablation: searching from the cover vertices alone finds the
	# same cliques under every branching, and every endpoint is the default.
	for branching in pivot new-edges rule-root rule-unbound; do \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-by-timestamp --framework=vertex --vertex-seeds=cover \
			--vertex-branching=$$branching --no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
			{ echo "vertex seeds=cover branching=$$branching by timestamp failed"; exit 1; }; \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--framework=vertex --vertex-seeds=cover --vertex-branching=$$branching \
			--no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=2, subsumed-maximal-cliques=0' || \
			{ echo "vertex seeds=cover branching=$$branching in one batch failed"; exit 1; }; \
	done
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=vertex \
		--no-output --stats 2>&1 | grep -q 'seeds=endpoints' || \
		{ echo "vertex seeds do not default to endpoints"; exit 1; }
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --vertex-seeds=some --no-output > /dev/null 2>&1
	for pipeline in batched streaming; do \
		$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-by-timestamp --framework=vertex --maintenance=full \
			--maintenance-pipeline=$$pipeline --change-output=both \
			--no-output --stats 2>&1 | \
			grep -q 'subsumed-maximal-cliques=6, current-maximal-cliques=2' || \
			{ echo "vertex full maintenance/$$pipeline failed"; exit 1; }; \
	done
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=both --no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0'
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=all --no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
		{ echo "framework=all failed"; exit 1; }
	# Reverse delete: the fixture's batches inserted, then deleted again, last
	# first.  Deleting undoes inserting, so the run destroys the four cliques
	# the insertions created, creates the six they subsumed and ends with no
	# edge and four single-vertex cliques -- and its change stream is the
	# insertions' with new and del swapped.
	for framework in star edge vertex; do \
		for algorithm in optimized optimized-bitmap residual-kernel; do \
			$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
				--batch-by-timestamp --framework=$$framework --maintenance=full \
				--update-mode=reverse-delete --subsumed-algorithm=$$algorithm \
				--no-output --stats 2>&1 | \
				grep -q 'new-maximal-cliques=6, subsumed-maximal-cliques=4, current-maximal-cliques=4' || \
				{ echo "reverse delete $$framework/$$algorithm failed"; exit 1; }; \
		done; \
	done
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=somei --maintenance=full --update-mode=reverse-delete \
		--no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=6, subsumed-maximal-cliques=4, current-maximal-cliques=4' || \
		{ echo "reverse delete somei failed"; exit 1; }
	# SOMEi's build in one batch (the measured deletions unchanged), and the
	# progress lines; --somei-build belongs to SOMEi alone.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=somei --maintenance=full --update-mode=reverse-delete \
		--somei-build=one-batch --progress --no-output --stats \
		> $(BUILD_DIR)/somei_one_batch.txt 2>&1 || exit 1
	grep -q 'new-maximal-cliques=6, subsumed-maximal-cliques=4, current-maximal-cliques=4' \
		$(BUILD_DIR)/somei_one_batch.txt || \
		{ echo "reverse delete somei with a one-batch build failed"; exit 1; }
	grep -q 'build=one-batch' $(BUILD_DIR)/somei_one_batch.txt || \
		{ echo "somei does not report its one-batch build"; exit 1; }
	grep -q '^progress: measured-begin' $(BUILD_DIR)/somei_one_batch.txt && \
		grep -q '^progress: measured-batches=' $(BUILD_DIR)/somei_one_batch.txt || \
		{ echo "--progress printed no progress lines"; exit 1; }
	rm -f $(BUILD_DIR)/somei_one_batch.txt
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=star \
		--somei-build=one-batch --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=somei \
		--maintenance=full --somei-build=maybe --no-output > /dev/null 2>&1
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=star --maintenance=full --maintenance-pipeline=batched \
		--subsumed-algorithm=residual-kernel --change-output=both \
		--output=$(BUILD_DIR)/forward.txt || exit 1
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=star --maintenance=full --update-mode=reverse-delete \
		--subsumed-algorithm=residual-kernel --change-output=both \
		--output=$(BUILD_DIR)/reverse.txt || exit 1
	awk -F'\t' -v OFS='\t' '{ print $$1, ($$2 == "new" ? "del" : "new"), $$3 }' \
		$(BUILD_DIR)/reverse.txt | LC_ALL=C sort > $(BUILD_DIR)/reverse.sorted
	LC_ALL=C sort $(BUILD_DIR)/forward.txt > $(BUILD_DIR)/forward.sorted
	cmp -s $(BUILD_DIR)/forward.sorted $(BUILD_DIR)/reverse.sorted || \
		{ echo "reverse delete is not the insertions' change stream reversed"; exit 1; }
	rm -f $(BUILD_DIR)/forward.txt $(BUILD_DIR)/reverse.txt \
		$(BUILD_DIR)/forward.sorted $(BUILD_DIR)/reverse.sorted
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--update-mode=reverse-delete --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--update-mode=reverse-delete --maintenance=full \
		--maintenance-pipeline=streaming --no-output > /dev/null 2>&1
	# Random labels: each edge labelled insert or delete; the delete-labelled
	# ones are built first, then every batch inserts and deletes at once.
	# Every producer and subsumed algorithm under both strategies, and SOMEi,
	# must write the star/residual-kernel change stream.
	for seed in 1 2 3; do \
		run() { $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-size=2 --maintenance=full --update-mode=random-labels \
			--label-seed=$$seed --change-output=both "$$@"; }; \
		run --framework=star --subsumed-algorithm=residual-kernel \
			--output=$(BUILD_DIR)/mixed_reference.txt || exit 1; \
		for framework in star edge vertex; do \
			for algorithm in optimized optimized-bitmap residual-kernel; do \
				for strategy in single-pass two-phase; do \
					run --framework=$$framework --subsumed-algorithm=$$algorithm \
						--mixed-strategy=$$strategy \
						--output=$(BUILD_DIR)/mixed.txt || exit 1; \
					cmp -s $(BUILD_DIR)/mixed_reference.txt $(BUILD_DIR)/mixed.txt || \
						{ echo "random labels seed $$seed $$framework/$$algorithm/$$strategy differs"; exit 1; }; \
				done; \
			done; \
		done; \
		run --framework=somei --output=$(BUILD_DIR)/mixed.txt || exit 1; \
		cmp -s $(BUILD_DIR)/mixed_reference.txt $(BUILD_DIR)/mixed.txt || \
			{ echo "random labels seed $$seed somei differs"; exit 1; }; \
		run --framework=somei --somei-build=one-batch \
			--output=$(BUILD_DIR)/mixed.txt || exit 1; \
		cmp -s $(BUILD_DIR)/mixed_reference.txt $(BUILD_DIR)/mixed.txt || \
			{ echo "random labels seed $$seed somei with a one-batch build differs"; exit 1; }; \
	done
	rm -f $(BUILD_DIR)/mixed_reference.txt $(BUILD_DIR)/mixed.txt
	# Every edge labelled insert is plain insertion; every edge labelled
	# delete builds the whole stream and empties it again, first batch first,
	# leaving the four vertices as single-vertex cliques.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=2 \
		--framework=star --maintenance=full --maintenance-pipeline=batched \
		--subsumed-algorithm=residual-kernel --change-output=both \
		--output=$(BUILD_DIR)/inserted.txt || exit 1
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=2 \
		--framework=star --maintenance=full --update-mode=random-labels \
		--insert-fraction=1 --subsumed-algorithm=residual-kernel \
		--change-output=both --output=$(BUILD_DIR)/labelled.txt || exit 1
	cmp -s $(BUILD_DIR)/inserted.txt $(BUILD_DIR)/labelled.txt || \
		{ echo "random labels with every edge inserted is not insertion"; exit 1; }
	rm -f $(BUILD_DIR)/inserted.txt $(BUILD_DIR)/labelled.txt
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=2 \
		--maintenance=full --update-mode=random-labels --insert-fraction=0 \
		--no-output --stats 2>&1 | grep -q 'current-maximal-cliques=4,' || \
		{ echo "random labels with every edge deleted did not empty the graph"; exit 1; }
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--update-mode=random-labels --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--update-mode=random-labels --maintenance=full \
		--maintenance-pipeline=streaming --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --maintenance=full \
		--mixed-strategy=two-phase --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --maintenance=full \
		--update-mode=random-labels --insert-fraction=1.5 --no-output > /dev/null 2>&1
	# Sliding window: the fixture's four edges in stream order are 10-20,
	# 20-30, 30-40, 10-30, and the window defaults to half of them.  The
	# first two are built; then each batch of one inserts the next edge and
	# deletes the oldest.  The first batch leaves 10 without an edge and gives
	# 40 its first; the second gives 10 an edge again and leaves 20 without
	# one.  Every producer and subsumed algorithm under both strategies, and
	# SOMEi, must write that change stream.
	printf '7\tnew\t30 40\n7\tnew\t10\n7\tdel\t10 20\n7\tdel\t40\n8\tnew\t10 30\n8\tnew\t20\n8\tdel\t20 30\n8\tdel\t10\n' \
		> $(BUILD_DIR)/window_expected.txt
	run() { $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--batch-size=1 --maintenance=full --update-mode=sliding-window \
		--change-output=both "$$@"; }; \
	for framework in star edge vertex; do \
		for algorithm in optimized optimized-bitmap residual-kernel; do \
			for strategy in single-pass two-phase; do \
				run --framework=$$framework --subsumed-algorithm=$$algorithm \
					--mixed-strategy=$$strategy \
					--output=$(BUILD_DIR)/window.txt || exit 1; \
				cmp -s $(BUILD_DIR)/window_expected.txt $(BUILD_DIR)/window.txt || \
					{ echo "sliding window $$framework/$$algorithm/$$strategy differs"; exit 1; }; \
			done; \
		done; \
	done; \
	run --framework=somei --output=$(BUILD_DIR)/window.txt || exit 1; \
	cmp -s $(BUILD_DIR)/window_expected.txt $(BUILD_DIR)/window.txt || \
		{ echo "sliding window somei differs"; exit 1; }; \
	run --framework=somei --somei-build=one-batch \
		--output=$(BUILD_DIR)/window.txt || exit 1; \
	cmp -s $(BUILD_DIR)/window_expected.txt $(BUILD_DIR)/window.txt || \
		{ echo "sliding window somei with a one-batch build differs"; exit 1; }
	rm -f $(BUILD_DIR)/window_expected.txt $(BUILD_DIR)/window.txt
	# One batch of two: the same window, one batch -- 10 keeps an edge, so
	# its single-vertex clique never forms.  An explicit window of three
	# deletes only 10-20; cut to its first three edges the run's window is one
	# edge and ends holding 30-40 and the single-vertex cliques of 10 and 20.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=2 \
		--maintenance=full --update-mode=sliding-window --no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=3, subsumed-maximal-cliques=3, current-maximal-cliques=3,' || \
		{ echo "sliding window in one batch failed"; exit 1; }
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--maintenance=full --update-mode=sliding-window --window-edges=3 \
		--no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=1, subsumed-maximal-cliques=1, current-maximal-cliques=3,' || \
		{ echo "sliding window of three failed"; exit 1; }
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--max-edges=3 --maintenance=full --update-mode=sliding-window \
		--no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=4, current-maximal-cliques=3,' || \
		{ echo "sliding window over the first three edges failed"; exit 1; }
	# --max-edges cuts the stream as --max-batches does when they agree.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--max-batches=3 --framework=star --maintenance=full \
		--change-output=both --output=$(BUILD_DIR)/cut_batches.txt || exit 1
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--max-edges=3 --framework=star --maintenance=full \
		--change-output=both --output=$(BUILD_DIR)/cut_edges.txt || exit 1
	cmp -s $(BUILD_DIR)/cut_batches.txt $(BUILD_DIR)/cut_edges.txt || \
		{ echo "--max-edges does not cut the stream as --max-batches does"; exit 1; }
	rm -f $(BUILD_DIR)/cut_batches.txt $(BUILD_DIR)/cut_edges.txt
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=2 \
		--maintenance=full --update-mode=sliding-window --window-edges=1 \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--maintenance=full --update-mode=sliding-window --window-edges=4 \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --maintenance=full \
		--window-edges=2 --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--maintenance=full --update-mode=sliding-window --warmup-edges=1 \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--maintenance=full --update-mode=sliding-window --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--maintenance=full --update-mode=sliding-window --insert-fraction=0.5 \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-size=1 \
		--max-edges=5 --no-output > /dev/null 2>&1
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=star --maintenance=full --maintenance-pipeline=streaming \
		--change-output=both --no-output --stats 2>&1 | \
		grep -q 'subsumed-maximal-cliques=6, current-maximal-cliques=2'
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --batch-by-timestamp \
		--framework=star --maintenance=full --maintenance-pipeline=batched \
		--change-output=both --no-output --stats 2>&1 | \
		grep -q 'subsumed-maximal-cliques=6, current-maximal-cliques=2'
	for algorithm in optimized optimized-bitmap residual-kernel; do \
		for pipeline in batched streaming; do \
			$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
				--batch-by-timestamp --framework=star --maintenance=full \
				--maintenance-pipeline=$$pipeline \
				--subsumed-algorithm=$$algorithm \
				--change-output=both --no-output --stats 2>&1 | \
				grep -q 'subsumed-maximal-cliques=6, current-maximal-cliques=2' || \
				{ echo "subsumed $$algorithm/$$pipeline failed"; exit 1; }; \
		done; \
	done
	# The residual-kernel size gate, on and off.
	for histogram in on off; do \
		for pipeline in batched streaming; do \
			$(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
				--batch-by-timestamp --framework=star --maintenance=full \
				--maintenance-pipeline=$$pipeline \
				--subsumed-algorithm=residual-kernel \
				--subsumed-size-histogram=$$histogram \
				--change-output=both --no-output --stats 2>&1 | \
				grep -q 'subsumed-maximal-cliques=6, current-maximal-cliques=2' || \
				{ echo "residual-kernel size-histogram=$$histogram/$$pipeline failed"; exit 1; }; \
		done; \
	done
	# The warm-up: cut into the run's own batches it leaves the change stream
	# exactly as it is without one and the statistics count only the batches
	# after it; cut differently it still leaves the same maximal cliques.  The
	# statistics go to a file: `grep -q` on a pipe would close it at the first
	# match, and the run it kills then never writes out its --output file.
	for pipeline in batched streaming; do \
		run() { $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
			--batch-size=1 --framework=star --maintenance=full \
			--maintenance-pipeline=$$pipeline \
			--subsumed-algorithm=residual-kernel --change-output=both "$$@"; }; \
		run --output=$(BUILD_DIR)/warmup_off.txt || exit 1; \
		run --warmup-edges=2 --output=$(BUILD_DIR)/warmup_on.txt --stats \
			2> $(BUILD_DIR)/warmup_on.log || exit 1; \
		grep -q 'batches=2, available-batches=2, max-batches=all, warmup-edges=2, warmup-batches=2, warmup-batch-size=1,' \
			$(BUILD_DIR)/warmup_on.log || \
			{ echo "warm-up statistics/$$pipeline failed"; exit 1; }; \
		cmp -s $(BUILD_DIR)/warmup_off.txt $(BUILD_DIR)/warmup_on.txt || \
			{ echo "warm-up change stream/$$pipeline failed"; exit 1; }; \
		plain=$$(run --no-output --stats 2>&1 | sed -n 's/^graph:.*current-maximal-cliques=\([0-9]*\).*/\1/p'); \
		coarse=$$(run --warmup-edges=3 --warmup-batch-size=3 --no-output --stats 2>&1 | sed -n 's/^graph:.*current-maximal-cliques=\([0-9]*\).*/\1/p'); \
		test -n "$$plain" && test "$$plain" = "$$coarse" || \
			{ echo "warm-up final state/$$pipeline failed ($$plain vs $$coarse)"; exit 1; }; \
		! run --warmup-edges=5 --no-output > /dev/null 2>&1 || \
			{ echo "warm-up longer than the stream was accepted"; exit 1; }; \
	done
	rm -f $(BUILD_DIR)/warmup_off.txt $(BUILD_DIR)/warmup_on.txt $(BUILD_DIR)/warmup_on.log
	# Unbound has been the default scope since 2026-10-08.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=star \
		--no-output --stats 2>&1 | grep -q 'branch-scope=unbound' || \
		{ echo "star branch scope does not default to unbound"; exit 1; }
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --star-branch-scope=maybe --no-output > /dev/null 2>&1
	# The root filter is on unless the ablation turns it off.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=star \
		--no-output --stats 2>&1 | grep -q 'root-filter=on' || \
		{ echo "star root filter is not on by default"; exit 1; }
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --star-root-filter=maybe --no-output > /dev/null 2>&1
	# The rule branches the star unless the 2026-10-31 ablation fixes a partition.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --framework=star \
		--no-output --stats 2>&1 | grep -q 'branching=rule' || \
		{ echo "star branching does not default to the rule"; exit 1; }
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --star-branching=maybe --no-output > /dev/null 2>&1
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --tie-order=file --no-output --stats 2>&1 | grep -q "tie-order=file"
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --tie-order=vertex --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --tie-order=shuffle --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --tie-seed=1 --no-output > /dev/null 2>&1
	# --measured-time-limit and --measured-time-extend (2026-11-06) are stated
	# on the graph line, and a run that finishes inside them ends normally;
	# zero, an extension without a limit and one below it are refused.
	$(BUILD_DIR)/immce_star tests/insertion_semantics.txt --no-output --stats 2>&1 | \
		grep -q 'result-mode=[a-z-]*, measured-time-limit=none, measured-time-extend=none$$'
	$(BUILD_DIR)/immce_star tests/somei_build_fixture.txt --batch-size=50 \
		--framework=somei --maintenance=full --update-mode=reverse-delete \
		--somei-build=one-batch --measured-time-limit=600 --progress \
		--no-output --stats 2>&1 | grep -q 'measured-time-limit=600, measured-time-extend=none$$'
	$(BUILD_DIR)/immce_star tests/somei_build_fixture.txt --batch-size=50 \
		--framework=somei --maintenance=full --update-mode=reverse-delete \
		--somei-build=one-batch --measured-time-limit=600 \
		--measured-time-extend=1200 --progress --no-output --stats 2>&1 | \
		grep -q 'measured-time-limit=600, measured-time-extend=1200$$'
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --measured-time-limit=0 \
		--no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt \
		--measured-time-extend=600 --no-output > /dev/null 2>&1
	! $(BUILD_DIR)/immce_star tests/insertion_semantics.txt --measured-time-limit=600 \
		--measured-time-extend=300 --no-output > /dev/null 2>&1
	@echo "check: ok"

sanitize:
	$(MAKE) BUILD_DIR=build-sanitize \
		CXXFLAGS='-std=c++20 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -fsanitize=address,undefined -fno-omit-frame-pointer' \
		LDFLAGS='-fsanitize=address,undefined' check

clean:
	rm -rf build build-sanitize

-include $(APP_OBJECTS:.o=.d) $(BUILD_DIR)/test_search.d \
	$(BUILD_DIR)/test_somei.d
