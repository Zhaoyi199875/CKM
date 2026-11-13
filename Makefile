CXX ?= g++
BUILD_DIR ?= build
CPPFLAGS ?= -Isrc
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDFLAGS ?=

CORE_SOURCES := src/insertion_stream.cpp src/incremental_graph.cpp \
	src/maximal_clique_index.cpp src/cover_seed.cpp src/vertex_seed.cpp \
	src/subsumed_search.cpp src/maintenance_engine.cpp
APP_SOURCES := $(CORE_SOURCES) src/main.cpp
APP_OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(APP_SOURCES))
CORE_OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(CORE_SOURCES))
CKM := $(BUILD_DIR)/ckm
FIXTURE := tests/insertion_semantics.txt

.PHONY: all check sanitize clean

all: $(CKM)

$(CKM): $(APP_OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/test_search: $(CORE_OBJECTS) $(BUILD_DIR)/test_search.o
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/%.o: tests/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

check: $(BUILD_DIR)/test_search $(CKM)
	# The unit tests: every search against brute force on all small graphs,
	# random streams and structured wide cases.
	$(BUILD_DIR)/test_search
	# The semantics fixture has four new maximal cliques when cut by
	# timestamp and two in one batch.  Both new-clique searches find them.
	for args in --new-clique-search=coverseed --new-clique-search=vertexseed; do \
		$(CKM) $(FIXTURE) --batch-by-timestamp --maintenance=new-only $$args \
			--no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=0' || \
			{ echo "$$args by timestamp failed"; exit 1; }; \
		$(CKM) $(FIXTURE) --maintenance=new-only $$args --no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=2, subsumed-maximal-cliques=0' || \
			{ echo "$$args in one batch failed"; exit 1; }; \
	done
	# Full maintenance by timestamp: four new and six subsumed cliques, two
	# left -- two cliques of two, and the single-vertex cliques {10}, {20},
	# {30}, {40}, each subsumed by its vertex's first edge.  VertexSeed in
	# place of CoverSeed writes the same change stream.
	$(CKM) $(FIXTURE) --batch-by-timestamp --change-output=both \
		--output=$(BUILD_DIR)/changes_ckm.txt --stats 2> $(BUILD_DIR)/changes_ckm.err
	grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=6, current-maximal-cliques=2' \
		$(BUILD_DIR)/changes_ckm.err
	$(CKM) $(FIXTURE) --batch-by-timestamp --change-output=both \
		--new-clique-search=vertexseed --output=$(BUILD_DIR)/changes.txt
	cmp -s $(BUILD_DIR)/changes_ckm.txt $(BUILD_DIR)/changes.txt || \
		{ echo "VertexSeed's change stream differs from CoverSeed's"; exit 1; }
	# Deletions: the fixture's batches inserted untimed, then deleted again,
	# last first.  Deleting undoes inserting, so the run destroys the four
	# cliques the insertions created, creates the six they subsumed and ends
	# with four single-vertex cliques.
	for args in --new-clique-search=coverseed --new-clique-search=vertexseed; do \
		$(CKM) $(FIXTURE) --batch-by-timestamp --update-mode=reverse-delete $$args \
			--no-output --stats 2>&1 | \
			grep -q 'new-maximal-cliques=6, subsumed-maximal-cliques=4, current-maximal-cliques=4' || \
			{ echo "reverse delete $$args failed"; exit 1; }; \
	done
	# Hybrid batches: the fixture's four edges in stream order are 10-20,
	# 20-30, 30-40, 10-30, and the window holds half of them.  The first two
	# are built; then each batch of one inserts the next edge and deletes
	# the oldest.  The first batch leaves 10 without an edge and gives 40 its
	# first; the second gives 10 an edge again and leaves 20 without one.
	printf '7\tnew\t30 40\n7\tnew\t10\n7\tdel\t10 20\n7\tdel\t40\n8\tnew\t10 30\n8\tnew\t20\n8\tdel\t20 30\n8\tdel\t10\n' \
		> $(BUILD_DIR)/window_expected.txt
	for args in --new-clique-search=coverseed --new-clique-search=vertexseed; do \
		$(CKM) $(FIXTURE) --batch-size=1 --update-mode=sliding-window \
			--change-output=both $$args --output=$(BUILD_DIR)/window.txt || exit 1; \
		cmp -s $(BUILD_DIR)/window_expected.txt $(BUILD_DIR)/window.txt || \
			{ echo "sliding window $$args differs"; exit 1; }; \
	done
	$(CKM) $(FIXTURE) --batch-size=2 --update-mode=sliding-window --no-output \
		--stats 2>&1 | \
		grep -q 'new-maximal-cliques=3, subsumed-maximal-cliques=3, current-maximal-cliques=3,'
	$(CKM) $(FIXTURE) --batch-size=1 --max-edges=3 --update-mode=sliding-window \
		--no-output --stats 2>&1 | \
		grep -q 'new-maximal-cliques=4, subsumed-maximal-cliques=4, current-maximal-cliques=3,'
	# The defaults, and what the statistics report.
	$(CKM) $(FIXTURE) --no-output --stats 2>&1 | grep -q \
		'new-clique-search=coverseed, maintenance=full'
	$(CKM) $(FIXTURE) --maintenance=new-only --new-clique-search=vertexseed \
		--no-output --stats 2>&1 | grep -q '^vertexseed: time-ms='
	$(CKM) $(FIXTURE) --update-mode=reverse-delete --measured-time-limit=600 \
		--progress --no-output --stats 2>&1 | grep -q 'measured-time-limit=600$$'
	# Refused: an unknown search, a hybrid window cut by timestamp, both
	# outputs, and a time limit of zero.
	! $(CKM) $(FIXTURE) --new-clique-search=unknown --no-output > /dev/null 2>&1
	! $(CKM) $(FIXTURE) --batch-by-timestamp --update-mode=sliding-window --no-output > /dev/null 2>&1
	! $(CKM) $(FIXTURE) --output=$(BUILD_DIR)/x.txt --no-output > /dev/null 2>&1
	! $(CKM) $(FIXTURE) --measured-time-limit=0 --no-output > /dev/null 2>&1
	rm -f $(BUILD_DIR)/changes*.txt $(BUILD_DIR)/changes_ckm.err $(BUILD_DIR)/window*.txt
	@echo "check: ok"

sanitize:
	$(MAKE) BUILD_DIR=build-sanitize \
		CXXFLAGS='-std=c++20 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -fsanitize=address,undefined -fno-omit-frame-pointer' \
		LDFLAGS='-fsanitize=address,undefined' check

clean:
	rm -rf build build-sanitize

-include $(APP_OBJECTS:.o=.d) $(BUILD_DIR)/test_search.d
