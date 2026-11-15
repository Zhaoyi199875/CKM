CXX ?= g++
BUILD_DIR ?= build
CPPFLAGS ?= -Isrc
CXXFLAGS ?= -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDFLAGS ?=
PYTHON ?= python3

SOURCES := src/insertion_stream.cpp src/incremental_graph.cpp \
	src/maximal_clique_index.cpp src/cover_seed.cpp src/subsumed_search.cpp \
	src/maintenance_engine.cpp src/main.cpp
OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
CKM := $(BUILD_DIR)/ckm
TEST_SEARCH := $(BUILD_DIR)/test_search
TEST_OBJECT := $(BUILD_DIR)/test_search.o

.PHONY: all check clean

all: $(CKM)

$(CKM): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TEST_SEARCH): $(filter-out $(BUILD_DIR)/main.o,$(OBJECTS)) $(TEST_OBJECT)
	$(CXX) $(LDFLAGS) $^ -o $@

$(TEST_OBJECT): tests/test_search.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

check: $(CKM) $(TEST_SEARCH)
	$(TEST_SEARCH)
	$(PYTHON) tests/test_cli.py $(CKM)

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

-include $(OBJECTS:.o=.d) $(TEST_OBJECT:.o=.d)
