CXX ?= g++
BUILD_DIR ?= build
CPPFLAGS ?= -Isrc
CXXFLAGS ?= -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDFLAGS ?=

SOURCES := src/insertion_stream.cpp src/incremental_graph.cpp \
	src/maximal_clique_index.cpp src/cover_seed.cpp src/subsumed_search.cpp \
	src/maintenance_engine.cpp src/main.cpp
OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
CKM := $(BUILD_DIR)/ckm

.PHONY: all clean

all: $(CKM)

$(CKM): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

-include $(OBJECTS:.o=.d)
