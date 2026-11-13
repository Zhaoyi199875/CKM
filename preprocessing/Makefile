CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Wconversion -Wshadow

TARGET := preprocess_dynamic_graph
SOURCE := preprocess_dynamic_graph.cpp

.PHONY: all clean sanitize

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) $(SOURCE) -o $(TARGET)

sanitize: CXXFLAGS := -std=c++20 -O1 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -fsanitize=address,undefined -fno-omit-frame-pointer
sanitize: $(SOURCE)
	$(CXX) $(CXXFLAGS) $(SOURCE) -fsanitize=address,undefined -o $(TARGET)_sanitize

clean:
	$(RM) $(TARGET) $(TARGET)_sanitize
