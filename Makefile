# for time measurements
CXX_CMD = clang++-12 -std=c++17 -O3 -Wall -Wextra -Wshadow -I include

# for debug
# CXX_CMD = clang++-10 -std=c++17 -O0 -g -fsanitize=address,undefined -Wall -Wextra -Wshadow -I include


BUILD_DIR = build

GENERATORS = $(patsubst src/gen/%.cpp, $(BUILD_DIR)/gen_%, $(wildcard src/gen/*.cpp))
ALGORITHMS = $(patsubst src/algorithms/%.cpp, $(BUILD_DIR)/algo_%, $(wildcard src/algorithms/*.cpp))

all: $(GENERATORS) $(ALGORITHMS)

$(GENERATORS): $(BUILD_DIR)/gen_%: src/gen/%.cpp | $(BUILD_DIR)
	$(CXX_CMD) $^ -o $@

$(ALGORITHMS): $(BUILD_DIR)/algo_%: src/algo_main.cpp src/algorithms/%.cpp | $(BUILD_DIR)
	$(CXX_CMD) $^ -o $@

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

.PHONY: all clean
