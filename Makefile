BINARIES=preprocess_data pseudoalign k_index_build

SRC_DIR=src
BIN_DIR=bin

.PHONY: $(BINARIES)
all: $(BINARIES)

preprocess_data:
	$(CXX) -std=c++17 -I ~/include -L ~/lib $(SRC_DIR)/preprocess_data.cpp -O3 -o $(BIN_DIR)/preprocess_data -lstdc++fs -g

pseudoalign:
	$(CXX) -std=c++17 -I ~/include -I r-index/internal -L ~/lib $(SRC_DIR)/pseudoalign.cpp -O3 -o $(BIN_DIR)/pseudoalign -lsdsl -ldivsufsort -ldivsufsort64 -g

k_index_build:
	$(CXX) -std=c++17 -I ~/include -I r-index/internal -L ~/lib $(SRC_DIR)/k_index_build.cpp -O3 -o $(BIN_DIR)/k_index_build -lsdsl -ldivsufsort -ldivsufsort64 -g