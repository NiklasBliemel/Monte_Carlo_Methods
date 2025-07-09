CXX = g++
CXXFLAGS = -O3 -std=c++11 -I.

SRC_DIR := .
OBJ_DIR := build
BIN_DIR := bin

PROGRAMS := sheet_1 sheet_2 sheet_3 sheet_4_1 sheet_4_2 sheet_4_3 sheet_5
OBJS_COMMON := Grid.o Ising.o XY.o Routines.o

# Full paths
OBJECTS_COMMON := $(addprefix $(OBJ_DIR)/,$(OBJS_COMMON))
BINARIES := $(addprefix $(BIN_DIR)/,$(PROGRAMS))

all: $(BINARIES)

# Rule for building each program
$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(OBJECTS_COMMON) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile each .cpp file into .o file in build/
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Ensure directories exist
$(OBJ_DIR):
	mkdir -p $@

$(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
