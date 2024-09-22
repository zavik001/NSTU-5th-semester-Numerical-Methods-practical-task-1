# Executable target files directory
BUILD_DIR = build
TARGET_FLOAT = $(BUILD_DIR)/ldlt_float.exe
TARGET_DOUBLE = $(BUILD_DIR)/ldlt_double.exe
TARGET_FLOAT_DOUBLE = $(BUILD_DIR)/ldlt_float_double.exe

# Compiler and flags
CXX = g++
CXXFLAGS_FLOAT = -DFF -Iinclude
CXXFLAGS_DOUBLE = -DDD -Iinclude
CXXFLAGS_LONG_DOUBLE = -DFD -Iinclude

# Directories with source and header files
SRC_DIR = src
INCLUDE_DIR = include

# Source and header files
SRC = $(SRC_DIR)/SLAUSolverLDLT.cpp Main.cpp

# Default build rule (create build directory and compile all versions)
all: $(BUILD_DIR) $(TARGET_FLOAT) $(TARGET_DOUBLE) $(TARGET_FLOAT_DOUBLE)
	@echo "Executables built. Use 'make runFloat', 'make runDouble', or 'make runFloatDouble' to run the respective version."

# Create build directory if it doesn't exist
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Rule for creating the float version
$(TARGET_FLOAT): $(SRC)
	@echo "Building float version..."
	$(CXX) $(CXXFLAGS_FLOAT) -o $@ $(SRC)

# Rule for creating the double version
$(TARGET_DOUBLE): $(SRC)
	@echo "Building double version..."
	$(CXX) $(CXXFLAGS_DOUBLE) -o $@ $(SRC)

# Rule for creating the float double version
$(TARGET_FLOAT_DOUBLE): $(SRC)
	@echo "Building float double version..."
	$(CXX) $(CXXFLAGS_LONG_DOUBLE) -o $@ $(SRC)

# Run the float version
runFloat: $(TARGET_FLOAT)
	@echo "Running float version..."
	./$(TARGET_FLOAT)

# Run the double version
runDouble: $(TARGET_DOUBLE)
	@echo "Running double version..."
	./$(TARGET_DOUBLE)

# Run the float double version
runFloatDouble: $(TARGET_FLOAT_DOUBLE)
	@echo "Running float double version..."
	./$(TARGET_FLOAT_DOUBLE)

# Clean up the build directory and executables
clean:
	@echo "Cleaning up..."
	@rm -rf $(BUILD_DIR)

.PHONY: all clean runFloat runDouble runFloatDouble
