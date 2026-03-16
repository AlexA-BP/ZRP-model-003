# Makefile

# Default build: debug
BUILD ?= debug

SRC_DIR = ./src/cpp_version
BUILD_DIR = ./build
BIN_DIR = ./bin
INCLUDE_DIR = ./include
LIB_DIR = ./lib

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
OBJS_WITHOUT_SRCS = $(BUILD_DIR/%.o)
TARGET_DEBUG = $(BIN_DIR)/debug
TARGET_RELEASE = $(BIN_DIR)/release

# predefined targets use all caps in their name
TARGET = $(TARGET_$(shell echo '$(BUILD)' | tr '[:lower:]' '[:upper:]'))

CXX = g++

CXXFLAGS_COMMON = -Wall -Wextra -std=c++20 -I$(INCLUDE_DIR) -L$(LIB_DIR)

# Linker Flags
LDLFLAGS =-Wl,-rpath=$(LIB_DIR) \
				-Wl,--enable-new-dtags -g -fsanitize=address

LIBS = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lz -ldl -lm 

CXXFLAGS_DEBUG = -g -O0
CXXFLAGS_RELEASE = -O2 -DNDEBUG

CXXFLAGS = $(CXXFLAGS_COMMON) $(CXXFLAGS_$(BUILD))

# create directories if they don't exist
$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

$(TARGET): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(LDLFLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean clean-all

clean:
	$(RM) -r $(BUILD_DIR)/*
	$(RM) -r $(BIN_DIR)/*

.PHONY: all
all: $(TARGET)
