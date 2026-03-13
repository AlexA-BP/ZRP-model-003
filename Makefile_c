# Makefile

# Default build: debug
BUILD ?= debug

SRC_DIR = ./src/C
BUILD_DIR = ./build
BIN_DIR = ./bin
INCLUDE_DIR = ./include
LIB_DIR = ./lib

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
OBJS_WITHOUT_SRCS = $(BUILD_DIR/%.o)
TARGET_DEBUG = $(BIN_DIR)/debug
TARGET_RELEASE = $(BIN_DIR)/release

# predefined targets use all caps in their name
TARGET = $(TARGET_$(shell echo '$(BUILD)' | tr '[:lower:]' '[:upper:]'))

CC = gcc

CFLAGS_COMMON = -Wall -Wextra -std=c99 -I$(INCLUDE_DIR) -L$(LIB_DIR)

# Linker Flags
LDLFLAGS =-Wl,-rpath=$(LIB_DIR) \
				-Wl,--enable-new-dtags -g -fsanitize=address

LIBS = -lhdf5_hl -lhdf5 -lz -ldl -lm -lpcg_random

CFLAGS_DEBUG = -g -O0
CFLAGS_RELEASE = -O2 -DNDEBUG

CFLAGS = $(CFLAGS_COMMON) $(CFLAGS_$(BUILD))

# create directories if they don't exist
$(shell mkdir -p $(BUILD_DIR) $(BIN_DIR))

$(TARGET): $(OBJS) $(BUILD_DIR)/entropy.o 
	$(CC) $(CFLAGS) $(LDLFLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean clean-all

clean:
	$(shell find $(BUILD_DIR) ! -name "entropy.o" -type f -exec rm -f {} + )
	$(RM) -r $(BIN_DIR)/*

.PHONY: all
all: $(TARGET)
