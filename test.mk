SRC_DIR = ./src
BUILD_DIR = ./build
BIN_DIR = ./bin
INCLUDE_DIR = ./include
LIB_DIR = ./lib

FILE = savetofile

RAW_SRCS = $(FILE).c simulation.c
SRCS = $(addprefix $(SRC_DIR)/,$(RAW_SRCS))
OBJS_WITH_SRC = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
RAW_OBJS_WITHOUT_SRC = entropy.o
OBJS_WITHOUT_SRCS = $(addprefix $(BUILD_DIR)/,$(RAW_OBJS_WITHOUT_SRC))
TARGET = $(BIN_DIR)/$(FILE)

CC = gcc
CFLAGS_COMMON = -Wall -Wextra -std=c99 -I$(INCLUDE_DIR) -L$(LIB_DIR)
CFLAGS_DEBUG = -g -O0
CFLAGS = $(CFLAGS_COMMON) $(CFLAGS_DEBUG)

LIBS = -lm -lpcg_random


$(TARGET): $(OBJS_WITH_SRC) $(OBJS_WITHOUT_SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJS_WITH_SRC): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: test