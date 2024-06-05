# Variables
CC = gcc
CFLAGS = -Wall -g
LDFLAGS = -lm

# Source files
SRCS = main.c
# Object files
OBJS = $(SRCS:.c=.o)
# Executable file
TARGET = main

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compiling source files to object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets
.PHONY: all clean

