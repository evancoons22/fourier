# Makefile to compile fft.c with math and portaudio libraries

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -O3

# Libraries to link
LIBS = -lm -lportaudio -lncurses

# Source file
SRC = synth.c

# Output executable
TARGET = synth

# Build rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# Clean rule
clean:
	rm -f $(TARGET)

