# Makefile to compile fft.c with math and portaudio libraries

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -O3

# Libraries to link
LIBS = -lm -lportaudio

# Source file
SRC = fft.c

# Output executable
TARGET = fft

# Build rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LIBS)

# Clean rule
clean:
	rm -f $(TARGET)

