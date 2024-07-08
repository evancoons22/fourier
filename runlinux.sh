#/bin/bash
gcc fft.c -lm -g -o fft && ./fft && python3 plot.py && feh plot.png
