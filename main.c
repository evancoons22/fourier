#include <math.h>
#include <stdio.h>
// #include "omp.h"

typedef struct {
    double real;
    double imag;
} complex;

#define PI 3.14159265

void DFT(complex *x, complex *X, int N, int K) {
    // the outer loop is the number of resulting frequencies we want
    for (int k = 0; k < K; k++) {
        X[k].real = 0.0;
        X[k].imag = 0.0;
        // the inner loop calculates each of those
        for (int n = 0; n < N; n++) {
            // use euler's formula to calculate exp(ix) = cos(x) + isin(x)
            // X[k].real += x[n].real * cos(angle);
            // X[k].imag -= x[n].imag * sin(angle);
            double angle = 2 * PI * k * n / N;
            double cos_angle = cos(angle);
            double sin_angle = sin(angle);
            X[k].real += x[n].real * cos_angle - x[n].imag * sin_angle;
            X[k].imag += x[n].real * sin_angle + x[n].imag * cos_angle;
        }
    }
}

void FFT(complex *x, complex *X, int N) { 
}

void iDFT(complex *X, complex *y, int N) {
    // computing the inverse DFT
    for (int k = 0; k < N; k++) {
        y[k].real = 0.0;
        y[k].imag = 0.0;
        // the inner loop calculates each of those
        for (int n = 0; n < N; n++) {
            // use euler's formula to calculate exp(ix) = cos(x) + isin(x)
            double angle = 2 * PI * k * n / N;
            y[k].real += (X[n].real * cos(angle) - X[n].imag * sin(angle));
            y[k].imag += (X[n].real * cos(angle) + X[n].imag * sin(angle));
        }
        y[k].real /= N;
        y[k].imag /= N;
    }
}

// assuming integer frequencies 1 to N
void print_desmos(complex *X, int N) {
    printf("printing desmos version\n");
    double period;
    for (int i = 0; i < N; i++) {
        period = 2 * PI / i;
        printf("+ %.5f * cos(%.5fx)", X[i].real, period);
    }
}

double test_func(double x) {
    // return 5 * cos(2*x / 5 + 2) + 3 * sin(2 *x);
    // return cos(PI * x) + cos(PI * 2 * x);
    // the period is 2pi/b, the frequency is b/2pi
    return cos(10 * 2 * PI * x);
}

void write_frequencies(complex *X, int K) {
    // write to data.txt
    FILE *f = fopen("data.txt", "w");
    for (int k = 0; k < K; k++) {
        // the indices are i * range / N
        fprintf(f, "%d %f %f\n", k, X[k].real, X[k].imag);
    }
}

int main() {
    int K = 30;
    int N = 100;
    int range = 100;
    complex X[K];
    complex x[N];

    for (int i = 0; i < N; i++) {
        //double index = (double)i / N * range;
        // double index = (double)i;
        x[i].real = test_func(i);
        x[i].imag = 0.0;
    }

    // x -> X -> y
    DFT(x, X, N, K);
    // iDFT(X, y, N);

    write_frequencies(X, K);

    // for (int i = 0; i < N; i++) {
    //     printf("original[%d]  = %.5f + %.5fi\n", i, x[i].real, x[i].imag);
    // }

    // for (int i = 0; i < N; i++) {
    //     printf("at frequency:[%d], real amplitude:   = %.5f + %.5fi\n", i,
    //     X[i].real, X[i].imag);
    // }

    /* desmos printing
       printf("printing desmos version\n");
       double period;
       for (int i = 0; i < N; i++) {
       period = 2 * PI / i;
       printf("+ %.5f * cos(%.5fx)", X[i].real, period);
       }
       */

    // reconstruction
    // for (int i = 0; i < N; i++) {
    //     printf("reconstructed sample[%d]  = %.5f + %.5fi\n", i, y[i].real,
    //     y[i].imag);
    // }

    // print_desmos(X, N);
    return 0;
}
