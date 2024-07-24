#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double real;
    double imag;
} complex;

#define PI 3.14159265

void DFT(complex *x, complex *X, int N, int K, int RANGE) {
    //  // M is the new length after zero-padding
    //  edit here for sampling frequency
    int M = 1 * N;
    
    // Allocate memory for the zero-padded signal
    complex *x_padded = (complex *)calloc(M, sizeof(complex));
    
    // Copy the original signal into the zero-padded signal
    for (int n = 0; n < N; n++) {
        x_padded[n] = x[n];
    }

    for (int k = 0; k < K; k++) {
        X[k].real = 0.0;
        X[k].imag = 0.0;
        // the inner loop calculates each of those
        for (int n = 0; n < N; n++) {
            // use euler's formula to calculate exp(ix) = cos(x) + isin(x)
            // in this convention, exp(-ix) for DFT
            double angle = 2 * PI * k * n / M;
            double cos_angle = cos(angle);
            double sin_angle = sin(angle);
            X[k].real += x_padded[n].real * cos_angle + x_padded[n].imag * sin_angle;
            X[k].imag += x_padded[n].real * sin_angle - x_padded[n].imag * cos_angle;
        }
        //X[k].real /= N;
        //X[k].imag /= N;
    }
    free(x_padded);
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
            double cos_angle = cos(angle);
            double sin_angle = sin(angle);
            y[k].real += (X[n].real * cos_angle - X[n].imag * sin_angle); 
            y[k].imag += (X[n].real * cos_angle + X[n].imag * sin_angle);
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
    // the period is 2pi/b, the frequency is b/2pi
    // return cos(10 * 2 * PI * x) + 2 * cos(20 * 2 * PI * x) + 3 * cos(30 * 2 * PI * x) + 4 * cos(40 * 2 * PI * x);  // frequency of 10
    return cos(10 * 2 * PI * x) + 2 * cos(20 * 2 * PI * x);  // frequency of 10
}

double test_func2(double x) { 
    int w = PI * 2 * 10;
    int total = 0.0; 
    for (int i = 1; i < 5; i ++) { 
        total += i*cos(w * x * i);
    } 
    return total;
} 

double test_func3(double x) { 
    double w = PI * 2;
    return 2 * sin(w*x) + sin(w*50*x);
}

void write_frequencies(complex *X, int K) {
    // write to data.txt
    FILE *f = fopen("data.txt", "w");
    for (int k = 0; k < K; k++) {
        // the indices are i * range / N
        fprintf(f, "%d %f %f\n", k, X[k].real, X[k].imag);
    }
}

void write_frequencies_csv(complex *X, int K) {
    // write to data.csv
    FILE *f = fopen("data.csv", "w");
    // write columns, index, real, imag
    fprintf(f, "index,real,imag\n");
    for (int k = 0; k < K; k++) {
        // the indices are i * range / N
        fprintf(f, "%d,%f,%f\n", k, X[k].real, X[k].imag);
    }
    
}

#define  K 64
#define  N 64
#define  RANGE  100

int main() {
    complex X[K];
    complex x[N];
    
    // malloc for dynamic size
    //complex *X = (complex *)malloc(K * sizeof(complex));
    //complex *x = (complex *)malloc(N * sizeof(complex));

    // gathering samples
    // printf("original signal \n");
    for (int i = 0; i < N; i++) {
        //double index = (double)i / N * range;
        // double index = (double)i;
        x[i].real = test_func((double)i/N);
        x[i].imag = 0.0;
        // printf("%f %f\n", x[i].real, x[i].imag);
    }

    // x -> X -> y
    DFT(x, X, N, K, RANGE);

    complex y[N];
    iDFT(X, y, N);

    // printf("reconstructed signal \n");
   // for (int i = 0; i < N; i++) {
   //     printf("%f %f\n", y[i].real, y[i].imag);
   // }

    write_frequencies_csv(X, K);

    // complex y[N];
    // iDFT(X, y, N);
    // print_desmos(y, N);


    return 0;
}
