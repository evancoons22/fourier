#include <stdio.h>
#include <math.h>

typedef struct {
    double real;
    double imag;
} complex;

#define PI 3.14159265

void DFT(complex *x, complex *X, int N) {
    // the outer loop is the number of resulting frequencies we want
    for (int k = 0; k < N; k++) {
        X[k].real = 0.0;
        X[k].imag = 0.0;
        // the inner loop calculates each of those
        for (int n = 0; n < N; n++) {
            // use euler's formula to calculate exp(ix) = cos(x) + isin(x)
            double angle = 2 * PI * k * n / N;
            X[k].real += x[n].real * cos(angle);
            X[k].imag -= x[n].imag * sin(angle);
        }
    }
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
            y[k].real += (X[n].real * cos(angle) - X[n].imag*sin(angle));
            y[k].imag +=  (X[n].real * cos(angle) + X[n].imag*sin(angle));
        }
        y[k].real /= N;
        y[k].imag /= N;
    }
} 

double test_func(double x) { 
    //return 5 * cos(2*x / 5 + 2) + 3 * sin(2 *x);
    return 5 * cos(x) + 3 * sin(2*x);
} 


int main() { 
    int N = 3;
    complex X[N];
    complex x[N];

    for (int i = 0; i < N; i ++) { 
        x[i].real += test_func(i);
        x[i].imag = 0.0;
    }

    complex y[N];

    DFT(x, X, N);
    iDFT(X,y,N);

    for (int i = 0; i < N; i++) { 
        printf("original[%d]  = %.5f + %.5fi\n", i, x[i].real, x[i].imag);
    } 

    for (int i = 0; i < N; i++) { 
        printf("sample[%d]  = %.5f + %.5fi\n", i, X[i].real, X[i].imag);
    } 

    /* desmos printing
    printf("printing desmos version\n");
    double period;
    for (int i = 0; i < N; i++) { 
        period = 2 * PI / i;
        printf("+ %.5f * cos(%.5fx)", X[i].real, period);
    } 
    */

    // reconstruction
    for (int i = 0; i < N; i++) { 
        printf("reconstructed sample[%d]  = %.5f + %.5fi\n", i, y[i].real, y[i].imag);
    } 

    printf("printing desmos version\n");
    double period;
    for (int i = 0; i < N; i++) { 
        period = 2 * PI / i;
        printf("+ %.5f * cos(%.5fx)", y[i].real, period);
    } 
} 
