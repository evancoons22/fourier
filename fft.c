#include <stdio.h> 
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

typedef struct {
    double real;
    double imag;
} complex;

// Utility functions for complex numbers
complex complex_add(complex a, complex b) {
    complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

complex complex_sub(complex a, complex b) {
    complex result;
    result.real = a.real - b.real;
    result.imag = a.imag - b.imag;
    return result;
}

complex complex_mul(complex a, complex b) {
    complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

complex complex_from_polar(double r, double theta) {
    complex result;
    result.real = r * cos(theta);
    result.imag = r * sin(theta);
    return result;
}

double test_func(double x) {
    // the period is 2pi/b, the frequency is b/2pi
    return cos(10 * 2 * PI * x) + 2 * cos(20 * 2 * PI * x) + 3 * cos(30 * 2 * PI * x) + 4 * cos(40 * 2 * PI * x);  // frequency of 10
}

void write_frequencies_csv(complex *X, int K) {
    // write to data.csv
    FILE *f = fopen("data.csv", "w");
    // write columns, index, real, imag
    fprintf(f, "index,real,imag\n");
    for (int k = 0; k < K/2; k++) {
        // the indices are i * range / N
        fprintf(f, "%d,%f,%f\n", k, sqrt(X[k].real*X[k].real + X[k].imag *  X[k].imag), 0.0);
        // fprintf(f, "%d,%f,%f\n", k, X[k].real, 0.0);
    }
}

unsigned int bit_reversal(unsigned int x,  int bits) {
    // want to go from 00011 to 11000, for example
    unsigned int y = 0;
    for (int i = 0; i < bits; i++) { 
        y = y << 1 | (x & 1);
        x >>= 1;
    } 
    return y;
} 

void bit_reverse_copy(complex *a, complex *b, int l) { 
    // for each number in the a array, reverse the bits
    // copy into array b
    int bits = (int)log2(l);
    for (unsigned int i = 0; i < l; i++) { 
        unsigned int revk = bit_reversal(i, bits);
        b[revk] = a[i];
    } 
} 

void print_binary(int x) { 
    // int is 16 bits? 
    for (int i = 1 << 15; i > 0; i = i >> 1) { 
        if ((i & x) != 0) { 
            printf("1");
        } else { 
            printf("0");
        } 
    } 
} 

void print_array(double *a, int n) { 
    for (int i = 0; i < n; i ++) { 
        printf("%f ", a[i]);
    } 

} 

void print_binary_array(int *a, int n) { 
    for (int i = 0; i < n; i ++) { 
        print_binary(a[i]);
        printf(" ");
    } 
} 

int is_power_of_two(int x) {
    return (x > 0) && ((x & (x-1)) == 0);
}  

void cooley_turkey_rec(complex *x, complex *X, int n) { 
    // if the length is one, return the term itself
    if (n == 1) {  
        X[0] = x[0];
        return;
    } 

    complex even[n/2];
    complex odd[n/2];
    
    int even_index = 0, odd_index = 0;

    // split into even and odd problems
    for (size_t i = 0; i < n; i++) { 
        if (i %2 == 0) { 
            even[even_index++] = x[i];
        } else { 
            odd[odd_index++] = x[i];
        } 
    } 

    // create 2 arrays for the transformed inputs

    complex even_transformed[n/2];
    complex odd_transformed[n/2];

    // call recursion again
    cooley_turkey_rec(even, even_transformed, n/2);
    cooley_turkey_rec(odd, odd_transformed, n/2);

    // combine 2 arrays
    for (size_t k = 0; k < n/2; k++) { 
        double angle = -2 * PI * k / n;
        complex twiddle = complex_from_polar(1.0, angle);
        complex t = complex_mul(twiddle, odd_transformed[k]);

        X[k] = complex_add(even_transformed[k], t);
        X[k + n / 2] = complex_sub(even_transformed[k], t);
    } 

} 

void cooley_turkey(complex *x, complex *X, int n) { 
    // y  is the bit reverse version of x
    //
    if (!is_power_of_two(n)) { 
        fprintf(stderr, "Error: FFT length must be a power of 2\n");
        return;
    } 

    complex y[n];
    // double *y = (double *)malloc(n * sizeof(double));
    bit_reverse_copy(x, y, n); 

    cooley_turkey_rec(y, X, n);
} 

int main() { 
    complex arr1[128];
    unsigned int n = 128;
    //// sample function for testing
    for (int i = 0; i < n; i++) {
        arr1[i].real= test_func((double)i/n);
        arr1[i].imag = 0.0;
    }

    complex X[n];

    cooley_turkey(arr1, X, n);

    // Print the results
    // for (int i = 0; i < n; i++) {
        // printf("X[%d] = %.2f + %.2fi\n", i, X[i].real, X[i].imag);
    // }
    
    write_frequencies_csv(X, n);

    return 0;

}
