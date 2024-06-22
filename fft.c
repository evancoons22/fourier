#include <stdio.h> 
#include <math.h>
#include <stdlib.h>

unsigned int bit_reversal(unsigned int x,  int bits) {
    // want to go from 00011 to 11000, for example
    unsigned int y = 0;
    for (int i = 0; i < bits; i++) { 
        y = y << 1 | (x & 1);
        x >>= 1;
    } 
    return y;
} 

void bit_reverse_copy(double *a, double *b, int l) { 
    // for each number in the a array, reverse the bits
    int bits = log2(l);
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

void print_array(int *a, int n) { 
    for (int i = 0; i < n; i ++) { 
        printf("%d ", a[i]);
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


int main() { 
    double arr1[8] = {0,1,2,3,4,5,6,7};

    unsigned int n = 8;

    double* arr2 = (double*)malloc(n * sizeof(double) * 8);
    bit_reverse_copy(arr1, arr2, n);


}
