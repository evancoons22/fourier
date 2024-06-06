#include <stdio.h> 
#include <math.h>

unsigned int bit_reversal(unsigned int x,  int bits) {
    // want to go from 00011 to 11000, for example
    int y = 0;
    for (int i = 0; i < bits; i++) { 
        y = y << 1 | (x & 1);
        x >>= 1;
    } 
    return y;
} 

void bit_reverse_copy(int *a, int *b, int l) { 
    // for each number in the a array, reverse the bits
    int bits = log(l);
    for (int i = 0; i < l; i++) { 
        b[i] = bit_reversal(a[i], bits);
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


int main() { 
    int num = 6;
    int y = bit_reversal(num, sizeof(int)*4);
    //printf("%lu\n", 4 * sizeof(int));
    // print_binary(num);
    // printf("\n");
    // print_binary(y);

    int arr1[6] = {1,2,3,4,5,7};
    int arr2[6] = {0,0,0,0,0,0};
    print_binary_array(arr1, 6);
    bit_reverse_copy(arr1, arr2, 6);
    printf("\n");
    print_binary_array(arr2, 6);

} 
