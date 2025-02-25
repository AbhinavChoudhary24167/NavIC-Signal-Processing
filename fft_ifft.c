#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RSP.h"

#define PI 3.14159265358979323846

// Swap two Complex numbers
void swap(Complex* a, Complex* b) {
    Complex temp = *a;
    *a = *b;
    *b = temp;
}

// Bit-reversal permutation
void bit_reversal(Complex* x, int n) {
    int j = 0;
    for (int i = 1; i < n; i++) {
        int bit = n >> 1;
        while (j >= bit) {
            j -= bit;
            bit >>= 1;
        }
        j += bit;
        if (i < j) {
            swap(&x[i], &x[j]);
        }
    }
}

// Complex FFT using Cooley-Tukey algorithm
void cfft(Complex *x, int n) {
    bit_reversal(x, n);
    for (int len = 2; len <= n; len *= 2) {
        double angle = -2.0 * PI / len;
        Complex wlen = {cos(angle), sin(angle)};
        for (int i = 0; i < n; i += len) {
            Complex w = {1.0, 0.0};
            for (int j = 0; j < len / 2; j++) {
                Complex u = x[i + j];
                Complex v;
                v.real = w.real * x[i + j + len / 2].real - w.imag * x[i + j + len / 2].imag;
                v.imag = w.real * x[i + j + len / 2].imag + w.imag * x[i + j + len / 2].real;
                x[i + j].real = u.real + v.real;
                x[i + j].imag = u.imag + v.imag;
                x[i + j + len / 2].real = u.real - v.real;
                x[i + j + len / 2].imag = u.imag - v.imag;
                Complex temp;
                temp.real = w.real * wlen.real - w.imag * wlen.imag;
                temp.imag = w.real * wlen.imag + w.imag * wlen.real;
                w = temp;
            }
        }
    }
}

// Complex IFFT (Inverse FFT)
void cifft(Complex *x, int n) {
    for (int i = 0; i < n; i++) {
        x[i].imag = -x[i].imag;
    }
    cfft(x, n);
    for (int i = 0; i < n; i++) {
        x[i].real /= n;
        x[i].imag = -x[i].imag / n;
    }
}
