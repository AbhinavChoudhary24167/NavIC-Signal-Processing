#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RSP.h"

#define PI 3.14159265358979323846

void RSP(double mag_out_rd[NFFT][NUM_FD], double* peak_mag, int* peak_range_idx, int* peak_doppler_idx,
         Complex fft_prn[NFFT], Complex in_signal[NFFT], int Nsamp, int nfft, double fd[NUM_FD], int num_fd, double delta) {
   
    // Temporary arrays for FFT processing
    Complex sgl_dp[NFFT];
    Complex fft_sgl_dp[NFFT];
    Complex product[NFFT];
    Complex ifft_result[NFFT];
   
    // Process Doppler bins
    for (int f = 0; f < num_fd; f++) {
        // Apply Doppler compensation
        for (int k = 0; k < Nsamp; k++) {
            double phase = -2 * PI * fd[f] * k * delta;
            Complex exp_factor = { cos(phase), sin(phase) };
            sgl_dp[k].real = in_signal[k].real * exp_factor.real - in_signal[k].imag * exp_factor.imag;
            sgl_dp[k].imag = in_signal[k].real * exp_factor.imag + in_signal[k].imag * exp_factor.real;
        }
       
        // Zero-padding for FFT
        for (int k = Nsamp; k < nfft; k++) {
            sgl_dp[k].real = 0.0;
            sgl_dp[k].imag = 0.0;
        }
       
        // Compute FFT
        for (int k = 0; k < nfft; k++) {
            fft_sgl_dp[k] = sgl_dp[k];
        }
        cfft(fft_sgl_dp, nfft);
       
        // Matched Filtering: Multiply FFT_sgl_dp with **conjugate** of fft_prn
        for (int k = 0; k < nfft; k++) {
            product[k].real = fft_sgl_dp[k].real * fft_prn[k].real + fft_sgl_dp[k].imag * fft_prn[k].imag;
            product[k].imag = fft_sgl_dp[k].imag * fft_prn[k].real - fft_sgl_dp[k].real * fft_prn[k].imag;
        }
       
        // Compute IFFT
        for (int k = 0; k < nfft; k++) {
            ifft_result[k] = product[k];
        }
        cifft(ifft_result, nfft);
       
        // Compute magnitude and store result
        for (int k = 0; k < nfft; k++) {
            mag_out_rd[k][f] = sqrt(ifft_result[k].real * ifft_result[k].real + ifft_result[k].imag * ifft_result[k].imag);
        }
    }
   
    // Peak search with LOCAL MAXIMUM validation
    *peak_mag = 0.0;
    *peak_range_idx = -1;
    *peak_doppler_idx = -1;

    for (int f = 0; f < num_fd; f++) {
        for (int k = 1; k < nfft - 1; k++) {  // Ensure it's not at the edges
            if (mag_out_rd[k][f] > mag_out_rd[k-1][f] && mag_out_rd[k][f] > mag_out_rd[k+1][f]) {
                if (mag_out_rd[k][f] > *peak_mag) {
                    *peak_mag = mag_out_rd[k][f];
                    *peak_range_idx = k;
                    *peak_doppler_idx = f;
                }
            }
        }
    }
}
