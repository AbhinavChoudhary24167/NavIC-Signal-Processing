#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdbool.h>
#include <float.h>
#include "navic_input.h"
#include "in_data5.h"

/* System parameters */
#define C 300000000        // Speed of light in m/s
#define FS_TX 1023000      // Transmitter sampling frequency
#define FS_RX 4000000      // Receiver sampling frequency
#define C_LEN 1023         // PRN code length
#define SATELLITE_NO 7     // Number of satellites
#define T 0.001            // Time period
#define NFFT 4096          // FFT size
#define DOPPLER_COUNT 25   // Number of Doppler frequency bins
#define MAX_CA_CODE_SIZE 4000  // Maximum size of resampled CA code
#define PI 3.14159265358979323846  // Pi value
#define detection_threshold 3.0  // Threshold for satellite detection

/* Doppler frequency axis */
const double fd[DOPPLER_COUNT] = {
    -6000, -5500, -5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500,
    -1000, -500, 0, 500, 1000, 1500, 2000, 2500, 3000, 3500,
    4000, 4500, 5000, 5500, 6000
};

/* Global variables */
int upfac = 4000;          // Upsampling factor
int downfac = 1023;        // Downsampling factor
int Nsamp_nu = 4000;       // Number of samples
double ca_code[MAX_CA_CODE_SIZE][SATELLITE_NO];  // Resampled CA code
int ca_code_size = 0;      // Size of resampled CA code
double prn_zp[NFFT][SATELLITE_NO];  // Zero-padded PRN
double complex phase[4000][DOPPLER_COUNT];  // Phase shifts for Doppler correction
double sat_magnitudes[SATELLITE_NO][NFFT][DOPPLER_COUNT];  // Correlation magnitudes

/* Function prototypes */
void fft(complex double *x, int n);
void RSP(const double fd[], double complex phase[][DOPPLER_COUNT], 
        complex double *fft_prn, complex double in[],
        double sat_magnitudes[SATELLITE_NO][NFFT][DOPPLER_COUNT], int current_sat);
        
/* Recursive FFT implementation */
void fft(complex double *x, int n) {
    // Base case: if n is 1 or less, return
    if(n <= 1) return;
    
    // Allocate memory for even and odd parts
    complex double *even = malloc(n/2 * sizeof(complex double));
    complex double *odd = malloc(n/2 * sizeof(complex double));
    
    // Split into even and odd indices
    for(int i = 0; i < n/2; i++) {
        even[i] = x[2*i];
        odd[i] = x[2*i + 1];
    }
    
    // Recursively compute FFT of even and odd parts
    fft(even, n/2);
    fft(odd, n/2);
    
    // Combine results
    for(int k = 0; k < n/2; k++) {
        complex double t = cexp(-2.0 * I * PI * k / n) * odd[k];
        x[k] = even[k] + t;
        x[k + n/2] = even[k] - t;
    }
    
    // Free allocated memory
    free(even);
    free(odd);
}

/* RSP (Replica Signal Processing) implementation */
void RSP(const double fd[], double complex phase[][DOPPLER_COUNT], 
        complex double *fft_prn, complex double *in,
        double sat_magnitudes[SATELLITE_NO][NFFT][DOPPLER_COUNT], int current_sat) {
    
    // Apply Doppler correction for each frequency bin
    for(int f = 0; f < DOPPLER_COUNT; f++) {
        complex double sgl_dp_zp[NFFT] = {0};  // Zero-padded signal with Doppler correction
        
        // Apply Doppler phase shift to input signal
        for(int i = 0; i < Nsamp_nu; i++) {
            sgl_dp_zp[i] = in[i] * phase[i][f];
        }

        // Perform FFT on Doppler-corrected signal
        fft(sgl_dp_zp, NFFT);

        // Multiply with PRN FFT (frequency domain correlation)
        complex double product[NFFT];
        for(int i = 0; i < NFFT; i++) {
            product[i] = sgl_dp_zp[i] * fft_prn[i];
        }

        // IFFT through conjugate FFT technique
        for(int i = 0; i < NFFT; i++) {
            product[i] = conj(product[i]);
        }
        fft(product, NFFT);
        
        // Store correlation magnitudes
        for(int i = 0; i < NFFT; i++) {
            double mag = cabs(product[i]) / NFFT;
            sat_magnitudes[current_sat][i][f] = mag;
        }
    }
}

int main() {
    extern complex double in[COMPLEX_DATA_SIZE];
    double time_axis[Nsamp_nu];
    complex double fft_prn[NFFT];

    // Initialize time axis
    double delt = 1.0 / FS_RX;
    for(int i = 0; i < Nsamp_nu; i++) {
        time_axis[i] = i * delt;
    }

    // Compute phase shifts for Doppler correction
    for(int i = 0; i < Nsamp_nu; i++) {
        for(int j = 0; j < DOPPLER_COUNT; j++) {
            phase[i][j] = cexp(-I * 2.0 * PI * time_axis[i] * fd[j]);
        }
    }

    // Resample PRN code (moved outside the satellite loop)
    int idx = 1;
    int count = 0;
    while(idx <= upfac * C_LEN && count < MAX_CA_CODE_SIZE) {
        int row = (int)ceil((double)idx / upfac) - 1;
        if(row >= C_LEN) break;

        for(int j = 0; j < SATELLITE_NO; j++) {
            ca_code[count][j] = navic_prn[row][j];
        }

        count++;
        idx += downfac;
    }
    ca_code_size = count;

    /* Initialize PRN matrix with zero-padding */
    for(int i = 0; i < Nsamp_nu; i++) {
        for(int j = 0; j < SATELLITE_NO; j++) {
            prn_zp[i][j] = ca_code[i][j];
        }
    }

    // Process each satellite
    for(int sat = 0; sat < SATELLITE_NO; sat++) {
        // Prepare PRN FFT for correlation
        complex double prn_signal[NFFT] = {0};
        for(int i = 0; i < NFFT; i++) {
            prn_signal[i] = prn_zp[i][sat];
        }

        // Compute FFT of PRN code
        fft(prn_signal, NFFT);
        for(int i = 0; i < NFFT; i++) {
            fft_prn[i] = conj(prn_signal[i]);
        }

        // Perform Replica Signal Processing
        RSP(fd, phase, fft_prn, in, sat_magnitudes, sat);
    }

    // Detection logic for satellites 
    double sat_max_values[SATELLITE_NO];
    int sat_max_nfft[SATELLITE_NO];
    int sat_max_doppler[SATELLITE_NO];
    bool detected[SATELLITE_NO] = {false};
    double ratios_left[SATELLITE_NO] = {0};
    double ratios_right[SATELLITE_NO] = {0};

    // Find maximum correlation peak for each satellite
    for (int s = 0; s < SATELLITE_NO; ++s) {
        double current_max = -DBL_MAX;
        for (int n = 0; n < NFFT; ++n) {
            for (int d = 0; d < DOPPLER_COUNT; ++d) {
                if (sat_magnitudes[s][n][d] > current_max) {
                    current_max = sat_magnitudes[s][n][d];
                    sat_max_nfft[s] = n;
                    sat_max_doppler[s] = d;
                }
            }
        }

        // Calculate range based on correlation peak position
        float range = (sat_max_nfft[s] * C) / FS_RX;

        printf("\nSatellite %d - Max magnitude: %f Range: %.3f  Doppler frequency: %.2f Hz\n",
               s + 1, current_max, range, fd[sat_max_doppler[s]]);

        // Extract correlation profile at maximum Doppler frequency
        double output[NFFT];
        for (int i = 0; i < NFFT; i++) {
            output[i] = sat_magnitudes[s][i][sat_max_doppler[s]];
        }

        // Find nearest local maxima for detection validation
        int left_idx = -1, right_idx = -1;
        for (int offset = 1; offset < NFFT; offset++) {
            int left = sat_max_nfft[s] - offset;
            int right = sat_max_nfft[s] + offset;
            
            // Find local maximum to the left
            if (left >= 0 && left_idx == -1) {
                bool is_max = true;
                if ((left > 0 && output[left] <= output[left - 1]) || 
                    (left < NFFT - 1 && output[left] <= output[left + 1])) {
                    is_max = false;
                }
                if (is_max) left_idx = left;
            }
            
            // Find local maximum to the right
            if (right < NFFT && right_idx == -1) {
                bool is_max = true;
                if ((right > 0 && output[right] <= output[right - 1]) || 
                    (right < NFFT - 1 && output[right] <= output[right + 1])) {
                    is_max = false;
                }
                if (is_max) right_idx = right;
            }
            
            // Break if both local maxima are found
            if (left_idx != -1 && right_idx != -1) break;
        }

        // Calculate peak-to-sidelobe ratios
        if (left_idx != -1) {
            ratios_left[s] = current_max / output[left_idx];
        }
        
        if (right_idx != -1) {
            ratios_right[s] = current_max / output[right_idx];
        }

        // Determine if satellite is detected based on peak-to-sidelobe ratios
        detected[s] = (ratios_left[s] > detection_threshold) || 
                       (ratios_right[s] > detection_threshold);
    }

    // Print detection summary
    printf("\nDetection Summary:\n==================\n");
    for (int s = 0; s < SATELLITE_NO; ++s) {
        printf("Satellite %d: %s\n", s + 1, detected[s] ? "DETECTED" : "Not detected");
    }

    return 0;
}