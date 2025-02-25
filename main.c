#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RSP.h"

#define C_LEN 1023
#define SATELLITE_NO 7
#define FS_TX 1023000.0
#define FS_RX 2000000.0
#define T 0.001
#define C 3.0e8
#define NFFT 2048
#define NUM_FD 25
#define N 4

// Function prototypes for data loading
void load_navic_prn(double navic_prn[C_LEN][SATELLITE_NO]);
void load_input_signal(const char *filename, double in_signal[NFFT]);

int main() {
    double navic_prn[C_LEN][SATELLITE_NO];
    load_navic_prn(navic_prn);
   
    double fracRepFactor = FS_RX / FS_TX;
    int upfac = (int)(fracRepFactor + 0.5);
    int downfac = 1;
   
    // Allocate arrays for resampled CA codes and their FFTs.
    // Here the first index is the satellite number.
    Complex ca_code[SATELLITE_NO][NFFT];
    Complex fft_prn[SATELLITE_NO][NFFT];
   
    // Resample PRN codes and compute FFT for each satellite.
    for (int s = 0; s < SATELLITE_NO; s++) {
        for (int i = 0; i < NFFT; i++) {
            int idx = (i * downfac) / upfac;
            if (idx < C_LEN) {
                ca_code[s][i].real = navic_prn[idx][s];
                ca_code[s][i].imag = 0.0;
            } else {
                ca_code[s][i].real = 0.0;
                ca_code[s][i].imag = 0.0;
            }
        }
        // Compute FFT of the PRN code for satellite s.
        cfft(ca_code[s], NFFT);
        for (int i = 0; i < NFFT; i++) {
            fft_prn[s][i] = ca_code[s][i];
        }
    }
   
    // Define Doppler frequency axis (from -6000 Hz to 6000 Hz in steps of 500 Hz)
    double fd[NUM_FD];
    double delta_f = 500.0;
    for (int i = 0; i < NUM_FD; i++) {
        fd[i] = -6000.0 + i * delta_f;
    }
   
    // Process each test case
    for (int j = 1; j <= N; j++) {
        printf("Test Case: %d\n", j);
       
        // Load input signal (captured data) from CSV
        double temp_signal[NFFT];
        char filename[20];
        sprintf(filename, "in_data%d.csv", j);
        load_input_signal(filename, temp_signal);
       
        // Convert the real input signal to a Complex array (imaginary parts = 0)
        Complex in_signal[NFFT];
        for (int i = 0; i < NFFT; i++) {
            in_signal[i].real = temp_signal[i];
            in_signal[i].imag = 0.0;
        }
       
        // This matrix will hold the range-Doppler output for one satellite.
        double mag_out_rd[NFFT][NUM_FD] = {0};
        double peak_mag;
        int peak_range_idx, peak_doppler_idx;
       
        // Array to store the local peak ratio for each satellite.
        double ratios[SATELLITE_NO] = {0};
       
        // For each satellite, run RSP and compute local peak ratio.
        for (int prnid = 0; prnid < SATELLITE_NO; prnid++) {
            RSP(mag_out_rd, &peak_mag, &peak_range_idx, &peak_doppler_idx,
                fft_prn[prnid], in_signal, NFFT, NFFT, fd, NUM_FD, 1.0 / FS_RX);
           
            // Local peak ratio: compare global peak with its immediate neighbors.
            double left_neighbor = (peak_range_idx > 0) ? mag_out_rd[peak_range_idx - 1][peak_doppler_idx] : 0.0;
            double right_neighbor = (peak_range_idx < NFFT - 1) ? mag_out_rd[peak_range_idx + 1][peak_doppler_idx] : 0.0;
            double neighbor_max = (left_neighbor > right_neighbor) ? left_neighbor : right_neighbor;
            double ratio = (neighbor_max > 0) ? (peak_mag / neighbor_max) : 1000.0;
            ratios[prnid] = ratio;
            // For debugging (remove if not needed):
            // printf("Satellite %d: Ratio = %lf (Peak: %lf at (%d, %d))\n", prnid + 1, ratio, peak_mag, peak_range_idx, peak_doppler_idx);
        }
       
        // Select the satellite with the maximum ratio
        int selected_sat = -1;
        double max_ratio = 0.0;
        for (int prnid = 0; prnid < SATELLITE_NO; prnid++) {
            if (ratios[prnid] > max_ratio) {
                max_ratio = ratios[prnid];
                selected_sat = prnid;
            }
        }
       
        // Now print the results in MATLAB style:
        // Only the satellite with maximum ratio is reported as detected.
        for (int prnid = 0; prnid < SATELLITE_NO; prnid++) {
            if (prnid == selected_sat) {
                printf("Satellite %d: Detection Threshold Met\n", prnid + 1);
            } else {
                printf("Satellite %d: Detection Threshold NOT Met\n", prnid + 1);
            }
        }
    }
    return 0;
}