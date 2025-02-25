#include <stdio.h>
#include <stdlib.h>
#include "RSP.h"  // Include this to access C_LEN, SATELLITE_NO, and NFFT

void load_navic_prn(double navic_prn[C_LEN][SATELLITE_NO]) {
    FILE *file = fopen("navic_prn.csv", "r");
    if (!file) {
        perror("Error opening navic_prn.csv");
        exit(EXIT_FAILURE);
    }

    // Read PRN codes from CSV
    for (int i = 0; i < C_LEN; i++) {
        for (int j = 0; j < SATELLITE_NO; j++) {
            fscanf(file, "%lf,", &navic_prn[i][j]);
        }
    }

    fclose(file);
}

void load_input_signal(const char *filename, double in_signal[NFFT]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening input signal file");
        exit(EXIT_FAILURE);
    }

    // Read input signal data
    for (int i = 0; i < NFFT; i++) {
        fscanf(file, "%lf,", &in_signal[i]);
    }

    fclose(file);
}

