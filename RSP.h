#ifndef RSP_H
#define RSP_H

#define C_LEN 1023
#define SATELLITE_NO 7
#define NFFT 2048
#define NUM_FD 25

// Complex number structure
typedef struct {
    double real;
    double imag;
} Complex;

// Function Prototypes
void RSP(double mag_out_rd[NFFT][NUM_FD], double* peak_mag, int* peak_range_idx, int* peak_doppler_idx,
         Complex fft_prn[NFFT], Complex in_signal[NFFT], int Nsamp, int nfft, double fd[NUM_FD], int num_fd, double delta);

void cfft(Complex *x, int n);
void cifft(Complex *x, int n);
void load_navic_prn(double navic_prn[C_LEN][SATELLITE_NO]);
void load_input_signal(const char *filename, double in_signal[NFFT]);

#endif