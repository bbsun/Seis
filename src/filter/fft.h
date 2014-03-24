#ifndef FFT_H_H
#define FFT_H_H
#include "fftw3.h"
void fft(int sign, int n, fftw_complex * signal);
#endif
