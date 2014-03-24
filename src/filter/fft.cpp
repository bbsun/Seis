#include "fftw3.h"
void fft(int sign, int n, fftw_complex * signal)
{
  fftw_plan p;
  p = fftw_plan_dft_1d(n,signal,signal,sign,FFTW_ESTIMATE);
  if(sign>0){
    for(int i=0;i<n;i++){
      signal[i][0] /=n;
      signal[i][1] /=n;
    }
  }
  fftw_execute(p);
  fftw_destroy_plan(p);
}
