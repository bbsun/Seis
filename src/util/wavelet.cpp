#include "wavelet.h"
#include "memory.h"
#include <math.h>
#include <complex>
#include "../filter/fft.h"
using namespace std;
float * rickerWavelet(float dt, float fr, int numDelay, int nt)
{
  int   wavLenth  = nt;
  float * wavelet = MyAlloc<float>::alc(nt);
  float tmp=3.1415926*3.1415926*fr*fr;
  for(int i=0;i<wavLenth; i++){
    float t = (i-numDelay)*1.0f*dt;
    float t2=t*t;
    wavelet[i]=(1.0f-2.0f*tmp*t2)*exp(-tmp*t2);}
  return wavelet;
}
float * wavenumber(float dt, int n)
{
  float * f = MyAlloc<float>::alc(n);
  float dw = 2*3.1415926f/n/dt;
  int nth =n/2;
  for(int iii = 1; iii <=n ; iii++){
    int i1 = iii;
    f[iii-1]=dw*(i1-1);
    if(iii>nth+1){
      i1 = n-iii+2;
      f[iii-1]=-dw*(i1-1);}
  }
  return f;
}
void corWavelet2D(float *data, float dt,int nt)
{
  fftw_complex * y=0; y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nt);
  for(int i=0; i < nt;i++)
    {
      y[i][0]=data[i];    y[i][1]=0.0;
    }
  fft(-1,nt,y);
  float * f = wavenumber(dt,nt);
  // correct phase and amplitude correct
  for(int i = 0; i<nt ; i++)
    {
      float as = sqrt(abs(f[i]));
      complex <float> amp(as,0.0);
      float signnum = (f[i]>0.0)?1.0:((f[i]<0.0)?-1.0:0.0);
      float ps = 0.5*3.1415926/4 * signnum;
      float axx = cos(ps)*1.0f;
      float bxx = sin(ps)*1.0f;
      complex<float>  pha(axx,bxx);
      complex<float>  tmp(y[i][0],y[i][1]);
      complex<float>  b = tmp*amp*pha;
      y[i][0]= b.real();
      y[i][1]= b.imag();
    }
  fft(1,nt,y);
  for(int i = 0; i<nt; i++)
    data[i] = y[i][0];
  fftw_free(y);
}
