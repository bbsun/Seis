#include "wavelet.h"
#include "memory.h"
#include <math.h>
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
