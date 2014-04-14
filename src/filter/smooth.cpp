#include <math.h>
#include "smooth.h"
#include "../util/memory.h"
void Smooth::smooth1d(float *x, float *y, int n, float sigma)
{
  float ss = sigma*sigma;
  float a  = (1.0f + ss -sqrt(1.0f + 2.0f*ss))/ss;
  float b = 1.0 -a;
  float sx = 1.0,sy=a;
  float yi=0.0;
  yi = sy*yi + sx*x[0];
  y[0] = yi;
  for(int i=1;i<n-1;i++){
    yi = a*yi + b*x[i];
    y[i] = yi;
  }
  sx /= (1.0+a);
  sy /= (1.0+a);
  yi = sy*yi + sx*x[n-1];
  y[n-1] = yi;
  for(int i=n-2;i>=0;i--){
    yi = a*yi + b*y[i];
    y[i] = yi;
  }
}

void Smooth::smooth2d1(float **x,float **y,int n1,int n2, float sigma)
{
  for(int i2=0;i2<n2;i2++)
    smooth1d(x[i2],y[i2],n1,sigma);
}

void Smooth::smooth2d2(float **x,float **y,int n1,int n2, float sigma)
{
  float * tx=MyAlloc<float>::alc(n2);
  float * ty=MyAlloc<float>::alc(n2);
  for(int i1=0;i1<n1;i1++)
    {
      for(int i2=0;i2<n2;i2++)
	tx[i2] = x[i2][i1];
      smooth1d(tx,ty,n2,sigma);
      for(int i2=0;i2<n2;i2++)
	y[i2][i1] = ty[i2];
    }
  MyAlloc<float>::free(tx);
  MyAlloc<float>::free(ty);
}
void Smooth::smooth2d12(float **x, float **y, int n1, int n2, float sigma1,float sigma2)
{
  smooth2d1(x,y,n1,n2,sigma1);
  smooth2d2(y,y,n1,n2,sigma2);
}
void Smooth::smooth2d21(float **x, float **y, int n1, int n2, float sigma1,float sigma2)
{
  smooth2d2(x,y,n1,n2,sigma2);
  smooth2d1(y,y,n1,n2,sigma1);
}
void Smooth::smooth1dSym(float *x, float *y, int n, float sigma)
{
  float ss = sigma*sigma;
  float a  = (1.0f + ss -sqrt(1.0f + 2.0f*ss))/ss;
  float b = 1.0f -a;
  float sx = b,sy=a;
  float yi=0.0;
  yi = sy*yi + sx*x[0];
  y[0] = yi;
  for(int i=1;i<n-1;i++){
    yi = a*yi + b*x[i];
    y[i] = yi;
  }
  sx /= (1.0+a);
  sy /= (1.0+a);
  yi = sy*yi + sx*x[n-1];
  y[n-1] = yi;
  for(int i=n-2;i>=0;i--){
    yi = a*yi + b*y[i];
    y[i] = yi;
  }
}
void Smooth::smooth2d1Sym(float **x,float **y,int n1,int n2, float sigma)
{
	for(int i2=0;i2<n2;i2++)
		smooth1dSym(x[i2],y[i2],n1,sigma);
}
void Smooth::smooth2d2Sym(float **x,float **y,int n1,int n2, float sigma)
{
  float * tx=MyAlloc<float>::alc(n2);
  float * ty=MyAlloc<float>::alc(n2);
  for(int i1=0;i1<n1;i1++)
    {
      for(int i2=0;i2<n2;i2++)
	tx[i2] = x[i2][i1];
      smooth1dSym(tx,ty,n2,sigma);
      for(int i2=0;i2<n2;i2++)
	y[i2][i1] = ty[i2];
    }
  MyAlloc<float>::free(tx);
  MyAlloc<float>::free(ty);
}
void Smooth::smooth2d12Sym(float **x, float **y, int n1, int n2, float sigma1,float sigma2)
{
  smooth2d1Sym(x,y,n1,n2,sigma1);
  smooth2d2Sym(y,y,n1,n2,sigma2);
}
void Smooth::smooth2d21Sym(float **x, float **y, int n1, int n2, float sigma1,float sigma2)
{
  smooth2d2Sym(x,y,n1,n2,sigma2);
  smooth2d1Sym(y,y,n1,n2,sigma1);
}
