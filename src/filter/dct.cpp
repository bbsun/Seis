#include "dct.h"
#include "fftw3.h"
#include "../util/memory.h"
#include "../util/debug.h"
#include "../util/arraymath.h"
DCT::DCT(int n, unsigned flags)
{
  sqrt2   = sqrt(2.0f);
  this->n = n;
  cd = 1.0f/(n*2)*sqrt2/2.0f;
  tmpf = MyAlloc<float>::alc(n);
  tmpi = MyAlloc<float>::alc(n);
  pf = fftwf_plan_r2r_1d(n,tmpf,tmpf,FFTW_REDFT10,flags);
  pi = fftwf_plan_r2r_1d(n,tmpi,tmpi,FFTW_REDFT01,flags);
  check(pf,"can not initialize pf in DCT(int n, unsigned flags) ");
  check(pi,"can not initialize pi in DCT(int n, unsigned flags) ");
}
DCT::~DCT()
{
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(pi);
  MyAlloc<float>::free(tmpf);
  MyAlloc<float>::free(tmpi);
}
void DCT::apply(float * in, float * out, int sign)
{
  static float n2 = n*2;
  if(sign>0){
    fftwf_execute_r2r(pf,in,out);
    for(int i=1;i<n;i++)
      out[i] *=sqrt2;
  }
  else{
    for( int i=1;i<n;i++)
      out[i] = in[i]*cd;
    out[0] =in[0]/n2;
    fftwf_execute_r2r(pi,out,out);
    
  }
}
