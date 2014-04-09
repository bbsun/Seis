/**
 *@file dct.h 
 *@brief discrete cosine transform for 1d array
 *@author Bingbing Sun
 *@version April-9-2014
 */
#ifndef DCT_H_H
#define DCT_H_H
#include "fftw3.h"
class DCT
{
 public:
  DCT(int n, unsigned flags);
  ~DCT();
  void apply(float *in , float *out, int sign);
 private:
  fftwf_plan pf;
  fftwf_plan pi;
  float * tmpf ;
  float * tmpi ;
  float cd;
  float sqrt2;
  int  n;
};
#endif
