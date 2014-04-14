/**
 *@file dct.h 
 *@brief discrete cosine transform for 1d array
 *@author Bingbing Sun
 *@version April-9-2014
 */
#ifndef DCT_H_H
#define DCT_H_H
#include "fftw3.h"
/**
 * Discrete Cosine Transform. 
 * DCT is the same with matlab and following show how to 
 * use.
 @code
 DCT dct;
 dct(n,FFTW_ESTIMATE); // do not do tests
 dct(n,FFTW_MEASURE);  // do tests to measure the time
 dct.apply(in,out,1)   // for forward dct transform
 dct.apply(in,out,-1)  // for backward dct transform
 @endcode
 */
class DCT
{
 public:
  /**
   * Constructor
   *@param n size of the array for DCT
   *@param flags FFTW_ESTIMATE or FFTW_MEASURE 
   */
  DCT(int n, unsigned flags);
  /**
   * Destructor 
   */
  ~DCT();
  /**
   * Apply the DCT to array
   *@param in  input array 
   *@param out output array
   *@param sign sign>0 for forward DCT and sign<0 for inverse DCT
   */
  void apply(float *in , float *out, int sign);
 private:
  int n        ; /**< size of the array */
  float cd     ; /**< temp value */
  float sqrt2  ; /**< temp value */
  fftwf_plan pf; /**< fftw_plan for DCT  */
  fftwf_plan pi; /**< fftw_paln for IDCT */
  float * tmpf ; /**< temp array */
  float * tmpi ; /**< temp array */
};
#endif
