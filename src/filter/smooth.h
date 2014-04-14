/**
 *@file smooth.h 
 *@brief recursive gaussian smoothing.
 *@author Bingbing Sun
 *@version April-14-2014
 */
#ifndef SMOOTH_H_H
#define SMOOTH_H_H
/**
 * recursive gaussian smoothing.
 * sigma=1 equals to have 3 sample half window for gaussian smoothing, following how to 
 * use
 @code
 int vector_length = 10;
 int sigma = 3;
 // smoothing using zero slope boundary condtion
 Smooth::smooth1d(in,out,vector_length,sigma);
 // smoothing using zero value boundary condtion
 Smooth::smooth1dSym(in, out, vector_length,sigma);
 @endcode
 */
class Smooth
{
 public:
  /**
   * smooth 1d array using zero slope boundary condition.
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n size of the array
   *@param   sigma smoothing parameters
   */
  static void smooth1d(float *x, float *y, int n, float sigma);
  /**
   * smooth 2d array along its first(fast) direction using zero slope boundary condition
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d1(float **x,float **y,int n1,int n2, float sigma);
    /**
   * smooth 2d array along its second direction using zero slope boundary condition.
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d2(float **x,float **y,int n1,int n2, float sigma);
  /**
   * smooth 2d array along its both direcions using zero slope boundary condition.
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d12(float **x, float **y, int n1, int n2, float sigma1,float sigma2);
  /**
   * smooth 2d array along its both directions using zero slope boundary condition
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d21(float **x, float **y, int n1, int n2, float sigma1,float sigma2);
    /**
   * smooth 1d array using zero value boundary condition.
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n size of the array
   *@param   sigma smoothing parameters
   */
  static void smooth1dSym(float *x, float *y, int n, float sigma);
    /**
   * smooth 2d array along its first(fast) direction using zero value boundary condition
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d1Sym(float **x,float **y,int n1,int n2, float sigma);
      /**
   * smooth 2d array along its second direction using zero value boundary condition.
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d2Sym(float **x,float **y,int n1,int n2, float sigma);
  /**
   * smooth 2d array along its both direcions using zero value boundary condition.
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d12Sym(float **x, float **y, int n1, int n2, float sigma1,float sigma2);
   /**
   * smooth 2d array along its both directions using zero value boundary condition
   *@param[in]   x input array
   *@param[out]  y output(smoothed) array
   *@param       n1 size of the array in first direction
   *@param       n2 size of the array in second direction
   *@param   sigma smoothing parameters
   */
  static void smooth2d21Sym(float **x, float **y, int n1, int n2, float sigma1,float sigma2);
  
};
#endif
