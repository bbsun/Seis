/**
 *@file hilbert.h 
 *@brief hilbert transform for 1d array.
 *@author Bingbing Sun
 *@version April-18-2014
 */
#ifndef HILBERT_H_H
#define HILBERT_H_H
#include <string.h>
template <class T>
/**
 *Hilbert transform.
 *This is the Hilbert transform applied in the time domain
 */
class Hilbert
{
 public:
  /**
   *Constructor.
   *@param lhhalf     half of the filter's size 
   */
  Hilbert(int lhhalf)
    {
      this->_lhhalf = lhhalf;
      int lh = _lhhalf*2+1;
      this->_h   = new T[lh];
      designFilter();
    }
  /**
   *Destructor
   */
  ~Hilbert()
    {
      delete [] _h;
    }
  /**
   *Apply the Hilbert transform.
   *@param in      input array
   *@param out    output array
   *@param n      array size
   */
  void apply(T* in, T* out,int n)
  {
    memset(out,0,sizeof(T)*n);
    int lh = _lhhalf*2+1;
    int nbegin = _lhhalf+1;
    T * h = _h;
    for(int ib = 1;ib<=n;ib++)
      for(int im = 1;im<=lh;im++){
	int id = ib +im -1;
	if(id>=nbegin && id<=nbegin+n-1)
	out[id-nbegin] += in[ib-1]* h[im-1];
      }
  }
 private:
  /**
   *Design the filter in time domain. 
   */
  void designFilter()
  {
    int lhhalf = _lhhalf;
    T * h  = _h;
    h[lhhalf] = 0;
    for(int i=1;i<=lhhalf;i++)
      {
      T taper = 0.54 + 0.46 * cos(3.1415926*i/lhhalf);
      h[lhhalf+i] = taper *(i%2)*2.0/(3.1415926*i);
      h[lhhalf-i] = -h[lhhalf+i]; 
      }
  }
  T * _h;         /**< filter of Hilbert transform */
  int _lhhalf;    /**< half window size of Hilbert transform */
};

#endif
