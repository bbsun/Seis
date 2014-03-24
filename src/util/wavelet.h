#ifndef WAVELET_H
#define WAVELET_H
/**
 * ricker wavelet generation 
 * the size of the wavelet is NT and NT = nt + nd1
 * here nd1 is used to shift the signal for calculation purpose
 * nd2 is the shift in the CSG
 * delay = nd1 + nd2 
 *@param dt         time sample interval
 *@param fr         main frequency 
 *@param delay      time dely of the signal in sample
 *@param NT         length of the wavelet
 *@return           ricker wavelet
 */
float * rickerWavelet(float dt, float fr, int delay, int NT);
/**
 * return the wavenumber for the signal in the frequency domain
 *@param dt         time sample interval
 *@param nt         length of the signal
 *@return           wavenumber for the signal
 */
float * wavenumber(float dt, int nt);
/**
 * compensate for the phase and amplitude change of the 2D propagation
 *@param data       signal
 *@param dt         time sample interval
 *@param nt         length of the signal
 */
void corWavelet2D(float *data, float dt,int nt);
#endif
