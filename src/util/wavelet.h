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
 */
float * rickerWavelet(float dt, float fr, int delay, int NT);
#endif
