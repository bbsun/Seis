#ifndef INIT_INVERSION_H
#define INIT_INVERSION_H
float ** setAbs(int nz, int nx, int pml);
float ** setPml(int nz, int nx,int pml);
float ** setVel(float **v, int nz, int nx, int pml);
float  ** setDv(float  **dv, int nz, int nx, int pml);
void SumSpray(float ** f, int n1, int n2);
/**
 *Shift a signal by m sample x_new=x[n-m]
 *@param in      input signal
 *@param out     output signal
 *@param m       shift in  sample 
 */
void shift(float * in, float * out, int nt, int m);
#endif
