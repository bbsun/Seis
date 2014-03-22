#ifndef INIT_INVERSION_H
#define INIT_INVERSION_H
float ** setAbs(int nz, int nx, int pml);
float ** setPml(int nz, int nx,int pml);
float ** setVel(float **v, int nz, int nx, int pml);
float  ** setDv(float  **dv, int nz, int nx, int pml);
void SumSpray(float ** f, int n1, int n2);
#endif
