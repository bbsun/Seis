#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H
float ** modeling(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel );
float ** forward (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** dv );
float ** adjoint (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** rec);
void processBoundary(float dt, float dx,float dz,int nx,int nz,int pml, float ** velpml,float ** u31,float ** u3,float ** u2,float ** u1);
void modeling2D_high(float **u0, float **u1, float **u2, float **vvzz,float **vvxx,int nzpml,int nxpml);
void modeling2D_low(float **u0, float **u1, float **u2, float **vvzz,float **vvxx,int nzpml,int nxpml);
#endif