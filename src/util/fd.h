#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H
static float ** modeling(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, float * wav, float **vel );
static float ** forward (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, float * wav, float **vel,float ** dv );
static float ** adjoint (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, float * wav, float **vel,float ** rec);
#endif
