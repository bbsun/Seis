#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H
float ** modeling(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel );
float ** forward (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** dv );
float ** adjoint (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** rec);
float ** forwordPlane(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sz,int gz, float **vel,float **sou,float **dv);
float ** adjointPlane(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sz,int gz, float **vel,float **sou,float **rec);
float ** rtm_true(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** rec);
float ** illum     (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz , float * wav, float **vel);
float ** illumPlane(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sz, int gz,  float **vel,float **sou);
void processBoundary(float dt, float dx,float dz,int nx,int nz,int pml, float ** velpml,float ** u31,float ** u3,float ** u2,float ** u1);
void SaveAtBoundary(float ***up,float ***down,float *** right, float *** left, float **u,int pml,int layer,int nz,int nx,bool adj,int it);
void modeling2D_high(float **u0, float **u1, float **u2, float **vvzz,float **vvxx,int nzpml,int nxpml);
void modeling2D_low(float **u0, float **u1, float **u2, float **vvzz,float **vvxx,int nzpml,int nxpml);
void FaQi(float **u, float **r, float **out, int nz, int nx, int pml);
#endif
