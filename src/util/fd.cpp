#include "global.h"
#include "fd.h"
#include "memory.h"
#include "debug.h"
#include "arraymath.h"
#include "InitInversion.h"
#include <iostream>
#include <omp.h>
#define LOOP  for(int ix=0;ix<nxpml;ix++) for(int iz=0;iz<nzpml;iz++)
using std::cout;
using std::endl;
extern int OMP_CORE;

float ** modeling(float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel )
{
  int nzpml = nz + 2*pml;
  int nxpml = nx + 2*pml;
  int szpml = sz + pml;
  int sxpml = sx + pml;
  float invdx2 = 1.0f/(dx*dx);
  float invdz2 = 1.0f/(dz*dz);
  int *  izpml  = MyAlloc<int> ::alc(nx);
  float ** w2d    = setAbs(nz,nx,pml);
  float ** velpml = setVel(vel,nz,nx,pml);
  float ** vt2    = MyAlloc<float>::alc(nzpml,nxpml); 
  float ** rec    = MyAlloc<float>::alc(nt,nx); 
  float ** u0     = MyAlloc<float>::alc(nzpml,nxpml);
  float ** u1     = MyAlloc<float>::alc(nzpml,nxpml);
  float ** u2     = MyAlloc<float>::alc(nzpml,nxpml);
  float ** u21    = MyAlloc<float>::alc(nzpml,nxpml);
  float ** vvzz   = MyAlloc<float>::alc(nzpml,nxpml);
  float ** vvxx   = MyAlloc<float>::alc(nzpml,nxpml);
  opern(vt2,velpml, velpml, COPY, nzpml, nxpml);
  opern(vt2,vt2   ,    vt2, SCAL, nzpml, nxpml, dt);
  opern(vt2,vt2   ,    vt2,  MUL, nzpml, nxpml);
  opern(vvzz,vt2,vt2,COPY,nzpml,nxpml);
  opern(vvxx,vt2,vt2,COPY,nzpml,nxpml);
  opern(vvzz,vvzz,vvzz,SCAL,nzpml,nxpml,invdz2);
  opern(vvxx,vvxx,vvxx,SCAL,nzpml,nxpml,invdx2);
  for(int ix =0 ; ix < nx ; ix++)
    izpml[ix] += pml;
  int NT = nt + ndel;
  for(int it=0;it<NT;it++)
    {
      if (it%1000==1)  cout<<"it="<<it<<endl ;
	
      modeling2D_high(u0,u1,u2,vvzz,vvxx,nzpml,nxpml);
      int itp = it - 1;
      if( itp >=0)
	u2[sxpml][szpml] += vt2[sxpml][szpml]*wav[itp];
      
      processBoundary(dt,dx,dz,nx,nz,pml,velpml,u21,u2,u1,u0);
      
      LOOP
	u2[ix][iz] = u2[ix][iz]*w2d[ix][iz] + u21[ix][iz]*(1.0f-w2d[ix][iz]); 
      
      itp = it - ndel;
      if(itp >=0)
	for(int ix=0;ix<nx;ix++){
	  int iz          = izpml[ix];
	  rec[ix][itp]    = u2[ix+pml][iz];
	}
      float ** tp = u0;
      u0 = u1;
      u1 = u2;
      u2 = tp;
    }	
  cout<<"Forwarding Modeling Finished"<<endl;
  MyAlloc<int>  ::free(izpml);
  MyAlloc<float>::free(u0);
  MyAlloc<float>::free(u1);
  MyAlloc<float>::free(u2);
  MyAlloc<float>::free(w2d);
  MyAlloc<float>::free(u21);
  MyAlloc<float>::free(vvzz);
  MyAlloc<float>::free(vvxx);
  MyAlloc<float>::free(velpml);
  MyAlloc<float>::free(vt2);
  return rec;
}
float ** forward (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** dv )
{

}
float ** adjoint (float dt, float dx, float dz, int nt, int ndel, int nx, int nz, int pml, int sx, int sz, int * gz, float * wav, float **vel,float ** rec)
{

}
void processBoundary(float dt, float dx,float dz,int nx,int nz,int pml, float ** velpml,float ** u31,float ** u3,float ** u2,float ** u1)
{
  int nxpml = nx + pml*2;
  int nzpml = nz + pml*2;
  float dt2 = dt*dt;
  float dx2 = (1.0f/dx)*(1.0f/dx)*dt2;
  float dz2 = (1.0f/dz)*(1.0f/dz)*dt2;
  for(int ix = 1;ix<=nxpml-2;ix++){
    u3[ix][0]    =velpml[ix][0]*velpml[ix][0]*dx2*(u2[ix+1][0]    + u2[ix-1][0]    - 2*u2[ix][0]) + 2*u2[ix][0] - u1[ix][0]; 
    u3[ix][nzpml-1] =velpml[ix][nzpml-1]*velpml[ix][nzpml-1]*dx2*(u2[ix+1][nzpml-1] + u2[ix-1][nzpml-1] - 2*u2[ix][nzpml-1])+ 2*u2[ix][nzpml-1] - u1[ix][nzpml-1];
  }
  for(int iz = 1;iz<=nzpml-2;iz++){
    u3[0][iz]    =velpml[0][iz]*velpml[0][iz]*dz2*(u2[0][iz+1]      + u2[0][iz-1]  - 2*u2[0][iz]) + 2*u2[0][iz] - u1[0][iz];
    u3[nxpml-1][iz] =velpml[nxpml-1][iz]*velpml[nxpml-1][iz]*dx2*(u2[nxpml-1][iz+1] + u2[nxpml-1][iz-1] - 2*u2[nxpml-1][iz])+ 2*u2[nxpml-1][iz] - u1[nxpml-1][iz];
  }
  u3[0][0]       =2*u2[0][0]      - u1[0][0];
  u3[0][nzpml-1]    =2*u2[0][nzpml-1]   - u1[0][nzpml-1];
  u3[nxpml-1][0]    =2*u2[nxpml-1][0]   - u1[nxpml-1][0];
  u3[nxpml-1][nzpml-1] =2*u2[nxpml-1][nzpml-1]- u1[nxpml-1][nzpml-1];
  for(int ix = 0;ix<nxpml;ix++){
    float yz1 = velpml[ix][1]*velpml[ix][1]*dz2*(u2[ix][2] + u2[ix][0] - 2* u2[ix][1]);
    float yz2 = velpml[ix][2]*velpml[ix][2]*dz2*(u2[ix][3] + u2[ix][1] - 2* u2[ix][2]);
    float yz0 = yz1*2 - yz2;
    u3[ix][0]     +=yz0;
    yz1 = velpml[ix][nzpml-2]*velpml[ix][nzpml-2]*dz2*(u2[ix][nzpml-3] + u2[ix][nzpml-1] - 2* u2[ix][nzpml-2]);
    yz2 = velpml[ix][nzpml-3]*velpml[ix][nzpml-3]*dz2*(u2[ix][nzpml-4] + u2[ix][nzpml-2] - 2* u2[ix][nzpml-3]);
    yz0 = yz1*2 - yz2;
    u3[ix][nzpml-1] +=yz0;
  }
  for(int iz = 0;iz<nzpml;iz++){
    float yz1 = velpml[1][iz]*velpml[1][iz]*dx2*(u2[2][iz] + u2[0][iz] - 2* u2[1][iz]);
    float yz2 = velpml[2][iz]*velpml[2][iz]*dx2*(u2[3][iz] + u2[1][iz] - 2* u2[2][iz]);
    float yz0 = yz1*2 - yz2;
    u3[0][iz]    +=yz0;
    yz1 = velpml[nxpml-2][iz]*velpml[nxpml-2][iz]*dx2*(u2[nxpml-3][iz] + u2[nxpml-1][iz] - 2* u2[nxpml-2][iz]);
    yz2 = velpml[nxpml-3][iz]*velpml[nxpml-3][iz]*dx2*(u2[nxpml-4][iz] + u2[nxpml-2][iz] - 2* u2[nxpml-3][iz]);
    yz0 = yz1*2 - yz2;
    u3[nxpml-1][iz] +=yz0;
  }
  for(int IK=1;IK<=pml;IK++){
    // up
    int K =  pml-IK;
    int b1  = pml+1-IK;
    int b2 =  pml+nx+IK-2;
    for (int ix =b1;ix<=b2;ix++){
      float a1 = 1.0f/dz/2.0f/dt;
      float a2 = 1.0f/2.0f/dt/dt/velpml[ix][K];
      float a3 = velpml[ix][K]/4.0f/dx/dx;
      float tmp = -a1*u3[ix][K+1] - a1*(u1[ix][K] - u1[ix][K+1]) 
	+a2*(u1[ix][K] - 2*u2[ix][K] + u3[ix][K+1] + u1[ix][K+1] - 2*u2[ix][K+1])
	-a3*(u3[ix+1][K+1] + u3[ix-1][K+1] - 2*u3[ix][K+1]+u1[ix+1][K] + u1[ix-1][K] - 2*u1[ix][K]);
      u31[ix][K] = -tmp/(a1 +a2);
    }
    // down
    K = nz + IK+pml-1;
    b1  = pml+1-IK;
    b2 =  pml+nx+IK-2;
    for (int ix =b1;ix<=b2;ix++){
      float a1 = 1.0f/dz/2.0f/dt;
      float a2 = 1.0f/2.0f/dt/dt/velpml[ix][K];
      float a3 = velpml[ix][K]/4.0f/dx/dx;
      float tmp = -a1*u3[ix][K-1] - a1*(u1[ix][K] - u1[ix][K-1]) 
	+a2*(u1[ix][K] - 2*u2[ix][K] + u3[ix][K-1] + u1[ix][K-1] - 2*u2[ix][K-1])
	-a3*(u3[ix+1][K-1] + u3[ix-1][K-1] - 2*u3[ix][K-1]+u1[ix+1][K] + u1[ix-1][K] - 2*u1[ix][K]);
      u31[ix][K] = -tmp/(a1 +a2);
    }
    // right 
    K = nx + IK+pml-1;
    b1 = pml+1-IK;
    b2 = pml+nz+IK-2;
    for (int iz = b1;iz<=b2;iz++){    
      float a1 = 1.0f/dx/2.0f/dt;
      float a2 = 1.0f/2.0f/dt/dt/velpml[K][iz];
      float a3 = velpml[K][iz]/4.0f/dz/dz;
      float tmp = -a1*u3[K-1][iz] - a1*(u1[K][iz] - u1[K-1][iz]) 
	+a2*(u1[K][iz] - 2*u2[K][iz] + u3[K-1][iz] + u1[K-1][iz] - 2*u2[K-1][iz]) 
	-a3*(u3[K-1][iz+1] + u3[K-1][iz-1] - 2*u3[K-1][iz]+u1[K][iz+1] + u1[K][iz-1] - 2*u1[K][iz]);
      u31[K][iz] = -tmp/(a1 +a2);
    }
    // left
    K =  pml+1 - IK-1;
    b1 = pml+1-IK;
    b2 = pml+nz+IK-2 ;
    for (int iz = b1;iz<=b2;iz++){    
      float a1 = 1.0f/dx/2.0f/dt;
      float a2 = 1.0f/2.0f/dt/dt/velpml[K][iz];
      float a3 = velpml[K][iz]/4.0f/dz/dz;
      float tmp = -a1*u3[K+1][iz] - a1*(u1[K][iz] - u1[K+1][iz]) 
	+a2*(u1[K][iz] - 2*u2[K][iz] + u3[K+1][iz] + u1[K+1][iz] - 2*u2[K+1][iz]) 
	-a3*(u3[K+1][iz+1] + u3[K+1][iz-1] - 2*u3[K+1][iz]+u1[K][iz+1] + u1[K][iz-1] - 2*u1[K][iz]);
      u31[K][iz] = -tmp/(a1 +a2);
    }
    // downright corner
    int zk = nz + pml+IK-1;
    int xk = nx + pml+IK-1;
    int zc0=zk;
    int xc0=xk;
    float ax1 = 1.0f/dz;
    float ax2 = 1.0f/dx;
    float ax3 = sqrt(2.0f)/velpml[xc0][zc0]/dt;
    float tmp = -ax1*u3[xc0][zc0-1] - ax2*u3[xc0-1][zc0]+ax3*(-u2[xc0][zc0]);
    u31[xc0][zc0] = -tmp/(ax1+ax2+ax3);
    // downleft corner
    zk = nz + pml+IK-1;
    xk = 1  + pml-IK-1;
    zc0 = zk;
    xc0 = xk;
    ax3 = sqrt(2.0f)/velpml[xc0][zc0]/dt;
    tmp = -ax1*u3[xc0][zc0-1] - ax2*u3[xc0+1][zc0]+ax3*(-u2[xc0][zc0]);
    u31[xc0][zc0] = -tmp/(ax1+ax2+ax3);
    // upright corner
    zk = 1 + pml-IK-1;
    xk = nx + pml+IK-1;
    zc0 = zk;
    xc0 = xk;
    ax3 = sqrt(2.0f)/velpml[xc0][zc0]/dt;
    tmp = -ax1*u3[xc0][zc0+1] - ax2*u3[xc0-1][zc0]+ax3*(-u2[xc0][zc0]);
    u31[xc0][zc0] = -tmp/(ax1+ax2+ax3);
    // upleft corner
    zk = 1  + pml-IK-1;
    xk = 1  + pml-IK-1;
    zc0 = zk;
    xc0 = xk;
    ax3 = sqrt(2)/velpml[xc0][zc0]/dt;
    tmp = -ax1*u3[xc0][zc0+1] - ax2*u3[xc0+1][zc0]+ax3*(-u2[xc0][zc0]);
    u31[xc0][zc0] = -tmp/(ax1+ax2+ax3);
  }
}
void modeling2D_low(float **u0, float **u1, float **u2, float **vvzz,float **vvxx,int nzpml,int nxpml)
{
#ifdef _OPENMP
  (void) omp_set_dynamic(0);
  if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
  (void) omp_set_num_threads(OMP_CORE);
#endif
#pragma omp parallel for
  for(int ix=1;ix<=nxpml-2;ix++)			
    for(int iz=1;iz<=nzpml-2;iz++)
      {
	float uzz  = u1[ix][iz-1] + u1[ix][iz+1] - 2*u1[ix][iz];
	float uxx  = u1[ix+1][iz] + u1[ix-1][iz] - 2*u1[ix][iz];
	u2[ix][iz] = vvzz[ix][iz]*uzz+vvxx[ix][iz]*uxx+2*u1[ix][iz]-u0[ix][iz];
      }
}
void modeling2D_high(float **u0, float **u1, float **u2, float **vvzz,float **vvxx,int nzpml,int nxpml)
{
  modeling2D_low(u0,u1,u2,vvzz,vvxx,nzpml,nxpml);
  register float c1 = -1.0f/560.0f;
  register float c2 =  8.0f/315.0f;
  register float c3 = -1.0f/5.0f  ;
  register float c4 =  8.0f/5.0f  ;
  register float c5 = -205.0f/72.0f;
  register float c6 =  8.0f/5.0f  ;
  register float c7 = -1.0f/5.0f  ;
  register float c8 =  8.0f/315.0f;
  register float c9 = -1.0f/560.0f;
#ifdef _OPENMP
  (void) omp_set_dynamic(0);
  if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
  (void) omp_set_num_threads(OMP_CORE);
#endif
#pragma omp parallel for
  for(int ix=4;ix<=nxpml-5;ix++)			
    for(int iz=4;iz<=nzpml-5;iz++)
      {
	float uzz  = u1[ix][iz-4]*c1 + u1[ix][iz-3]*c2 +  u1[ix][iz-2]*c3 + 
	  u1[ix][iz-1]*c4 + u1[ix][iz  ]*c5 +  u1[ix][iz+1]*c6 +
	  u1[ix][iz+2]*c7 + u1[ix][iz+3]*c8 +  u1[ix][iz+4]*c9 ;
	float uxx  = u1[ix-4][iz]*c1 + u1[ix-3][iz]*c2 +  u1[ix-2][iz]*c3 +
	  u1[ix-1][iz]*c4 + u1[ix  ][iz]*c5 +  u1[ix+1][iz]*c6 +
	  u1[ix+2][iz]*c7 + u1[ix+3][iz]*c8 +  u1[ix+4][iz]*c9 ;
	u2[ix][iz] = vvzz[ix][iz]*uzz+vvxx[ix][iz]*uxx+2*u1[ix][iz]-u0[ix][iz];
      }
}