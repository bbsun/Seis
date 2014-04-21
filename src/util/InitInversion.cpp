#include "InitInversion.h"
#include "arraymath.h"
#include "memory.h"
#include "../filter/dct.h"
#include "../filter/hilbert.h"
#include "../filter/fft.h"
#include "wavelet.h"
#include <iostream>
using std::cout;
using std::endl;
float ** setAbs(int nz, int nx, int pml)
{
  int nzpml = nz + 2*pml;
  int nxpml = nx + 2*pml;
  float ** w=0; 
  w=MyAlloc<float>::alc(nzpml,nxpml);
  opern(w,VALUE,nzpml,nxpml,1.0f);
  float * weight=0; 
  weight = MyAlloc<float>::alc(pml);
  for(int i=0;i<pml;i++)
    weight[i] = 1.0f - (i+1)*1.0f/pml;
  for(int iz=0;iz<nzpml;iz++){
    for(int ix=0;ix<pml;ix++){
      w[ix][iz] = weight[pml-ix-1];
    }
    for(int ix=pml+nx;ix<nxpml;ix++){
      w[ix][iz] = weight[ix-pml-nx];
    }
  }
  for(int ix=pml;ix<pml+nx;ix++){
    for(int iz=0;iz<pml;iz++){
      w[ix][iz] = weight[pml-iz-1];
    }
    for(int iz=pml+nz;iz<nzpml;iz++){
      w[ix][iz] = weight[iz-pml-nz];
    }
  }
  MyAlloc<float>::free(weight);
  return w;
}
float ** setPml(int nz, int nx, int pml)
{
  int nzpml = nz+2*pml;
  int nxpml = nx+2*pml;
  float* damp = MyAlloc<float>::alc(pml);
  float** damp_pml =0;
  damp_pml = MyAlloc<float>::alc(nzpml,nxpml);
  opern(damp_pml,VALUE,nzpml,nxpml,1.0f);
  for (int ix=0; ix<pml; ++ix) {
    float ax = 0.015*(pml-1-ix)/(pml-1)*19;
    //ax= 0.0008*(pml-1-ix);
    damp[ix] = exp(-ax*ax);
  }
  /*for (int ix=0; ix<pml;++ix){
    damp[ix] = 1.0f-0.08f*(pml-ix-1)/(pml-1);
    }*/
  int az=pml,ax=pml,bz=pml,bx=pml+nx-1;
  int cz=pml+nz-1,cx=pml,dz0=pml+nz-1,dx0=pml+nx-1;
  for(int ip=1;ip<=pml;++ip)
    {
      // shang
      int az1 = az-ip,ax1 = ax-ip,bz1=bz-ip,bx2=bx+ip;
      int cz1 = cz+ip,cx1 = cx-ip,dz1=dz0+ip,dx1=dx0+ip;
      for(int ix =ax1;ix<=bx2;ix++)
	{
	  damp_pml[ix][az1] = damp[pml-ip];
	  damp_pml[ix][cz1] = damp[pml-ip];
	}
      for(int iz = az1;iz<=cz1;iz++)
	{
	  damp_pml[ax1][iz] = damp[pml-ip];
	  damp_pml[dx1][iz] = damp[pml-ip];
	}
      
    }
  /*  for (int ix=0; ix<pml; ++ix) {
      for (int iz=0; iz<nzpml; ++iz) {
      damp_pml[pml-ix-1    ][iz] = damp[pml-1-ix];
      damp_pml[nxpml-pml+ix][iz] = damp[pml-1-ix];
      }
      }
      
      for (int iz=0; iz<pml; ++iz) {
      for (int ix=pml; ix<nxpml-(pml); ++ix) {
      damp_pml[ix][pml-iz-1    ] = damp[pml-1-iz];
      damp_pml[ix][nzpml-pml+iz] = damp[pml-1-iz];
      }
      }*/
  MyAlloc<float>::free(damp);
  return damp_pml;
}
float ** setVel(float **v, int nz, int nx, int pml)
{
  int nzpml = nz+2*pml;
  int nxpml = nx+2*pml;
  float** vel = MyAlloc<float>::alc(nzpml,nxpml);
  for (int ix=0; ix<nx; ++ix)
    for (int iz=0; iz<nz; ++iz)
      vel[ix+pml][iz+pml] = v[ix][iz];
  
    for (int iz=0; iz<pml; ++iz) 
      for (int ix=0; ix<nxpml; ++ix) {
	vel[ix][iz] = vel[ix][pml];
	vel[ix][nzpml-pml+iz] = vel[ix][nzpml-pml-1];
      }
    
    for (int ix=0; ix<pml; ++ix) 
      for (int iz=0; iz<nzpml; ++iz) {
	vel[ix][iz] = vel[pml][iz];
	vel[nxpml-pml+ix][iz] = vel[nxpml-pml-1][iz];
      }
    return vel;
}


float ** setDv(float **dv, int nz, int nx, int pml)
{
  int nzpml = nz+2*pml;
  int nxpml = nx+2*pml;
  float** vel = MyAlloc<float>::alc(nzpml,nxpml);
  
  for (int ix=0; ix<nx; ++ix)
	for (int iz=0; iz<nz; ++iz)
	  vel[ix+pml][iz+pml] = dv[ix][iz];
  return vel;
}
void SumSpray(float ** f, int n1, int n2)
{
	for(int i1=0;i1<n1;i1++){
		for(int i2=1;i2<n2;i2++)
		{
			f[i2][i1] +=f[i2-1][i1];
		}
		for(int i2=0;i2<n2-1;i2++)
		{
			f[i2][i1] = f[n2-1][i1];
		}
	}

}
void shift(float * in, float * out, int nt, int shift)
{
  DCT dct(nt, FFTW_ESTIMATE);
  float * cf = MyAlloc<float>::alc(nt);
  float * cf1= MyAlloc<float>::alc(nt);
  float * cf2= MyAlloc<float>::alc(nt);
  float *  y1= MyAlloc<float>::alc(nt);
  float *  y2= MyAlloc<float>::alc(nt);
  dct.apply(in,cf,1);
  opern(cf1,cf,cf,COPY,nt);
  opern(cf2,cf,cf,COPY,nt);
  shift = -shift;
  float tmp = shift*3.1415296/nt;
  for(int i=0;i<nt;i++)
    {
      float cs = cos(tmp*i);
      float ss = sin(tmp*i);
      cf1[i] *= cs;
      cf2[i] *= ss;
    }
  dct.apply(cf1,y1,-1);
  dct.apply(cf2,y2,-1);
  opern(out,y1,y1,COPY,nt);
  Hilbert<float> hl(100);
  hl.apply(y2,y1,nt);
  opern(out,out,y1,SUB,nt);
  MyAlloc<float>::free(cf);
  MyAlloc<float>::free(cf1);
  MyAlloc<float>::free(cf2);
  MyAlloc<float>::free(y1);
  MyAlloc<float>::free(y2);
}
void shift(float **in0, float **out, int nt,int nx, int shift)
{
  float *  fr= wavenumber(1,nt);
  float * cf = MyAlloc<float>::alc(nt);
  float * cf1= MyAlloc<float>::alc(nt);
  float * cf2= MyAlloc<float>::alc(nt);
  float *  y1= MyAlloc<float>::alc(nt);
  float *  y2= MyAlloc<float>::alc(nt);
  float *  y3= MyAlloc<float>::alc(nt);
  fftwf_plan pf;
  fftwf_plan pi; 
  float * tmpf=MyAlloc<float>::alc(nt) ; 
  float * tmpi=MyAlloc<float>::alc(nt) ; 
  pf = fftwf_plan_r2r_1d(nt,tmpf,tmpf,FFTW_REDFT10,FFTW_ESTIMATE);
  pi = fftwf_plan_r2r_1d(nt,tmpi,tmpi,FFTW_REDFT01,FFTW_ESTIMATE);
  Hilbert<float> hl(150);			   
  shift = -shift;
  float tmp = shift*3.1415296/nt;
  float tmp1 = 1.0f/2/nt;
  float tmp2 = 2.0f/nt/2.0f/nt;
  for(int i = 0; i< nx ; i++){
    fftwf_execute_r2r(pf,in0[i],cf);
    for(int it = 0; it<nt;it++){
      float cs = cos(tmp*it);
      float ss = sin(tmp*it);
      float cfv = cf[it];
      cf1[it] = cfv*cs;
      cf2[it] = cfv*ss;
    }
    fftwf_execute_r2r(pi,cf1,y1);
    fftwf_execute_r2r(pi,cf2,y2);
    hl.apply(y2,y3,nt);
    for(int it = 0; it<nt;it++){
      out[i][it] = (y1[it] - y3[it]) *tmp1;
    }
  }
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(pi);
 
  MyAlloc<float>::free(fr);
  MyAlloc<float>::free(cf);
  MyAlloc<float>::free(cf1);
  MyAlloc<float>::free(cf2);
  MyAlloc<float>::free(y1);
  MyAlloc<float>::free(y2);
  MyAlloc<float>::free(y3);
  MyAlloc<float>::free(tmpf);
  MyAlloc<float>::free(tmpi);
}
void shiftFFT(float ** in0, float ** out, int nt,int nx, int shift)
{
  float *  fr= wavenumber(1,nt);
  float * cf = MyAlloc<float>::alc(nt);
  float * cf1= MyAlloc<float>::alc(nt);
  float * cf2= MyAlloc<float>::alc(nt);
  float *  y1= MyAlloc<float>::alc(nt);
  float *  y2= MyAlloc<float>::alc(nt);
  fftwf_complex * yc=(fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*nt);
  fftwf_plan p1; 
  fftwf_plan p2;
  fftwf_plan pf;
  fftwf_plan pi; 
  float * tmpf=MyAlloc<float>::alc(nt) ; 
  float * tmpi=MyAlloc<float>::alc(nt) ; 
  p1 = fftwf_plan_dft_1d(nt,yc,yc,-1,FFTW_ESTIMATE);
  p2 = fftwf_plan_dft_1d(nt,yc,yc,1,FFTW_ESTIMATE);
  pf = fftwf_plan_r2r_1d(nt,tmpf,tmpf,FFTW_REDFT10,FFTW_ESTIMATE);
  pi = fftwf_plan_r2r_1d(nt,tmpi,tmpi,FFTW_REDFT01,FFTW_ESTIMATE);
  shift = -shift;
  float tmp = shift*3.1415296/nt;
  float tmp1 = 1.0f/2/nt;
  float tmp2 = 2.0f/nt/2.0f/nt;
  for(int i = 0; i< nx ; i++){
    fftwf_execute_r2r(pf,in0[i],cf);
    for(int it = 0; it<nt;it++){
      float cs = cos(tmp*it);
      float ss = sin(tmp*it);
      float cfv = cf[it];
      cf1[it] = cfv*cs;
      cf2[it] = cfv*ss;
    }
    fftwf_execute_r2r(pi,cf1,y1);
    fftwf_execute_r2r(pi,cf2,y2);
    for(int it=0;it<nt;it++){ 
      yc[it][0] = y2[it];
      yc[it][1] = 0.0f;
    }
    fftwf_execute(p1);
    for(int it=1;it<nt;it++)
      if(fr[it]<0 && it!=nt/2){
	yc[it][0] = 0.0;
	yc[it][1] = 0.0;
      }
    yc[0][0] /=2;
    yc[0][1] /=2;
    yc[nt/2][0] /=2;
    yc[nt/2][1] /=2;
    fftwf_execute(p2);
    for(int it = 0; it<nt;it++){
      out[i][it] = 0.0f;
      out[i][it] = y1[it]*tmp1 - yc[it][1] *tmp2;
    }
  }
  fftwf_free(yc);
  fftwf_destroy_plan(p1);
  fftwf_destroy_plan(p2);
  fftwf_destroy_plan(pf);
  fftwf_destroy_plan(pi);
 
  MyAlloc<float>::free(fr);
  MyAlloc<float>::free(cf);
  MyAlloc<float>::free(cf1);
  MyAlloc<float>::free(cf2);
  MyAlloc<float>::free(y1);
  MyAlloc<float>::free(y2);
  MyAlloc<float>::free(tmpf);
  MyAlloc<float>::free(tmpi);
}
void shiftFFT(float * in, float * out, int nt, int shift)
{
  DCT dct(nt, FFTW_ESTIMATE);
  float * cf = MyAlloc<float>::alc(nt);
  float * cf1= MyAlloc<float>::alc(nt);
  float * cf2= MyAlloc<float>::alc(nt);
  float *  y1= MyAlloc<float>::alc(nt);
  float *  y2= MyAlloc<float>::alc(nt);
  float *  fr= wavenumber(1,nt);
  fftw_complex * yc=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*nt);
  dct.apply(in,cf,1);
  opern(cf1,cf,cf,COPY,nt);
  opern(cf2,cf,cf,COPY,nt);
  shift = -shift;
  float tmp = shift*3.1415296/nt;
  for(int i=0;i<nt;i++){
    float cs = cos(tmp*i);
    float ss = sin(tmp*i);
    cf1[i] *= cs;
    cf2[i] *= ss;
  }
  dct.apply(cf1,y1,-1);
  dct.apply(cf2,y2,-1);
  opern(out,y1,y1,COPY,nt);
  for(int i=0;i<nt;i++){ 
    yc[i][0] = y2[i];
    yc[i][1] = 0.0f;
  }
  fft(-1,nt,yc);
  for(int i=1;i<nt;i++)
    if(fr[i]<0 && i!=nt/2){
      yc[i][0] = 0.0;
      yc[i][1] = 0.0;
    }
  yc[0][0] /=2;
  yc[0][1] /=2;
  yc[nt/2][0] /=2;
  yc[nt/2][1] /=2;
  fft(1,nt,yc);
  for(int i=0;i<nt;i++)
    out[i] = out[i] - yc[i][1]*2.0f;
  fftw_free(yc);
  MyAlloc<float>::free(cf);
  MyAlloc<float>::free(cf1);
  MyAlloc<float>::free(cf2);
  MyAlloc<float>::free(y1);
  MyAlloc<float>::free(y2);
  MyAlloc<float>::free(fr);
}
void shiftSimple(float **in, float **out, int nt, int nx, int m)
{
	opern(out,VALUE,nt,nx,0.0f);
	for(int ix = 0; ix<nx;ix++){
		for(int it=0;it<nt;it++){
			int index = it - m;
			if(index>=0 && index < nt)
			out[ix][it] = in[ix][index];
		}
	}
}
