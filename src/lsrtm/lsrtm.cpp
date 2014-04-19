#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "../util/memory.h"
#include "../util/arraymath.h"
#include "../util/coord.h"
#include "../util/parser.h"
#include "../util/wavelet.h"
#include "../util/global.h"
#include "../util/fd.h"
#include "../util/io.h"
#include "../util/InitInversion.h"
#include "datastr.h"
#include "inversion.h"
#include "../filter/dct.h"
#include "../filter/smooth.h"
using std::string;
using std::ifstream;
using std::cout;
using std::endl;
void test();
void sigsbee();
int main(int argc, char *argv[], char *envp[])
{
  int rank;
  int nprocs;
  //--initial MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  //--read param 
  
  Inversion inv(rank,nprocs);
  inv.readParamFile(argv[1]);
  inv.getConfig();
  inv.test();
  //--stop mp
  MPI_Finalize();
  
  return 0;
}
void test()
{
  if(true){
    cout<<"test of the time shift"<<endl;
    float * wav=rickerWavelet(0.0001,20,0,10000);
    writeSu("wav.su",10000,wav);
    float * shiftwav=MyAlloc<float>::alc(10000);
    shift(wav,shiftwav,10000,2000);
    writeSu("wavshift.su",10000,shiftwav);
    exit(0);
  }
  if(false){
  cout<<"test of the DCT and IDCT"<<endl;
  float * signal=MyAlloc<float>::alc(7);
  for(int i=0;i<7;i++)
    signal[i] = 5+i;
  dump("s",signal,7);
  DCT dct(7,FFTW_ESTIMATE);
  dct.apply(signal,signal,1);
  dump("sf",signal,7);
  dct.apply(signal,signal,-1);
  dump("sfi",signal,7);
  exit(0);
  float a=2.0; float b=3.0; float c = 4.0; float d= 5.0;
  float minf = min(a,b,c,d);
  float maxf = max(a,b,c,d);
  cout<<minf<<" min of "<< a << b << c << d<<endl;
  cout<<maxf<<" max of "<< a << b << c << d<<endl;
  if(true)
    {
      cout<<" test of modeling " << endl;
    }
  }
  int nt = 3500*2*2; // 14000
  int nx=150*2;
  int nz = 400;
  int sx = nx/2;
  int sz = 0;
  int pml = 30;
  float dt = 0.0003f;
  float dx = 7.62f*3.5567;
  float dz = 7.62f*2;
  int * igz = MyAlloc<int> :: alc(nx);
  float ** v0=MyAlloc<float>::alc(nz,nx);
  float **  v=MyAlloc<float>::alc(nz,nx); 
  float **dv =MyAlloc<float>::alc(nz,nx);
  float    fr = 20.0f;
  int delay = 100;
  int delaycal = 100;
  int NT = nt + delaycal;
  int numdelay = delay + delaycal;
  float * wav = rickerWavelet(dt,fr,numdelay,NT);
  //corWavelet2D(wav,dt,NT);
  float ** rec1 = MyAlloc<float>::alc(nt,nx);
  OMP_CORE = 3;
  opern(igz,VALUE,nx,6);
  opern(v,VALUE,nz,nx,4500.0f);
  opern(v0,VALUE,nz,nx,4000.0f);
  for(int ix=0;  ix<nx;ix++)
  for(int iz=0;iz<nz/2;iz++)
    {
      v[ix][iz]=4000.0f;
    }
  for(int ix=0;  ix<nx;ix++)
  for(int iz=nz/2;iz<nz;iz++)
    {
      dv[ix][iz]=1000.0f;
    }
  /*for(int ix=0;ix<nx;ix++)
    for(int iz=nz*2/3;iz<nz;iz++)
      {
	v[ix][iz]=4500.0f;
	}*/
  rec1=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v );
  float **rec4=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v0 );
  opern(rec1,rec1,rec4,SUB,nt,nx);
  write("rec.bin",nt,nx,rec1);
  read("rec.bin",nt,nx,rec1);
  //writeSu("rec.su",nt,nx,rec1);
  //writeSu("rec1.su",nt,rec1[sx]);
  Smooth::smooth2d1(v,v0,nz,nx,10);
  opern(dv,v,v0,SUB,nz,nx);
  writeSu("v0.su",nz,nx,v0);
  writeSu("v.su",nz,nx,v);
  writeSu("dv.su",nz,nx,dv);
  /*float ** born = forward ( dt, dx,  dz,  nt, delaycal, nx, nz, pml, sx, sz, igz, wav, v0, dv );
  write("born.bin",nt,nx,born);
  writeSu("born.su",nt,nx,born);*/
  float ** img  = adjoint ( dt, dx,  dz,  nt, delaycal, nx, nz, pml, sx, sz, igz, wav, v0, rec1);
  writeSu("img0.su",nz,nx,img);
  SumSpray(img,nz,nx);
  writeSu("img.su",nz,nx,img);
  exit(0);
  corWavelet2D(wav, dt,NT);
  float **rec2=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v );
  float **rec5=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v0 );
  opern(rec2,rec2,rec5,SUB,nt,nx);
  writeSu("rec2.su",nt,rec2[sx]);
  // MyAlloc<float>::free(rec1);
  //MyAlloc<float>::free(rec2);
  // initi
}
void sigsbee()
{
  int nt = 3500*2*2*2*2*2; // 14000
  int nx=  1000;
  int nz = 1201;
  int sx = 500;
  int sz = 10;
  int pml = 30;
  float dt = 0.00015f;
  float dx = 7.25;
  float dz = 7.25;
  int * igz = MyAlloc<int> :: alc(nx);
  float ** v0=MyAlloc<float>::alc(nz,nx);
  float **  v=MyAlloc<float>::alc(nz,nx); 
  float **dv =MyAlloc<float>::alc(nz,nx);
  float    fr = 10.0f;
  int delay = 100;
  int delaycal = 100;
  int NT = nt + delaycal;
  int numdelay = delay + delaycal;
  float * wav = rickerWavelet(dt,fr,numdelay,NT);
  //float **rec1=MyAlloc<float>::alc(nt,nx);
  read("/home/sunbb/sh/velx_1201_3201.dat",nz,nx,v);
  opern(v0,VALUE,nz,nx,v[0][0]);
  OMP_CORE = 10;
  opern(igz,VALUE,nx,6);
  float **rec1=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v );
  float **rec4=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v0 );
  
  //opern(rec1,rec1,rec4,SUB,nt,nx);
  write("rec_sig.dat",nt,nx,rec1);
  //read("rec_sig.dat",nt,nx,rec1);
  //writeSu("rec_sig.su",nt,nx,rec1);
  //writeSu("rec1_sig.su",nt,rec1[sx]);
  
}
