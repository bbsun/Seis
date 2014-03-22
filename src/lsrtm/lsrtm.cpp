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
#include "datastr.h"
#include "inversion.h"
using std::string;
using std::ifstream;
using std::cout;
using std::endl;
void test();
int main(int argc, char *argv[], char *envp[])
{
  test();
  return;
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
  //--stop mpi
  MPI_Finalize();
  
  return 0;
}
void test()
{
  float a=2.0; float b=3.0; float c = 4.0; float d= 5.0;
  float minf = min(a,b,c,d);
  float maxf = max(a,b,c,d);
  cout<<minf<<" min of "<< a << b << c << d<<endl;
  cout<<maxf<<" max of "<< a << b << c << d<<endl;
  if(true)
    {
      cout<<" test of modeling " << endl;
    }
  
  int nt = 3500*2*2; // 14000
  int nx=150;
  int nz = 400;
  int sx = 75;
  int sz = 0;
  int pml = 30;
  float dt = 0.0003f;
  float dx = 7.62f*3.5567;
  float dz = 7.62f*2;
  float ** v0=MyAlloc<float>::alc(nz,nx);
  float **  v=MyAlloc<float>::alc(nz,nx); 
  float **dv =MyAlloc<float>::alc(nz,nx);
  float    fr = 20.0f;
  int delay = 100;
  int delaycal = 100;
  int NT = nt + delaycal;
  int numdelay = delay + delaycal;
  float * wav = rickerWavelet(dt,fr,numdelay,NT);
  OMP_CORE = 4;
  opern(v,VALUE,nz,nx,3000.0f);
  opern(v0,VALUE,nz,nx,3000.0f);
  for(int ix=0;  ix<nx;ix++)
  for(int iz=nz/2;iz<nz;iz++)
    {
      v[ix][iz]=4000.0f;
    }
  float **rec1=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  wav, v );
  float **rec4=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  wav, v0 );
  opern(rec1,rec1,rec4,SUB,nt,nx);
  writeSu("rec1.su",nt,rec1[sx]);
  corWavelet2D(wav, dt,NT);
  float **rec2=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  wav, v );
  float **rec5=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  wav, v0 );
  opern(rec2,rec2,rec5,SUB,nt,nx);
  writeSu("rec2.su",nt,rec2[sx]);
  MyAlloc<float>::free(rec1);
  MyAlloc<float>::free(rec2);
  // initi
}
