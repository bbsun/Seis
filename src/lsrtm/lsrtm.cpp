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
#include "datastr.h"
#include "inversion.h"
using std::string;
using std::ifstream;
using std::cout;
using std::endl;
void test();
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
  
  int nt = 3500; // 14000
  int nx=150;
  int nz = 256;
  int sx = 75;
  int sz = 0;
  int pml = 30;
  float dt = 0.0006f*2.0f;
  float dx = 7.62f*3.5567;
  float dz = 7.62f*2;
  float ** v0=MyAlloc<float>::alc(nz,nx);
  float **  v=MyAlloc<float>::alc(nz,nx); 
  float **dv =MyAlloc<float>::alc(nz,nx);
  float    fr = 20.0f;
  int delay = 1.0f/fr/dt;
  int delaycal = 100;
  int NT = nt + delaycal;
  int numdelay = delay + delaycal;
  float * wav = rickerWavelet(dt,fr,numdelay,NT);
  OMP_CORE = 4;
  opern(v,VALUE,nz,nx,3000);
  for(int ix=0;  ix<nx;ix++)
  for(int iz=128;iz<nz;iz++)
    {
      v[ix][iz]=4000;
    }
  
  // initi
}
