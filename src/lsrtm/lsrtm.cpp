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
#include "../filter/taup.h"
#include "../filter/hilbert.h"
extern int OMP_CORE;
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
	//Inversion inv(rank,nprocs);
	//inv.readParamFile(argv[1]);
	//inv.getConfig();
	//inv.test();
  //--stop mp
   test();
	int bbb;
  MPI_Finalize();
  return 0;
}
void test()
{
  if(true){
  cout<<"test of the imaging condition for RTM"<<endl;
  float dt  = 0.0005;
  int   nt  = 14000;
  float dx  = 7.5f;
  float dz  = 7.5f;
  int   nx  = 1001;
  int   nz  = 600;
  int   sx  =  500;
  int   sz  =    0;
  int   pml =   30;
  int   nzpml = nz + pml*2;
  int   nxpml = nx + pml*2;
  float ** v=    MyAlloc<float>::alc(nz,nx);
  float **v0=    MyAlloc<float>::alc(nz,nx);
  float **rec=   MyAlloc<float>::alc(nt,nx);
  float **rech=  MyAlloc<float>::alc(nt,nx);
  float **vc=    MyAlloc<float>::alc(nz,nx);
  int * igz=MyAlloc<int>::alc(nx);
   cout<<"aaa"<<endl;
  for(int ix=0;ix<nx;ix++)
  igz[ix]=0;
  cout<<"bbb"<<endl;
  int delay = 0;
  int delaycal = 500;
  int NT = nt + delaycal;
  int numdelay = delay + delaycal;
  float fr = 25.0f;
  cout<<"bbb"<<endl;
  float *   wav = rickerWavelet(dt,fr,numdelay,NT);
  float *  wavh = MyAlloc<float>::alc(NT)         ;
   float *** imgX = MyAlloc<float>::alc(nzpml,nxpml,4);
  OMP_CORE = 16;
  cout<<"read begin"<<endl;
  read("/home/sunbb/sh/model/velx_600_1001.dat",nz,nx, v);
  read("/home/sunbb/sh/model/velsm_600_1001.dat",nz,nx,v0);
  cout<<"read finish"<<endl;
	for(int is = 0;is<nx;is+=4){
		sx = is;
		cout<<"shot "<< sx <<endl;
  rec         =      modeling( dt,  dx,  dz,  nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v);
  //write("rec.dat",nt,nx,rec);
  //read("rec0.dat", nt,nx,rec);
   opern(vc,VALUE,nz,nx,1500.0f);
   float **rec0         =      modeling( dt,  dx,  dz,  nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, vc);
   opern(rec,rec0,rec,SUB,nt,nx);
  //write("rec0.dat",nt,nx,rec);
  //read("rec0.dat",nt,nx,rec);
 
  Hilbert<float> hl(100);
  //hilbert transform of the sources 
  hl.apply(wav,wavh,nt);
  // hilbert of the record
  for(int ix=0;ix<nx;ix++)
  hl.apply(rec[ix],rech[ix],nt);
  float ** img = adjointUpDown( dt, dx, dz, nt, delaycal, nx, nz, pml, sx, sz, igz,  wav, wavh,v0, rec,  rech, imgX);
  writeSu("img_bing.su",nz,nx,img);
  writeSu("img00.su",nzpml,nxpml,imgX[0]);
  writeSu("imgzz.su",nzpml,nxpml,imgX[1]);
  writeSu("imgzt.su",nzpml,nxpml,imgX[2]);
  writeSu("imgtz.su",nzpml,nxpml,imgX[3]);
  MyAlloc<float>::free(rec);
  MyAlloc<float>::free(rec0);
  MyAlloc<float>::free(img);
}
  exit(0);
  }
  if(false){
    cout<<"test of the radon transform "<<endl;
    cout<<"large memory required "<<endl;
    int nt = 2000;
    int nx = 401;
    int np = 80;
	float dx = 10.0f/1000;
	float dt = 0.001;
	float p0 = -0.362f;
	float dp = 0.01f;
	float x0 = -200*dx;
    float ** rec = MyAlloc<float>::alc(nt,nx);
    float ** A   = MyAlloc<float>::alc(nt,np);
    read("/home/sunbb/sh/WEMVA/rec.bin",nt,nx,rec);
    for(int ix = 0;ix<200;ix++)
		opern(rec[ix],VALUE,nt,0.0f);
    Taup taup(dt,  nt, x0, dx,  nx, p0,  dp,  np );
    //Taup::apply(rec,A,false,dt,nt,x0,dx,nx,p0,dp,np);
    taup.apply(rec,A,false,false);
    write("/home/sunbb/sh/WEMVA/rec_mute.bin",nt,nx,rec);
    write("/home/sunbb/sh/WEMVA/a_2000_80.bin",nt,np,A);
    for(int ix=0;ix<nx;ix++)
		opern(rec[ix],VALUE,nt,0.0f);
    //Taup::apply(rec,A,true,dt,nt,x0,dx,nx,p0,dp,np);
    taup.apply(rec,A,true,false);
    write("/home/sunbb/sh/WEMVA/at_2000_80.bin",nt,nx,rec);
    exit(0);
  }
  if(false){
    cout<<"test of the time shift"<<endl;
    float * wav=rickerWavelet(0.0001,20,0,10000);
    writeSu("wav.su",10000,wav);
    float * shiftwav=MyAlloc<float>::alc(10000);
    //shift(wav,shiftwav,10000,2000);
    writeSu("wavshift.su",10000,shiftwav);
    //shiftFFT(wav,shiftwav,10000,2000);
    writeSu("wavshiftFFT.su",10000,shiftwav);
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
  int nt = 20000; // 14000
  int nx=  2807;
  int nz = 1500;
  int sx = nx/2;
  int sz = 300;
  int pml = 20;
  float dt = 0.0004f;
  float dx = 15;
  float dz = 10;
  int * igz = MyAlloc<int> :: alc(nx);
  float ** v0=MyAlloc<float>::alc(nz,nx);
  float **  v=MyAlloc<float>::alc(nz,nx); 
  float **dv =MyAlloc<float>::alc(nz,nx);
  float    fr = 15.0f;
  int delay =   0;
  int delaycal = 500;
  int NT = nt + delaycal;
  int numdelay = delay + delaycal;
  float * wav = rickerWavelet(dt,fr,numdelay,NT);
  //corWavelet2D(wav,dt,NT);
  float ** rec1 = MyAlloc<float>::alc(nt,nx);
  OMP_CORE = 30;
  read("/home/sunbb/sh/ch/veltest.dat",nz,nx,v);
  
  opern(igz,VALUE,nx,200);
  rec1=modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v );
  write("rec.bin",nt,nx,rec1);
  exit(0);

  
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
  rec1         =  modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v );
  float **rec4 =  modeling( dt,  dx,  dz, nt, delaycal, nx, nz, pml, sx,  sz,  igz, wav, v0 );
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
  writeSu( "img0.su",nz,nx,img  );
  SumSpray(img,nz,nx);
  writeSu( "img.su", nz,nx,img  );
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
