#include "mpi.h"
#include "inversion.h"
#include "../util/parser.h"
#include "../util/coord.h"
#include "../util/arraymath.h"
#include "../util/debug.h"
#include "../util/memory.h"
#include "../util/io.h"
#include <string>
#include <iostream>
using std::cout;
using std::endl;
using std::string;
#define verbose      false
Inversion::Inversion(int rank, int nprocs)
{
  this->rank   = rank;
  this->nprocs = nprocs;
}
Inversion::~Inversion()
{
  MyAlloc<int>::free(ng);
  MyAlloc<float>::free(sc);
  MyAlloc<float>::free(gc);
}
void Inversion::readParamFile(string file)
{
    ParserFromFile pff(file);
    // int parameters
    pff.getInt   ( param.ns        );
    pff.getInt   ( param.ngmax     );
    pff.getInt   ( param.nx        );
    pff.getInt   ( param.nz        );
    pff.getInt   ( param.nt        );
    pff.getInt   ( param.npml      );
    pff.getInt   ( param.delay     );
    pff.getInt   ( param.mask      );
    pff.getInt   ( param.maxiter   );
    pff.getInt   ( param.lpad      );
    pff.getInt   ( param.rpad      );
    // float parameters
    pff.getFloat ( param.dx        );
    pff.getFloat ( param.dz        );
    pff.getFloat ( param.dt        );
    pff.getFloat ( param.fr        );
    pff.getFloat ( param.velfx     );
    // string parameters
    pff.getString( param.precsg    );
    pff.getString( param.wdir      );
    pff.getString( param.v0file    );
    pff.getString( param.coordfile );
    // dependent 
    if(param.mask.val) {
      param.maskfile.opt=0; 
      pff.getString( param.maskfile);
    }
    
    if(rank==0 )
    param.print();
    // allocate data
    int ns    = param.ns.val;
    int ngmax = param.ngmax.val; 
    this->ng = MyAlloc<int>::alc(ns);
    this->sc = MyAlloc<float>::alc(2,ns);
    this->gc = MyAlloc<float>::alc(ngmax,2,ns);
}
void Inversion::getConfig()
{
  if(rank==0){
    string CoordFileName = param.coordfile.val;
    Coords * cord = new Coords();
    int ns = cord->countShots(CoordFileName);
    if(ns !=param.ns.val) {
      error(" error : ns ");
      cout<<endl;
      cout<<"ns = "<<param.ns.val <<" from parameter files "<<endl;
      cout<<"ns = "<<         ns  <<" from scanning "<<CoordFileName<<endl;
      exit(0);
    }
    int * ng=0;
    cord->countReceivers(CoordFileName,ns,ng);
    int ngmax = max(ng,ns);
    if(ngmax!=param.ngmax.val) {
      error(" error : ngmax ");
      cout<<endl;
      cout<<"ngmax = "<<param.ngmax.val <<" from parameter files "<<endl;
      cout<<"ngmax = "<<      ngmax     <<" from scanning "<<CoordFileName<<endl;
      exit(0);
    }
    float *** gc=0;
    float ** sc=0;
    cord->getSources(CoordFileName,ns,sc);
    cord->getReceivers(CoordFileName,ns,ng,gc);
    opern(this->ng,ng,ng,COPY,ns);
    opern(this->sc,sc,sc,COPY,2,ns);
    opern(this->gc,gc,gc,COPY,ngmax,2,ns);
    delete cord;
    MyAlloc<int>::free(ng);
    MyAlloc<float>::free(sc);
    MyAlloc<float>::free(gc);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  int ns    = param.ns.val;
  int ngmax = param.ngmax.val; 
  
  MPI_Bcast(&ng[0],       ns,        MPI_INT,   0, MPI_COMM_WORLD); 
  MPI_Bcast(&sc[0][0],    ns*2,      MPI_FLOAT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&gc[0][0][0], ns*2*ngmax,MPI_FLOAT, 0, MPI_COMM_WORLD); 
  if(rank==2 && verbose)
    {
      Coords::printShotInfo(1,10, ns,ngmax,this->ng,this->sc,this->gc);
      Coords::printShotInfo(ns-10,ns,ns,ngmax,this->ng,this->sc,this->gc);
    }
}

float ** Inversion::getVel(float &velfxsub, int &nzsub, int &nxsub, float ** vel,int nz, int nx, float velfx, int is)
{
  int lpad    = param.lpad.val;
  int rpad    = param.rpad.val;
  float dx    = param.dx.val;
  int ngt  = ng[is];           
  float sx = sc[is][0];
  float gx_min = min(gc[is][0],ngt);
  float gx_max = max(gc[is][0],ngt);
  float x_min  = min(sx, gx_min);
  float x_max  = max(sx, gx_max);
  float x_min_pad = x_min - lpad * dx;
  float x_max_pad = x_max + rpad * dx;
  float x_min_vel = velfx + dx;
  float x_max_vel = velfx + (nx-1) * dx -dx;
  velfxsub = x_min_pad;
  nzsub = nz;
  nxsub = (x_max_pad - x_min_pad)/dx + 1;
  float **velsub =  MyAlloc<float>::alc(nzsub,nxsub);
  for(int ix = 0; ix < nxsub ; ix++)
    {
      float xcor = velfxsub + ix * dx;
      if( (xcor >= x_min_vel) && (xcor <= x_max_vel ) )
	{
	  int xloc = (xcor-velfx)/dx; 
	  check(xloc>=0 && xloc<nx,"error: getVel in inversion.cpp xloc must be within range [0,nx)");
	    memcpy(velsub[ix], vel[xloc],sizeof(float)*nz); 
	}
      else if (xcor < x_min_vel)
	{
	  int xloc = (x_min_vel-velfx)/dx;
	  check(xloc>=0 && xloc<nx,"error: getVel in inversion.cpp xloc must be within range [0,nx)");
	  memcpy(velsub[ix], vel[xloc],sizeof(float)*nz); 
	}
      else if (xcor > x_max_vel)
	{
	  int xloc = (x_max_vel-velfx)/dx;
	  check(xloc>=0 && xloc<nx,"error: getVel in inversion.cpp xloc must be within range [0,nx)");
	  memcpy(velsub[ix], vel[xloc],sizeof(float)*nz); 
	}
      else
	{
	  error("error: getVel() in inversion.cpp") ;
	}
    }
  return velsub;
}

void Inversion::swapModel(float ** m, float ** d, float mfx, float dfx, int nzm,int nxm,int nzd,int nxd, bool add)
{
  check(nzm==nzd,"error: SwapModel in inversion.cpp nzm and nzd must be the same. ");
  if(!add) opern(d,VALUE,nzd,nxd,0.0f);
  float dx    = param.dx.val;
  float x_min = mfx ;
  float x_max = mfx + (nxm - 1) * dx ;
  for(int ix = 0; ix < nxd ; ix ++)
    {
      float xcor = dfx + ( ix - 1 ) * dx ;
      if( xcor >= x_min && xcor <= x_max)
	{
	  int xloc = (xcor - x_min) / dx;
	  check(xloc>=0 && xloc<nxm,"error: SwapModel in inversion.cpp xloc must be within range [0,nxm)");
	  for(int iz = 0; iz < nzd; iz ++)
	    d[ix][iz] += m[xloc][iz];
	}
    }
}

void Inversion::swapReord(float **CSG, int is, int ng, int nt, float ** rec, float fx, int nx, bool adj)
{
  adj? opern(CSG,VALUE,nt,ng,0.0f) : opern(rec,VALUE,nt,nx,0.0f);
  float dx = param.dx.val;
  float x_min = fx;
  float x_max = fx + ( nx - 1 ) * dx;
  for(int ig = 0; ig < ng  ; ig ++)
    {
      float gx = gc[is][0][ig];
      if( x_min<=gx && gx<=x_max)
	{
	  int xloc = (gx - x_min)/dx ;
	  check( xloc>=0 && xloc<nx , "error: getRecord in inversion.cpp xloc must be within range [0,nx)") ;
	  adj ? memcpy(CSG[ig],rec[xloc],sizeof(float)*nt) : memcpy(rec[xloc],CSG[ig],sizeof(float)*nt) ;
	}
    }
}
string Inversion::obtainCSGName  (int is)
{
  char name[256];
  sprintf(name,"%s%d.dat",param.precsg.val,is);
  return string(name);
}
string Inversion::obtainBornName (int is)
{
  char name[256];
  sprintf(name,"%sBORN%d.dat",param.wdir.val,is);
  return string(name);
}
string Inversion::obtianImageName(int is)
{
  char name[256];
  sprintf(name,"%sIMAGE%d.dat",param.wdir.val,is);
  return string(name);
}
string Inversion::obtainNameDat(string dir,string filename,int index)
{
  char name[256];
  sprintf(name,"%s%s%d.dat",dir,filename,index);
  return string(name);
}
string Inversion::obtainNameSu(string dir,string filename,int index)
{
  char name[256];
  sprintf(name,"%s%s%d.su",dir,filename,index);
  return string(name);
}
void Inversion::test()
{
  int nx = param.nx.val;
  int nz = param.nz.val;
  float velfx = param.velfx.val;
  float ** vel=MyAlloc<float>::alc(nz,nx);read(param.v0file.val,nz,nx,vel);
  float ** grid=MyAlloc<float>::alc(nz,nx);
  int ns = param.ns.val;
  int nt = param.nt.val;
  int ngmax = param.ngmax.val;
  float ** CSG=MyAlloc<float>::alc(nt,ngmax);
  // check the velocity model used for each shot
  string wdir = param.wdir.val;
  for(int is=400;is<ns;is++)
    {
      float velfxsub;
      int  nzsub;
      int  nxsub;
      float ** velsub = getVel(velfxsub, nzsub, nxsub, vel,nz, nx, velfx, is);
      float ** recsub = MyAlloc<float>::alc(nt,nxsub);
      string tmp = obtainNameSu(param.wdir.val,"velsub",is);
      writesu(tmp,nzsub,nxsub,velsub);
       cout<<"subvel for shot :" << is << " nzsub: "<<nzsub<<" nxsub : "<<nxsub<<" velfx : "<<velfxsub<<endl;
       /*char tmpf1[256];
      char tmpf2[256];
      sprintf(tmpf1,"%s%d.dat",param.precsg.val,is);
      sprintf(tmpf2,"%srecsub%d.dat",wdir,is);
      read(string(tmpf1),nt,ng[is],CSG);
      swapReord(CSG, is, ng[is], nt, recsub, velfxsub, nxsub, false);
      write(string(tmpf2),nt,nxsub,recsub);
      cout<<"recsub for shot :" << is << " nzsub: "<<nt<<" nxsub : "<<nxsub<<" velfx : "<<velfxsub<<endl;
      char tmpf3[256];
      sprintf(tmpf3,"%sgrid%d.dat",wdir,is);
      swapModel(velsub, grid, velfxsub, velfx, nzsub,nxsub, nz, nx, true);
      write(string(tmpf3),nz,nx,grid);
      char tmpf4[256];
      sprintf(tmpf4,"%sborn%d.dat",wdir,is);
      opern(CSG,RANDOM,nt,ng[is]);
      swapReord(CSG, is, ng[is], nt, recsub, velfxsub, nxsub, true);
      write(string(tmpf4),nt,ng[is],CSG);
      
      MyAlloc<float>::free(velsub);
      MyAlloc<float>::free(recsub);*/
    }
  MyAlloc<float>::free(vel);
  MyAlloc<float>::free(grid);
}
