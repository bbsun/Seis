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
/**
 * Return the velocity model for finite difference calcuation
 * returned velocity model velsub will be used in calculation
 * the returned parameters velfxsub and nxsub is useful for
 * insert the record into the calculation model.
 *@param [out] velfxsub       first sample of the output in x-axis
 *@param [out] nzsub          depth sample number of the output
 *@param [out] nxsub          distance sample number of the output
 *@param [in ] vel            velocity model
 *@param [in ] nz             depth sample number of the input
 *@param [in ] nx             distance sample number of the input
 *@param [in ] velfx          first sample of the input in x-axis
 *@param [in ] is             shot number index, in range [0, ns-1]
 *@return                     sub velocty model for calculation
 */
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
/**
 * cast the grid for calculation to the grid for inversion.
 * swap the two model, the entrance of d will not be zero if its 
 * coordinate is within the range of m.
 *@param [in]  m                    grid for calculation
 *@param [out] d                    grid for inversion
 *@param [in]  mfx                  first sample of m in x-axis
 *@param [in]  dfx                  first sample of d in x-axis
 *@param [in]  nzm                  depth sample number of m
 *@param [in]  nxm                  distance sample number of m
 *@param [in]  nzd                  depth sample number of d
 *@param [in]  nxd                  distance sample number of d
 *@param [in]  add                  clear d or not
 */
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
/**
 * get the record for calcualtion or push back the calculated record to CSG
 * adjoint = false, CSG --> rec
 * adojint = true , rec --> CSG
 *@param [in] CSG                   commmon shot gather 
 *@param [in] is                    shot index
 *@param [in] ng                    number of gathers of CSG
 *@param [in] nt                    number of sample 
 *@param [in] rec                   rec for calculation 
 *@param [in] fx                    first sample of rec in x-axis
 *@param [in] nx                    number of gather for rec
 *@param [in] adj                   adjoint or not 
 */
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
      char tmpf[256];
      sprintf(tmpf,"%svelsub%d.dat",wdir,is);
      write(string(tmpf),nzsub,nxsub,velsub);
      cout<<"subvel for shot :" << is << " nzsub: "<<nzsub<<" nxsub : "<<nxsub<<" velfx : "<<velfxsub<<endl;
      char tmpf1[256];
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
      MyAlloc<float>::free(recsub);
    }
  MyAlloc<float>::free(vel);
  MyAlloc<float>::free(grid);
}
