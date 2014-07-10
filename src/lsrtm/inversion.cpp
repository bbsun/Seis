#include "mpi.h"
#include "inversion.h"
#include "../util/parser.h"
#include "../util/coord.h"
#include "../util/arraymath.h"
#include "../util/debug.h"
#include "../util/memory.h"
#include "../util/io.h"
#include "../util/fd.h"
#include "../util/wavelet.h"
#include "../util/global.h"
#include "../util/InitInversion.h"
#include "../filter/hilbert.h"
#include "../filter/smooth.h"
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <omp.h>
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::to_string;
#define verbose      true
#define FALSE false
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
    pff.getInt   ( param.nthread   );
    pff.getInt   ( param.ns        ); 
    pff.getInt   ( param.ngmax     );
    pff.getInt   ( param.nx        );
    pff.getInt   ( param.nz        );
    pff.getInt   ( param.nt        );
    pff.getInt   ( param.npml      );
    pff.getInt   ( param.delay     );
    pff.getInt   ( param.delaycal  );
    pff.getInt   ( param.mask      );
    pff.getInt   ( param.maxiter   );
    pff.getInt   ( param.lpad      );
    pff.getInt   ( param.rpad      );
    pff.getInt   ( param.planeTag  );
    pff.getInt   ( param.sbegin    );
    pff.getInt   ( param.send      );
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
    pff.getString( param.vfile    );
    pff.getString( param.coordfile );
    // dependent 
    if(param.mask.val) {
      param.maskfile.opt=0; 
      pff.getString( param.maskfile);
    }
    if(param.planeTag.val) {
	  pff.getInt( param.np );
	  pff.getFloat( param.dp);
	  pff.getFloat( param.pmin);
	}
    if(rank==0 )
      param.print();
    // allocate data
    int ns    = param.ns.val;
    int ngmax = param.ngmax.val; 
    this->ng  = MyAlloc<int>::alc(ns);
    this->sc  = MyAlloc<float>::alc(2,ns);
    this->gc  = MyAlloc<float>::alc(ngmax,2,ns);
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
    int ngmin = min(ng,ns);
    cout<<"ngmax = " << ngmax << " from scanning "<<endl;
    cout<<"ngmin = " << ngmin << " from scanning "<<endl;
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
  /*if(rank==2 && verbose)
    {
      Coords::printShotInfo(1,10, ns,ngmax,this->ng,this->sc,this->gc);
      Coords::printShotInfo(ns-10,ns,ns,ngmax,this->ng,this->sc,this->gc);
      }*/
}
void Inversion::modeling_MPI( float ** v )
{
  if(rank == 0)
    cout<< " modeling "<<endl;
  
  int is=0;
  int ns   = param.ns.val;
  int idone = 0;
  int itotal= 0;
  MPI_Status  status;
  int myid = rank;
  if(myid==0)
    masterRun(ns);
  else
    {
      is = -1;
      MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);    // ask for a shot
      while ( is < ns)
	{
	  MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);  // receiver a new shot
	  if( is>=0 && is <ns ) {
	    idone ++; // begin
	    int nx   = param.nx.val;
	    int nz   = param.nz.val;
	    int nt   = param.nt.val;
	    int npml     = param.npml.val;
	    int delay    = param.delay.val;
	    int delaycal = param.delaycal.val;
	    int ngmax = param.ngmax.val;
	    float dt  = param.dt.val;
	    float dx  = param.dx.val;
	    float dz  = param.dz.val;
	    float fr  = param.fr.val;
	    float ** CSG = MyAlloc<float>::alc(nt,ngmax);
	    // check the velocity model used for each shot
	    string wdir = param.wdir.val;
	    int NT = nt + delaycal;
	    int numdelay = delay + delaycal;
	    float * wav = rickerWavelet(dt,fr,numdelay,NT);
	    float velfx  = param.velfx.val;
	    float velfxsub;
	    int  nzsub;
	    int  nxsub;
	    float ** velsub = getVel(velfxsub, nzsub, nxsub, v, nz, nx, velfx, is);
	    int * igz = getIgz(velfxsub,nxsub,is);
	    int sx = (int)((this->sc[is][0] - velfxsub)/dx+0.0001f);
	    int sz = (int)((this->sc[is][1] - 0.0     )/dz+0.0001f);
	    cout<<"shot "<< is << " source x: "<<sc[is][0] << " source z: "<<sc[is][1]<<endl;
	    float ** recsub  = modeling( dt, dx,  dz, nt, delaycal, nxsub, nzsub, npml, sx, sz,  igz, wav, velsub );
	    swapReord( CSG, is, this->ng[is], nt, recsub, velfxsub, nxsub, true );
	    removeDirect(CSG,v,is);
	    write( obtainCSGName(is), nt, this->ng[is], CSG );
	    MyAlloc<float>::free(CSG);
	    MyAlloc<float>::free(wav);
	    MyAlloc<float>::free(velsub);
	    MyAlloc<int>::free(igz);
	    MyAlloc<float>::free(recsub);
	    sleep(2); // end 
	    MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
	  }
	  else
	    MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	}
      MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}
void  Inversion::illum_MPI( float ** ss)
{
 
  string v0file = param.v0file.val;
  float ** v0 = MyAlloc<float>::alc(param.nz.val,param.nx.val);
  read(v0file,param.nz.val,param.nx.val,v0);
  if(rank==0)
    cout<<"illum"<<endl;
  int is=0;
  int idone = 0;
  int itotal= 0;
  int myid = rank;
  int ns  =  param.ns.val;
  MPI_Status  status;
  if(rank==0)
    masterRun(ns);
  else
    {
      is = -1;
      MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            //ask for a shot 
      while (is < ns) 
	{
	  MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);   //receive a new shot
	  if(is >=0 && is < ns){
	    idone++;
	    int nx   = param.nx.val;
	    int nz   = param.nz.val;
	    int nt   = param.nt.val;
	    int npml     = param.npml.val;
	    int delay    = param.delay.val;
	    int delaycal = param.delaycal.val;
	    int ngmax = param.ngmax.val;
	    float dt  = param.dt.val;
	    float dx  = param.dx.val;
	    float dz  = param.dz.val;
	    float fr  = param.fr.val;
	    // check the velocity model used for each shot
	    string wdir = param.wdir.val;
	    int NT = nt + delaycal;
	    int numdelay = delay + delaycal;
	    float * wav = rickerWavelet(dt,fr,numdelay,NT);
	    float velfx  = param.velfx.val;
	    float velfxsub;
	    int  nzsub;
	    int  nxsub;
	    float ** velsub = getVel(velfxsub, nzsub, nxsub, v0, nz, nx, velfx, is);
	    int sx = (int)((this->sc[is][0] - velfxsub)/dx+0.0001f);
	    int sz = (int)((this->sc[is][1] - 0.0     )/dz+0.0001f);
	    cout<<"shot "<< is << " source x: "<<sc[is][0] << " source z: "<<sc[is][1]<<endl;
	    float ** sssub = illum( dt,  dx , dz, nt, delaycal, nxsub, nzsub, npml,  sx, sz, wav, velsub);
	    float ** sstmp = MyAlloc<float>::alc(nz,nx);
	    swapModel(sssub, sstmp, velfxsub, velfx, nzsub,nxsub, nz,nx, false);
	    string ssfile = obtainNameDat(param.wdir.val,"ILLUM",is);
	    write(ssfile,nz,nx,sstmp);
	    MyAlloc<float>::free(wav);
	    MyAlloc<float>::free(velsub);
	    MyAlloc<float>::free(sssub);
	    MyAlloc<float>::free(sstmp);
	    sleep(2); 
	    MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
	  } else {
	    MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	    // fprintf(stderr,"Slave=%d Finshed, waiting for exit \n", myid);
	  }
	}
      MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  MyAlloc<float>::free(v0);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){	
    int nx   = param.nx.val;
    int nz   = param.nz.val;
    opern(ss,VALUE, nz , nx,0.0f);
    float ** sstmp = MyAlloc<float>::alc(nz,nx);
    for( int is=0; is<ns; is++){
      string ssfile = obtainNameDat(param.wdir.val,"ILLUM",is);
      read(ssfile,nz,nx,sstmp);
      opern(ss,sstmp,ss,ADD,nz,nx);
    }
    MyAlloc<float>::free(sstmp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&ss[0][0], param.nz.val*param.nx.val,MPI_FLOAT,0,MPI_COMM_WORLD); 
}
void  Inversion::illumPlane_MPI( float ** ss)
{
 
  string v0file = param.v0file.val;
  float ** v0 = MyAlloc<float>::alc(param.nz.val,param.nx.val);
  read(v0file,param.nz.val,param.nx.val,v0);
  if(rank==0)
    cout<<"illumPlane"<<endl;
  int is=0;
  int idone = 0;
  int itotal= 0;
  int myid = rank;
  int ns  =  param.np.val;
  MPI_Status  status;
  if(rank==0)
    masterRun(ns);
  else
    {
      is = -1;
      MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            //ask for a shot 
      while (is < ns) 
	{
	  MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);   //receive a new shot
	  if(is >=0 && is < ns){
	    idone++;
	    int nx   = param.nx.val;
	    int nz   = param.nz.val;
		int nt   = param.nt.val;
		int npml     = param.npml.val;
		int delaycal = param.delaycal.val;
		float dt  = param.dt.val;
		float dx  = param.dx.val;
		float dz  = param.dz.val;
		float fr  = param.fr.val;
		string wdir = param.wdir.val;
		int NT = nt + delaycal;
		float ** sou = MyAlloc<float>::alc(NT,nx);
		string soufile;
		soufile  = obtainNameDat(param.wdir.val,"PLANE_S_NEW",is);
		read(soufile,NT,nx,sou);
		int sz = (sc[0][1]-0.0f)/dz + 0.0001f;
		int gz = (gc[0][1][0]-0.0f)/dz + 0.0001f;
		cout<< " p  "<< param.pmin.val + is * param.dp.val <<" sz "<< sc[0][1] << " gz "<< gc[0][1][0] <<endl;
	    float ** illumX = illumPlane( dt, dx, dz,  nt,  delaycal,  nx,  nz,  npml, sz, gz, v0, sou);
	    string ssfile = obtainNameDat(param.wdir.val,"ILLUM_PLANE",is);
	    write(ssfile,nz,nx,illumX);
	    MyAlloc<float>::free(sou);
	    MyAlloc<float>::free(illumX);
	    sleep(2); 
	    MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
	  } else {
	    MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	    // fprintf(stderr,"Slave=%d Finshed, waiting for exit \n", myid);
	  }
	}
      MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  MyAlloc<float>::free(v0);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){	
    int nx   = param.nx.val;
    int nz   = param.nz.val;
    opern(ss,VALUE, nz , nx,0.0f);
    float ** sstmp = MyAlloc<float>::alc(nz,nx);
    for( int is=0; is<ns; is++){
      string ssfile = obtainNameDat(param.wdir.val,"ILLUM_PLANE",is);
      read(ssfile,nz,nx,sstmp);
      opern(ss,sstmp,ss,ADD,nz,nx);
    }
    MyAlloc<float>::free(sstmp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&ss[0][0], param.nz.val*param.nx.val,MPI_FLOAT,0,MPI_COMM_WORLD); 
} 
void  Inversion::adjoint_MPI( float ** img, int migTag)
{
 
  string v0file = param.v0file.val;
  float ** v0 = MyAlloc<float>::alc(param.nz.val,param.nx.val);
  read(v0file,param.nz.val,param.nx.val,v0);
  if(rank==0)
    cout<<"adjoint"<<endl;
  
  int is=0;
  int idone = 0;
  int itotal= 0;
  int myid = rank;
  int ns  =  param.ns.val;
  int sbegin = param.sbegin.val;
  int send   = param.send.val;
  MPI_Status  status;
  if(rank==0)
    masterRun(ns);
  else
    {
      is = -1;
      MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            //ask for a shot 
      while (is < ns) 
	{
	  MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);   //receive a new shot
	  if(is >=0 && is < ns){
	    idone++;
	    if(is>=sbegin && is<=send)
	      {
		int nx   = param.nx.val;
		int nz   = param.nz.val;
		int nt   = param.nt.val;
		int npml     = param.npml.val;
		int delay    = param.delay.val;
		int delaycal = param.delaycal.val;
		int ngmax = param.ngmax.val;
		float dt  = param.dt.val;
		float dx  = param.dx.val;
		float dz  = param.dz.val;
		float fr  = param.fr.val;
		// check the velocity model used for each shot
		string wdir = param.wdir.val;
		int NT = nt + delaycal;
		int numdelay = delay + delaycal;
		float * wav = rickerWavelet(dt,fr,numdelay,NT);
		float ** CSG = MyAlloc<float>::alc(nt,ngmax);
		float ** mask= MyAlloc<float>::alc(nz,nx);
		float velfx  = param.velfx.val;
		float velfxsub;
		int  nzsub;
		int  nxsub;
		float ** velsub = getVel(velfxsub, nzsub, nxsub, v0, nz, nx, velfx, is);
		float ** recsub = MyAlloc<float>::alc(nt,nxsub);
		string recfile;
		if(migTag==RTM_IMG)
		  recfile = obtainCSGName(is);
		else
		  recfile = obtainBornName(is);
		read(recfile,nt,ng[is],CSG);
		removeDirect(CSG,v0,is);
		swapReord( CSG, is, this->ng[is], nt, recsub, velfxsub, nxsub, false );
		int * igz = getIgz(velfxsub,nxsub,is);
		int sx = (int)((this->sc[is][0] - velfxsub)/dx+0.0001f);
		int sz = (int)((this->sc[is][1] - 0.0     )/dz+0.0001f);
		cout<<"shot "<< is << " source x: "<<sc[is][0] << " source z: "<<sc[is][1]<<endl;
		float ** imgsub = adjoint ( dt,  dx , dz, nt, delaycal, nxsub, nzsub, npml,  sx, sz, igz, wav, velsub, recsub);
		float ** imgtmp = MyAlloc<float>::alc(nz,nx);
		swapModel(imgsub, imgtmp, velfxsub, velfx, nzsub,nxsub, nz,nx, false);
		string imgfile = obtainNameDat(param.wdir.val,"CSG_IMAGE",is);
		if(param.mask.val){
		  read(param.maskfile.val,nz,nx,mask);
		  opern(imgtmp,imgtmp,mask,MUL,nz,nx);
		}
		write(imgfile,nz,nx,imgtmp);
		MyAlloc<int>::free(igz);
		MyAlloc<float>::free(CSG);
		MyAlloc<float>::free(wav);
		MyAlloc<float>::free(velsub);
		MyAlloc<float>::free(recsub);
		MyAlloc<float>::free(imgsub);
		MyAlloc<float>::free(imgtmp);
		MyAlloc<float>::free(mask);
	      }
	    sleep(2); 
	    MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
	  } else {
	    MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	    // fprintf(stderr,"Slave=%d Finshed, waiting for exit \n", myid);
	  }
	}
      MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  MyAlloc<float>::free(v0);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){	
    int nx   = param.nx.val;
    int nz   = param.nz.val;
    opern(img,VALUE, nz , nx,0.0f);
    float ** imgtmp = MyAlloc<float>::alc(nz,nx);
    for( int is=0; is<ns; is++){
      string imgfile = obtainNameDat(param.wdir.val,"CSG_IMAGE",is);
      read(imgfile,nz,nx,imgtmp);
      opern(img,imgtmp,img,ADD,nz,nx);
    }
    MyAlloc<float>::free(imgtmp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&img[0][0], param.nz.val*param.nx.val,MPI_FLOAT,0,MPI_COMM_WORLD); 
}
void  Inversion::forward_MPI( float ** dv )
{
	
  string v0file = param.v0file.val;
  float ** v0 = MyAlloc<float>::alc(param.nz.val,param.nx.val);
  read(v0file,param.nz.val,param.nx.val,v0);
  if(rank==0)
    cout<<"forward"<<endl;
  
  int is=0;
  int idone = 0;
  int itotal= 0;
  int myid = rank;
  int ns  =  param.ns.val;
  MPI_Status  status;
  if(rank==0)
    masterRun(ns);
  else
    {
      is = -1;
      MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            //ask for a shot 
      while (is < ns) 
	{
	  MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);   //receive a new shot
	  if(is >=0 && is < ns){
	    idone++;
	    int nx   = param.nx.val;
	    int nz   = param.nz.val;
	    int nt   = param.nt.val;
	    int npml     = param.npml.val;
	    int delay    = param.delay.val;
	    int delaycal = param.delaycal.val;
	    int ngmax = param.ngmax.val;
	    float dt  = param.dt.val;
	    float dx  = param.dx.val;
	    float dz  = param.dz.val;
	    float fr  = param.fr.val;
	    // check the velocity model used for each shot
	    string wdir = param.wdir.val;
	    int NT = nt + delaycal;
	    int numdelay = delay + delaycal;
	    float * wav = rickerWavelet(dt,fr,numdelay,NT);
	    float ** CSG = MyAlloc<float>::alc(nt,ngmax);
	    float velfx  = param.velfx.val;
	    float velfxsub;
	    int  nzsub;
	    int  nxsub;
	    float ** velsub = getVel(velfxsub, nzsub, nxsub, v0, nz, nx, velfx, is);
	    float ** dvsub  = getVel(velfxsub, nzsub, nxsub, dv, nz, nx, velfx, is);
	    int * igz = getIgz(velfxsub,nxsub,is);
	    int sx = (int)((this->sc[is][0] - velfxsub)/dx+0.0001f);
	    int sz = (int)((this->sc[is][1] - 0.0     )/dz+0.0001f);
	    if(verbose)	
	      cout<<"shot "<< is << " source x: "<<sc[is][0] << " source z: "<<sc[is][1]<<endl;
	    float ** recsub = forward(dt, dx, dz, nt, delaycal,nxsub,nzsub,npml, sx, sz, igz , wav, velsub, dvsub);
	    string bornfile;
	    bornfile = obtainBornName(is);
	    swapReord( CSG, is, this->ng[is], nt, recsub, velfxsub, nxsub, true );
	    removeDirect(CSG,v0,is);
	    write(bornfile,nt,ng[is],CSG);
	    MyAlloc<int>::free(igz);
	    MyAlloc<float>::free(CSG);
	    MyAlloc<float>::free(wav);
	    MyAlloc<float>::free(dvsub);
	    MyAlloc<float>::free(velsub);
	    MyAlloc<float>::free(recsub);
	    sleep(2); 
	    MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
	  } else {
	    MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	    // fprintf(stderr,"Slave=%d Finshed, waiting for exit \n", myid);
	  }
	}
      MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
  MyAlloc<float>::free(v0);
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

int*  Inversion::getIgz ( float velfx, int nx ,int is)
{
  int * igz = MyAlloc<int>::alc(nx);
  
  float dx = param.dx.val;
  float dz = param.dz.val;
  float x_min = velfx -0.001*dx;
  float x_max = velfx + (nx -1 ) * dx +0.001*dx;
  int NZ = param.nz.val;
  for( int ix = 0; ix < nx ; ix++ )
      igz[ix] = 0;
  
  for( int ig = 0; ig < ng[is]; ig++){
    float x = gc[is][0][ig];
    float z = gc[is][1][ig];
    if( x>= x_min && x<= x_max){
      int xloc = (x - x_min )/dx + 0.1f;
      int zloc = (z - 0.0f) /dz+0.1f;
      //cout<<"is "<< is <<" ig "<<ig<<" xloc "<<xloc<<" zloc "<< zloc << endl;
      check(xloc>=0 && xloc<nx,"error: getIgz in inversion.cpp xloc must be within range [0,nx)");
      check(zloc>=0 && zloc<NZ,"error: getIgz in inversion.cpp zloc musg be within range [0,NZ)");
      igz[xloc] = zloc;
    }
  }
  return igz;
}
void Inversion::swapReord(float **CSG, int is, int ng, int nt, float ** rec, float fx, int nx, bool adj)
{
  adj? opern(CSG,VALUE,nt,ng,0.0f) : opern(rec,VALUE,nt,nx,0.0f);
  float dx = param.dx.val;
  float x_min = fx - 0.001*dx;
  float x_max = fx + ( nx - 1 ) * dx +0.001*dx;
  for(int ig = 0; ig < ng  ; ig ++)
    {
      float gx = gc[is][0][ig];
      if( x_min<=gx && gx<=x_max)
	{
	  int xloc = (gx - x_min)/dx+0.1f ;
	  //cout<<"is "<< is << "ig "<< ig << "xloc " << xloc <<endl;
	  check( xloc>=0 && xloc<nx , "error: getRecord in inversion.cpp xloc must be within range [0,nx)") ;
	  adj ? memcpy(CSG[ig],rec[xloc],sizeof(float)*nt) : memcpy(rec[xloc],CSG[ig],sizeof(float)*nt) ;
	}
    }
}
void Inversion::swapModel(float ** m, float ** d, float mfx, float dfx, int nzm,int nxm,int nzd,int nxd, bool add)
{
  int layer = 4;
  check(nzm==nzd,"error: SwapModel in inversion.cpp nzm and nzd must be the same. ");
  if(!add)
    opern(d,VALUE,nzd,nxd,0.0f);
  float dx    = param.dx.val;
  float x_min = mfx -0.001*dx  ;
  float x_max = mfx + (nxm - 1) * dx +0.001*dx;
  for(int ix = layer; ix < nxd-layer ; ix ++)
    {
      float xcor = dfx + ( ix - 1 ) * dx ;
      if( xcor >= x_min && xcor <= x_max)
	{
	  int xloc = (xcor - x_min) / dx+0.0001f;
	  check(xloc>=0 && xloc<nxm,"error: SwapModel in inversion.cpp xloc must be within range [0,nxm)");
	  if(xloc>=4 && xloc<=nxm-4)
	  for(int iz = layer; iz < nzd-layer; iz ++)
	    d[ix][iz] += m[xloc][iz];
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
string Inversion::obtainImageName(int is)
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

void Inversion::removeDirect(float ** csg, float ** v, int is)
{
  int ngg  = ng[is];
  int nt   = param.nt.val;
  float sx =  sc[is][0];
  float sz =  sc[is][1];
  float dx = param.dx.val;
  float dz = param.dz.val;
  float dt = param.dt.val;
  float fr = param.fr.val;
  float velfx = param.velfx.val;
  int delay = param.delay.val;
  int  isz = (sz-0.0f)/dz + 0.0001f;
  int  isx = (sx-velfx)/dx + 0.0001f;

  for(int ig = 0; ig < ngg ; ig++)
    {
      float t  = 0.0f;
      float gx = gc[is][0][ig];
      float gz = gc[is][1][ig];
      int igx = (gx -velfx)/dx + 0.0001f;
      int igz = (gz -0.0f )/dz + 0.0001f;
      int adddelay = abs((gz - sz)/dz+0.0001f);
      int  nx = param.nx.val;
      int  nz = param.nz.val;
      adddelay +=10;
      check(igz >= 0 && igz < nz, "igz in removeDirect must be within range [0 nz)");
      if(igx<=isx){
	for(int i=igx;i<isx;i++){
	  if(i>=0 && i<nx)
	    t +=dx/v[i][igz];
	  else if(igx <0)
	    t +=dx/v[0][igz];
	  else if(igx >=nx)
	    t +=dx/v[nx-1][igz];
	}
      }
      else{
	for(int i=igx;i>isx;i--){
	  if(i>=0 && i<nx)
	    t +=dx/v[i][igz];
	  else if(igx <0)
	    t +=dx/v[0][igz];
	  else if(igx >=nx)
	    t +=dx/v[nx-1][igz];
	}
      }
      int it = (t/dt + delay + adddelay+1.0f/fr/dt*0.5f*5.0) + 0.0001f;
      it = it<=(nt-1)?it:nt-1;
      it = it>= 1    ?it:1;
      opern(csg[ig],VALUE,it,0.0f);
  
    }
}
void Inversion::masterRun(int ns)
{
  int is=0; int ii=0;
  int i=0;  int nn=0;
  int idone = 0;
  int itotal= 0;
  int * done = MyAlloc<int>::alc(ns);
  MPI_Status  status;
  int nprocs = this->nprocs;
  int myid = rank;
  int recvd_tag,recvd_source, send_tag ;
  
  for(is=0; is < ns; is++) done[is] = 0;  //done=0, not issued, done=1, issued but not finished, done=2, finished
  
  is = 0;
  while (is < ns   ) {
    MPI_Recv(&is, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    recvd_tag = status.MPI_TAG;
    recvd_source = status.MPI_SOURCE;
    
    
    if( recvd_tag  == 1 ) { assert(is >=0 && is < ns);  done[is] = 2; }     //mark is_th source finished  
    
    //find not issued shot 
    is=0;  while( is < ns && done[is] != 0) is++;     
    
    if(is < ns) {  //find more to do
      done[is] = 1;   //mark as issued
      send_tag = 0;
      MPI_Send(&is, 1, MPI_INT, recvd_source, send_tag, MPI_COMM_WORLD);
      //fprintf(stderr,"send shot %d to slave %d \n", is, recvd_source);
    }
    
    //find not finished shot when all shots were issued already
    if(is >= ns) { is=0; while(is < ns && done[is] != 1) is++; }   //find not finished shot
    
    if(is >=  ns) { //all finished 
      for(is=0; is < ns; is++) {
	//assert(done[is] == 2);
	//fprintf(stderr, "is=%d done=%d \n", is, done[is]);
      }
      nn = 0;
      for(i=1; i < nprocs; i++) {
	MPI_Send(&ns, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	MPI_Recv(&ii, 1, MPI_INT, i,  2, MPI_COMM_WORLD, &status);
	nn += ii; 
	//fprintf(stderr, "node %d has done =%d  total done=%d \n", i, ii, nn );
      }
      // fprintf(stderr,"master=%d Finshed, waiting for exit\n", myid);
    }
    
  } 
  MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MyAlloc<int>::free(done);
}
void Inversion::planeWavePrepare()
{
  if(rank==0)
    {
       
      int nx   = param.nx.val;
      int nz   = param.nz.val;
      int nt   = param.nt.val;
      int delay    = param.delay.val;
      int delaycal = param.delaycal.val;
      int ngmax = param.ngmax.val;
      float dt  = param.dt.val;
      float dx  = param.dx.val;
      float dz  = param.dz.val;
      float fr  = param.fr.val;
      float NT  = nt + delaycal;
      float velfx = param.velfx.val;
      float ** sou=MyAlloc<float>::alc(NT,nx);
      float ** rec=MyAlloc<float>::alc(nt,nx);
      float ** tr =MyAlloc<float>::alc(nt,nx);
      float **CSGN=MyAlloc<float>::alc(nt,ngmax);
      int ns    = param.ns.val;
      int ds    = ns;
      float *** data = MyAlloc<float>::alc(nt,ngmax,ds);
      for(int i = 0; i<1;i++){
	  cout<<"read"<<endl;
	  for(int j=0; j<ds ;j++){
	    int index = j + i*ds;
	    if(index<ns);{
	     cout<<index<<" ,";
	     string file = obtainCSGName(index);
	     read(file,nt,ng[index],data[j]);
	   }
	 }
	 cout<<"read finish"<<endl;
	 float sxmax=sc[0][0];
	 for(int iss=1;iss<param.ns.val;iss++)
	   {
	    sxmax = sxmax<sc[iss][0] ? sc[iss][0]: sxmax; 
	   }
	cout<<"sxmax "<<sxmax <<" ix:" << sxmax/dx << endl;
	  for(int is=0;is<param.np.val;is++)
	    {
	      opern(sou,VALUE,NT,nx,0.0f);
	      opern(rec,VALUE,nt,nx,0.0f);
	      int nx   = param.nx.val;
	      int nz   = param.nz.val;
	      int nt   = param.nt.val;
	      int delay    = param.delay.val;
	      int delaycal = param.delaycal.val;
	      int ngmax = param.ngmax.val;
	      float dt  = param.dt.val;
	      float dx  = param.dx.val;
	      float dz  = param.dz.val;
	      float fr  = param.fr.val;
	      float NT  = nt + delaycal;
	      float velfx = param.velfx.val;
	      float p     = param.pmin.val + is * param.dp.val;
	      float velmin = velfx-0.0001*dx;
	      float velmax = sxmax;
	      float * wav  = rickerWavelet(dt,fr,delay+delaycal,NT);
	      float * tw   = MyAlloc<float>::alc(NT);
	    cout<<"  ik "<< is << "  k_x : "<< p << endl; 
	    for(int ishot = 0; ishot < ds ; ishot ++){
		float sx = sc[ishot][0];
		check((sx-velmin)>=0 && (velmax-sx)>=0,"sx must be in the range [velmin,velmax]");
		float distance = (p>0)? (sx -velmin): (velmax -sx);
		int idts = distance*sqrt(p*p)/dt +delay + delaycal;
		int idtr = distance*sqrt(p*p)/dt ;
		//shiftFFT(wav,tw,NT,idts);
		harmonic(wav, tw, NT, p, sx);
		int isx = (sx - velfx)/dx + 0.001f;
		check(isx>=0 && isx<nx, "isx must be in the range [0 nx) in planeWavePrepare");
		opern(sou[isx],sou[isx],tw,ADD,NT);
		// process receivers
		#ifdef _OPENMP
  (void) omp_set_dynamic(0);
  if (omp_get_dynamic()) {printf("Warning: dynamic adjustment of threads has been set\n");}
  (void) omp_set_num_threads(20);
#endif
#pragma omp parallel for
  for(int igg=0;igg<ng[ishot];igg++){
    harmonic(data[ishot][igg],CSGN[igg],nt,p,sx);
  }
		swapReord( CSGN, ishot, this->ng[ishot], nt, tr, velfx, nx, false );
		opern(rec,tr,rec,ADD,nt,nx);
	      }
	      cout<<"finishe"<<endl;
	      string sfile = obtainNameDat(param.wdir.val,"PLANE_S_NEW",is);
	      string rfile = obtainNameDat(param.wdir.val,"PLANE_R_NEW",is);
	      write(sfile,NT,nx,sou);
	      write(rfile,nt,nx,rec);		
	      MyAlloc<float>::free(wav);
	      MyAlloc<float>::free(tw);
	    }
	 
	}
	  MyAlloc<float>::free(sou);
	  MyAlloc<float>::free(rec);
	  MyAlloc<float>::free(tr) ;
	  MyAlloc<float>::free(CSGN);
	MyAlloc<float>::free(data);
    }
  MPI_Barrier(MPI_COMM_WORLD);
}
 void Inversion::planeWavePrepare_MPI()
 {
	 check(param.planeTag.val,"in param planeTag must be true(none zero) to use subroutine planeWavePrepare");
	 int np     = param.np.val;
	 float dp   = param.dp.val;
	 float pmin = param.pmin.val;
     if(rank==0)
     cout<<" planewavePrepare "<<endl;
     int is=0;
     int idone = 0;
     int itotal= 0;
     int myid = rank;
     int ns  =  np;
     MPI_Status  status;
     if(rank==0)
     masterRun(np);
     else
     {
      is = -1;
      MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            //ask for a shot 
      while (is < ns) 
	{
	  MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);   //receive a new shot
	  if(is >=0 && is < ns){
		  
	    idone++;
	    int nx   = param.nx.val;
	    int nz   = param.nz.val;
	    int nt   = param.nt.val;
	    int npml     = param.npml.val;
	    int delay    = param.delay.val;
	    int delaycal = param.delaycal.val;
	    int ngmax = param.ngmax.val;
	    float dt  = param.dt.val;
	    float dx  = param.dx.val;
	    float dz  = param.dz.val;
	    float fr  = param.fr.val;
	    float NT  = nt + delaycal;
	    float velfx = param.velfx.val;
	    float ** sou=MyAlloc<float>::alc(NT,nx);
	    float ** rec=MyAlloc<float>::alc(nt,nx);
	    float ** tr =MyAlloc<float>::alc(nt,nx);
	    float ** CSG=MyAlloc<float>::alc(nt,ngmax);
		float **CSGN=MyAlloc<float>::alc(nt,ngmax);
	    float p     = pmin + is * dp;
	    float velmin = velfx-0.0001*dx;
	    float velmax = velfx + (nx-1)*dx+0.0001*dx;
	    float * wav  = rickerWavelet(dt,fr,0,NT);
	    float * tw   = MyAlloc<float>::alc(NT);
	    cout<<"  ip "<< is << "  p : "<< p << endl; 
	    for(int ishot = 0; ishot < param.ns.val ; ishot ++)
	    {
			if(rank==1)
			cout<<ishot<<endl;
			// process source
			float sx = sc[ishot][0];
			
			check((sx-velmin)>=0 && (velmax-sx)>=0,"sx must be in the range [velmin,velmax]");
			float distance = (p>0)? (sx -velmin): (velmax -sx);
			int idts = distance*sqrt(p*p)/dt +delay + delaycal;
			int idtr = distance*sqrt(p*p)/dt ;
			//shiftFFT(wav,tw,NT,idts);
			int isx = (sx - velmin)/dx + 0.0001f;
			check(isx>=0 && isx<nx, "isx must be in the range [0 nx) in planeWavePrepare");
			opern(sou[isx],sou[isx],tw,ADD,NT);
			// process receivers
			string CSGfile = obtainCSGName(ishot);
			read(CSGfile,nt,this->ng[ishot],CSG);
			int check = (idtr % (nt*2));
			if( check>=0 && check<nt) ;
			  //shiftSimple(CSG,CSGN,nt,ng[1],check);
			else
				opern(CSGN,VALUE,nt,ng[1],0.0f);
			//shift(CSG,CSG,nt,this->ng[ishot],idtr);
			//for(int ig=0; ig<this->ng[ishot];ig++)
			//{
			//	shiftFFT(CSG[ig],CSG[ig],nt,idtr);
			//}
			swapReord( CSGN, ishot, this->ng[ishot], nt, tr, velfx, nx, false );
			opern(rec,tr,rec,ADD,nt,nx);
		}
		string sfile = obtainNameDat(param.wdir.val,"PLANE_S_NEW",is);
		string rfile = obtainNameDat(param.wdir.val,"PLANE_R_NEW",is);
		//string sfile_su = obtainNameSu(param.wdir.val,"PLANE_S",is);
		//string rfile_su = obtainNameSu(param.wdir.val,"PLANE_R",is);
		write(sfile,NT,nx,sou);
		write(rfile,nt,nx,rec);
		//writeSu(sfile_su,NT,nx,sou);
		//writeSu(rfile_su,nt,nx,rec);
		MyAlloc<float>::free(sou);
		MyAlloc<float>::free(rec);
		MyAlloc<float>::free(tr) ;
		MyAlloc<float>::free(CSG);
		MyAlloc<float>::free(CSGN);
		MyAlloc<float>::free(wav);
		MyAlloc<float>::free(tw);
	    sleep(2); 
	    MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
	  } else {
	    MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	    // fprintf(stderr,"Slave=%d Finshed, waiting for exit \n", myid);
	  }
	}
      MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
 }
void  Inversion::adjointPlane_MPI( float ** img, int migTag){
  string v0file = param.v0file.val;
  float ** v0 = MyAlloc<float>::alc(param.nz.val,param.nx.val);
  read(v0file,param.nz.val,param.nx.val,v0);
  if(rank==0)
    cout<<"adjointPlane"<<endl;
  
  int is=0;
  int idone = 0;
  int itotal= 0;
  int myid = rank;
  int ns  =  param.np.val;
  MPI_Status  status;
  if(rank==0)
    masterRun(ns);
  else{
    is = -1;
    MPI_Send(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);            //ask for a shot 
    while (is < ns) {
      MPI_Recv(&is, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);   //receive a new shot
      if(is >=0 && is < ns){
	idone++;
	int nx   = param.nx.val;
	int nz   = param.nz.val;
	int nt   = param.nt.val;
	int npml     = param.npml.val;
	int delaycal = param.delaycal.val;
	float dt  = param.dt.val;
	float dx  = param.dx.val;
	float dz  = param.dz.val;
	float fr  = param.fr.val;
	string wdir = param.wdir.val;
	int NT = nt + delaycal;
	float ** sou = MyAlloc<float>::alc(NT,nx);
	float ** rec = MyAlloc<float>::alc(nt,nx);
	float ** mask= MyAlloc<float>::alc(nz,nx);
	string recfile;
	string soufile;
	soufile  = obtainNameDat(param.wdir.val,"PLANE_S_NEW",is);
	if(migTag==RTM_IMG)
	  recfile =  obtainNameDat(param.wdir.val,"PLANE_R_NEW",is);
	else
	  recfile = obtainBornName(is);
	read(soufile,NT,nx,sou);
	read(recfile,nt,nx,rec);
	int sz = (sc[0][1]-0.0f)/dz + 0.0001f;
	int gz = (gc[0][1][0]-0.0f)/dz + 0.0001f;
	cout<< " p  "<< param.pmin.val + is * param.dp.val <<" sz "<< sc[0][1] << " gz "<< gc[0][1][0] <<endl;
	float ** imgX = adjointPlane( dt, dx, dz,  nt,  delaycal,  nx,  nz,  npml, sz, gz, v0, sou, rec);
	string imgfile = obtainImageName(is);
	if(param.mask.val){
	  read(param.maskfile.val,nz,nx,mask);
	  opern(imgX,imgX,mask,MUL,nz,nx);
	}
	write(imgfile,nz,nx,imgX);
	MyAlloc<float>::free(sou);
	MyAlloc<float>::free(rec);
	MyAlloc<float>::free(mask);
	MyAlloc<float>::free(imgX);
	sleep(2); 
	MPI_Send(&is, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);    //send finished shot, and ask for another new shot
      } else {
	MPI_Send(&idone, 1, MPI_INT, 0, 2, MPI_COMM_WORLD); 
	// fprintf(stderr,"Slave=%d Finshed, waiting for exit \n", myid);
      }
    }
    MPI_Reduce(&idone, &itotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  MyAlloc<float>::free(v0);
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0){	
    int nx   = param.nx.val;
    int nz   = param.nz.val;
    opern(img,VALUE, nz , nx,0.0f);
    float ** imgtmp = MyAlloc<float>::alc(nz,nx);
    for( int is=0; is<ns; is++){
      string imgfile = obtainImageName(is);
      read(imgfile,nz,nx,imgtmp);
      opern(img,imgtmp,img,ADD,nz,nx);
    }
    MyAlloc<float>::free(imgtmp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&img[0][0], param.nz.val*param.nx.val,MPI_FLOAT,0,MPI_COMM_WORLD); 
}
void Inversion::test()
{
  bool modeling_test=true;
  bool planewavemigrationtest = false;
  // test of the plane wave migration
  if(planewavemigrationtest)
  {
	  int nt = 3000;
	  int nx = 1800*5;
	  int nz = 300;
	  int sz = 4;
	  int gz = 4;
	  float dt = 0.0016;
	  float dx = 0.01;
	  float dz = 0.01;
	  int delay = 10;
	  int delaycal = 100;
	  float fr = 20;
	  int NT = nt + delaycal;
	  float ** sou = MyAlloc<float>::alc(NT,nx);
	  float theta = 80;
	  float **    v= MyAlloc<float>::alc(nz,nx);
	  float **   v0= MyAlloc<float>::alc(nz,nx);
	  float **   dv= MyAlloc<float>::alc(nz,nx);
	  float **  rec= MyAlloc<float>::alc(nt,nx);
	  float * wav = rickerWavelet( dt,  fr, 0, NT);
	  for(int ix = 0; ix<nx;ix++)
	  for(int iz = 0; iz<nz;iz++)
	  {
		  if(iz>=0 && iz<1.0f/2*nz)
		  v[ix][iz] = 2.0;
		  if(iz>=nz/2 && iz <nz*4.0f/5)
		  v[ix][iz] = 2.5f;
		  if(iz>=nz*2.0f/3)
		  v[ix][iz] = 2.5f;
	  }
	  writeSu("v.su",nz,nx,v);
	  Smooth::smooth2d1(v,v0,nz,nx,10);
	  opern(dv,v,v0,SUB,nz,nx);
	  writeSu("v0.su",nz,nx,v0);
	  writeSu("v1.su",nz,nx,dv);
	  float p     = sin(theta*3.1415926/180)/v[0][0];
	  cout<<"p"<< p <<endl;
	  for(int ix = 0; ix<nx;ix++){
		   int mxx = (int)(ix*dx*p/dt) + delay + delaycal;
		   // shiftFFT(wav,sou[ix],NT,mxx);
		   int myyX = (ix*dx*p/dt + delay+delaycal);
		   int myy = myyX %(2*NT);
		   cout<<myyX <<" % " <<myy <<endl;
		   //if(myy>=NT)
		   //opern(sou[ix],VALUE,NT,0.0f);
		   
		  }
	   writeSu("source.su",NT,nx,sou);
	   //OMP_CORE = nthread;
	  rec= forwordPlane(dt, dx, dz, nt, delaycal, nx, nz, 20, sz, gz, v0, sou, dv);
	  writeSu("rec.su",nt,nx,rec);
	  write("rec.dat",nt,nx,rec);
	  read("rec.dat",nt,nx,rec);
	  float ** img = adjointPlane(dt,dx,dz,nt,delaycal,nx,nz,20,sz,gz,v0,sou,rec);
	  writeSu("img.su",nz,nx,img);
  }
  if(modeling_test)
    {
      int nx = param.nx.val;
      int nz = param.nz.val;
      int nt = param.nt.val;
      string vfile = param.vfile.val;
      float ** v = MyAlloc<float>::alc(nz,nx);
      float ** dv= MyAlloc<float>::alc(nz,nx);
      float **img= MyAlloc<float>::alc(nz,nx);
      read(vfile,nz,nx,v);
      Smooth::smooth2d1(v,dv,nz,nx,10);
      opern(dv,v,dv,SUB,nz,nx);
      OMP_CORE = param.nthread.val;
      
      
    /*  if(false){
      cout<<" test of shift "<<endl;
      float ** CSG = MyAlloc<float>::alc(nt,ng[1]);
      float ** CSGN= MyAlloc<float>::alc(nt,ng[1]);
      float *  a= MyAlloc<float>::alc(nt);
      float *  b= MyAlloc<float>::alc(nt);
      read(obtainCSGName(1),nt,ng[1],CSG);
      for(int i=0;i<30000;i+=1000){
	cout<< i <<" "<<ng[1]<<endl;
      //shift(CSG,CSGN,nt,ng[1],i);
      int check = (i % (nt*2));
      if( check>=0 && check<nt);
      //shiftSimple(CSG,CSGN,nt,ng[1],check);
      else
      opern(CSGN,VALUE,nt,ng[1],0.0f);
      string file = obtainNameDat(param.wdir.val,"xx",i);
      write(file,nt,ng[1],CSGN);
      }
  }
	if(false)
	{
		
  string v0file = param.v0file.val;
  float ** v0 = MyAlloc<float>::alc(param.nz.val,param.nx.val);
  read(v0file,param.nz.val,param.nx.val,v0);
			int is = 500;
		int nx   = param.nx.val;
	    int nz   = param.nz.val;
	    int nt   = param.nt.val;
	    int npml     = param.npml.val;
	    int delay    = param.delay.val;
	    int delaycal = param.delaycal.val;
	    int ngmax = param.ngmax.val;
	    float dt  = param.dt.val;
	    float dx  = param.dx.val;
	    float dz  = param.dz.val;
	    float fr  = param.fr.val;
	    // check the velocity model used for each shot
	    string wdir = param.wdir.val;
	    int NT = nt + delaycal;
	    int numdelay = delay + delaycal;
	    float * wav = rickerWavelet(dt,fr,numdelay,NT);
	    float ** CSG = MyAlloc<float>::alc(nt,ngmax);
	    float ** mask= MyAlloc<float>::alc(nz,nx);
	    float velfx  = param.velfx.val;
	    float velfxsub;
	    int  nzsub;
	    int  nxsub;
	    float ** velsub = getVel(velfxsub, nzsub, nxsub, v0, nz, nx, velfx, is);
	    float ** recsub = MyAlloc<float>::alc(nt,nxsub);
	    string recfile;
	      recfile = obtainCSGName(is);
	    read(recfile,nt,ng[is],CSG);
	    removeDirect(CSG,v0,is);
	    swapReord( CSG, is, this->ng[is], nt, recsub, velfxsub, nxsub, false );
	    int * igz = getIgz(velfxsub,nxsub,is);
	}
*/
	  
      
      //forward_MPI(dv);
      //exit(0);
      //modeling_MPI(v);
       //planeWavePrepare();
      //planeWavePrepare_MPI();
       
      // illumPlane_MPI(img);
      // write("illumPlane_simple.dat",nz,nx,img);
       //  adjointPlane_MPI(img,RTM_IMG);
       //writeSu("imagPlane_simple.su",nz,nx,img);
      adjoint_MPI(img,RTM_IMG);
      writeSu("imgShot.su",nz,nx,img);
      // test of moving
      if(false){
      string csgfile = obtainCSGName(1);
      float ** csg = MyAlloc<float>::alc(param.nt.val,ng[1]);
      read(csgfile,param.nt.val,ng[1],csg);
      removeDirect(csg,v,1);
      writeSu("csg_dr.su",param.nt.val,ng[1],csg);
      }
      writeSu("img.su",nz,nx,img);
      MyAlloc<float>::free(v);
      MyAlloc<float>::free(img);
    }
  if(false)
    {
  int nx = param.nx.val;
  int nz = param.nz.val;
  int delaycal = param.delaycal.val;
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
      writeSu(tmp,nzsub,nxsub,velsub);
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
}
