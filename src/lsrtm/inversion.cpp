#include "mpi.h"
#include "inversion.h"
#include "../util/parser.h"
#include "../util/coord.h"
#include "../util/arraymath.h"
#include "../util/debug.h"
#include "../util/memory.h"
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
  cout<<ns<<endl;
  cout<<ngmax<<endl;
  
  MPI_Bcast(&ng[0],       ns,        MPI_INT,   0, MPI_COMM_WORLD); 
  MPI_Bcast(&sc[0][0],    ns*2,      MPI_FLOAT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&gc[0][0][0], ns*2*ngmax,MPI_FLOAT, 0, MPI_COMM_WORLD); 
  if(rank==2 && verbose)
    {
      Coords::printShotInfo(1,10, ns,ngmax,this->ng,this->sc,this->gc);
      Coords::printShotInfo(ns-10,ns,ns,ngmax,this->ng,this->sc,this->gc);
    }
}
