#include "taup.h"
#include <math.h>
#include <iostream>
#include "../util/arraymath.h"
#include "../util/debug.h"
using std::cout;
using std::endl;
Taup::Taup(float dt, int nt, float x0,float dx, int nx, float p0,float dp, int np)
{
  this->dt = dt;
  this->nt = nt;
  this->x0 = x0;
  this->dx = dx;
  this->nx = nx;
  this->p0 = p0;
  this->dp = dp;
  this->np = np;
  v_it  = new vector<int>[nx];  check(v_it!=0,  "can not init v_it  in Taup::Taup()");
  v_ip  = new vector<int>[nx];  check(v_ip!=0,  "can not init v_ip  in Taup::Taup()");
  v_iit = new vector<int>[nx];  check(v_iit!=0, "can not init v_iit in Taup::Taup()");
  v_aa  = new vector<float>[nx];check(v_aa !=0, "can not init v_aa  in Taup::Taup()");
  init();
}
Taup::~Taup()
{
	for(int ix=0;ix<nx;ix++){
		 v_it[ix].clear();
		 v_ip[ix].clear();
		v_iit[ix].clear();
		 v_aa[ix].clear();
	}
	delete [] v_it;
	delete [] v_ip;
	delete [] v_iit;
	delete [] v_aa;
}
void Taup::init()
{
  for(int ip=0;ip<np;ip++){
    float p = p0 + ip*dp;
    for(int ix=0;ix<nx;ix++){
      float x = x0+ix*dx;
      for(int it=0;it<nt;it++){
	float t  = it*dt;
	float tt = (t + p*x)/dt;
	int iit  = floor(tt);
	float aa = tt - iit;
	if(iit>=0 && iit<nt-1){
	v_it[ix].push_back(it);
	v_ip[ix].push_back(ip);
	v_iit[ix].push_back(iit);
	v_aa[ix].push_back(aa);
	}
      }
    }
  }
}
void Taup::apply(float ** tx,float ** tp,bool adj, bool add)
{
  if(!add){
  if(!adj)
    opern(tp,VALUE,nt,np,0.0f);
  else
    opern(tx,VALUE,nt,nx,0.0f);
  }
  if(!adj){
	  cout<<"forward taup transform"<<endl;
    for(int ix=0;ix<nx;ix++){
      int c_size = v_iit[ix].size();
      for(int ic = 0;ic < c_size; ic++){
	int it   =  v_it[ix][ic];
	int ip   =  v_ip[ix][ic];
	int iit  =  v_iit[ix][ic];
	float aa =  v_aa[ix][ic];
	tp[ip][it]  += tx[ix][iit]*(1.0f-aa) + tx[ix][iit+1]*aa;    
      }
    }
  }
  else{
	  cout<<"adjoint taup transform"<<endl;
    for(int ix=0;ix<nx;ix++){
      int c_size = v_iit[ix].size();
      for(int ic = 0;ic < c_size; ic++){
	int it   =  v_it[ix][ic];
	int ip   =  v_ip[ix][ic];
	int iit  =  v_iit[ix][ic];
	float aa =  v_aa[ix][ic];
	tx[ix][iit]   += tp[ip][it]*(1.0f - aa);
	tx[ix][iit+1] += tp[ip][it]*aa;    
      }
    }
  }
}
void Taup::apply(float ** tx, float **tp, bool adj, bool add, float dt, int nt, float x0, float dx, int nx, float p0, float dp, int np)
{
  if(!add){
   if(!adj)
    opern(tp,VALUE,nt,np,0.0f);
  else
    opern(tx,VALUE,nt,nx,0.0f);
  }
    if(!adj){
		cout<<"forward taup transform"<<endl;
		  for(int ip=0;ip<np;ip++){
    float p = p0 + ip*dp;
    for(int ix=0;ix<nx;ix++){
      float x = x0+ix*dx;
      for(int it=0;it<nt;it++){
	float t  = it*dt;
	float tt = (t + p*x)/dt;
	int iit  = floor(tt);
	float aa = tt - iit;
	if(iit>=0 && iit<nt-1){
		tp[ip][it]  += tx[ix][iit]*(1.0f-aa) + tx[ix][iit+1]*aa;}
		}
		}
		}
		}
	else{
		cout<<"adjoint taup transform"<<endl;
		for(int ip=0;ip<np;ip++){
    float p = p0 + ip*dp;
    for(int ix=0;ix<nx;ix++){
      float x = x0+ix*dx;
      for(int it=0;it<nt;it++){
	float t  = it*dt;
	float tt = (t + p*x)/dt;
	int iit  = floor(tt);
	float aa = tt - iit;
	if(iit>=0 && iit<nt-1){
		tx[ix][iit]   += tp[ip][it]*(1.0f - aa);
		tx[ix][iit+1] += tp[ip][it]*aa;}
		}
		}
		}
		}
}
