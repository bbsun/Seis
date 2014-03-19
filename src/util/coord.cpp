#include "coord.h"
#include <iostream>
#include <fstream>
#include "memory.h"
#include "arraymath.h"
#include <stdio.h>
#include <stdlib.h>
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
Coords::Coords()
{
}
Coords::Coords(const Coords & other)
{
}
Coords::~Coords()
{
}
void Coords::getSources(char * CoordFileName, int ns, float ** &sc)
{
  if(sc==0)
    sc=MyAlloc<float>::alc(2,ns);
  char buffer[256];
  ifstream myfile(CoordFileName);
  if(!myfile){
    cout << "Unalbe to open myfile";
    exit(1);
  }
  int cs=0;
  float hs0,vs0,hg0,vg0;
  float hs1,vs1,hg1,vg1;
  myfile.getline(buffer,256);
  sscanf(buffer,"%f",&hs0);
  sscanf(buffer,"%*s%f",&vs0);
  sscanf(buffer,"%*s%*s%f",&hg0);
  sscanf(buffer,"%*s%*s%*s%f",&vg0);
  sc[0][0]=hs0; sc[0][1]=vs0;
  while(!myfile.eof())
    {
      myfile.getline(buffer,256);
      sscanf(buffer,"%f",&hs1);
      sscanf(buffer,"%*s%f",&vs1);
      sscanf(buffer,"%*s%*s%f",&hg1);
      sscanf(buffer,"%*s%*s%*s%f",&vg1);
      if( hs1==hs0)
	cs +=0;
      else{
	cs++;
	sc[cs][0]=hs1;sc[cs][1]=vs1;
      }
      hs0 = hs1; vs0 = vs1;
    }
  myfile.close();
}
void Coords::getReceivers(char * CoordFileName, int ns, int * ng,float *** &gc)
{
  if(gc==0)
    gc=MyAlloc<float>::alc(max(ng,ns),2,ns);
  char buffer[256];
  ifstream myfile(CoordFileName);
  if(!myfile){
    cout << "Unalbe to open myfile";
    exit(1);
  }
  int cs=0,ig=0;
  float hs0,vs0,hg0,vg0;
  float hs1,vs1,hg1,vg1;
  myfile.getline(buffer,256);
  sscanf(buffer,"%f",&hs0);
  sscanf(buffer,"%*s%f",&vs0);
  sscanf(buffer,"%*s%*s%f",&hg0);
  sscanf(buffer,"%*s%*s%*s%f",&vg0);
  gc[cs][0][ig]=hg0;gc[cs][1][ig]=vg0;
  while(!myfile.eof())
    {
      myfile.getline(buffer,256);
      sscanf(buffer,"%f",&hs1);
      sscanf(buffer,"%*s%f",&vs1);
      sscanf(buffer,"%*s%*s%f",&hg1);
      sscanf(buffer,"%*s%*s%*s%f",&vg1);
      if( hs1==hs0){
	cs +=0;
	ig++;
	gc[cs][0][ig]=hg1;gc[cs][1][ig]=vg1;
      }
      else{
	cs++;
	ig=0;
	gc[cs][0][ig]=hg1;gc[cs][1][ig]=vg1;
      }
      hs0 = hs1; vs0 = vs1;
    }
  myfile.close();
}
void Coords::countReceivers(char * CoordFileName, int ns, int * &ng)
{
  if(ng==0)
    ng=MyAlloc<int>::alc(ns);
  
  ifstream myfile(CoordFileName);
  if(!myfile){
    cout << "Unalbe to open myfile";
    exit(1);
  }
  int cs=0;int ig=0;
  float hs0,vs0,hg0,vg0;
  float hs1,vs1,hg1,vg1;
  char buffer[256];
  myfile.getline(buffer,256);
  sscanf(buffer,"%f",&hs0);
  sscanf(buffer,"%*s%f",&vs0);
  sscanf(buffer,"%*s%*s%f",&hg0);
  sscanf(buffer,"%*s%*s%*s%f",&vg0);
  ig = 1;
  while(!myfile.eof())
    {	
      myfile.getline(buffer,256);
      sscanf(buffer,"%f",&hs1);
      sscanf(buffer,"%*s%f",&vs1);
      sscanf(buffer,"%*s%*s%f",&hg1);	
      sscanf(buffer,"%*s%*s%*s%f",&vg1);	
      if( hs1==hs0){
	cs +=0;
	ig++;
      }
      else{
	cs++;
	ig = 1;
      }
      ng[cs]=ig;
      hs0 = hs1; vs0 = vs1;hg0=hg1;vg0=vg1;
    }
  if(hg0==hg1 && vg0==vg1)
    ng[cs]=ig-1;	
  myfile.close();
}
int Coords::countShots(char * CoordFileName)
{
  char buffer[256];
  ifstream myfile(CoordFileName);
  if(!myfile){
    cout << "Unalbe to open myfile";
    exit(1);
  }
  int cs=1;
  float hs0,vs0,hg0,vg0;
  float hs1,vs1,hg1,vg1;
  myfile.getline(buffer,256);
  sscanf(buffer,"%f",&hs0);
  sscanf(buffer,"%*s%f",&vs0);
  sscanf(buffer,"%*s%*s%f",&hg0);
  sscanf(buffer,"%*s%*s%*s%f",&vg0);
  while(!myfile.eof())
    {
      myfile.getline(buffer,256);
      sscanf(buffer,"%f",&hs1);
      sscanf(buffer,"%*s%f",&vs1);
      sscanf(buffer,"%*s%*s%f",&hg1);
      sscanf(buffer,"%*s%*s%*s%f",&vg1);
      if( hs1==hs0)
	cs +=0;
      else{
	cs++;
      }
      hs0 = hs1; vs0 = vs1;
    }
  myfile.close();
  return cs;
}
void Coords::printShotInfo(int BeginShot,int EndShot,char *CoordFileName,char* OutputFileName)
{
  ifstream myfile(CoordFileName);
  Coords * cord = new Coords();
  int ns = cord->countShots(CoordFileName);
  cout << ns<<endl;
  int * ng=0;
  cord->countReceivers(CoordFileName,ns,ng);
  float *** gc=0;
  float ** sc=0;
  cord->getSources(CoordFileName,ns,sc);
  cord->getReceivers(CoordFileName,ns,ng,gc);
  ofstream *outfile=0;
  if(OutputFileName !=0)
    outfile=new ofstream(OutputFileName);
  if(OutputFileName!=0){
    if(!outfile){
      cout << "Unalbe to open myfile";
      exit(1);
    }
  }
  if(OutputFileName!=0){
    *outfile<<"Total Shot :"<<ns<<endl;
    for(int i=BeginShot-1;i<=EndShot-1;i++){
      *outfile<<"Shot Num: "<<i+1<<" x: "<<sc[i][0]<<" z: "<<sc[i][1]<<"num of traces"<<ng[i]<<endl;
      for(int j=0;j<ng[i];j++)
	{
	  *outfile<<"receiver Num: "<<j+1<<" x: "<<gc[i][0][j]<<" z: "<<gc[i][1][j]<<endl;
	}
    }
  }
  else
    {
      cout <<"Total Shot :"<<ns<<endl;
      for(int i=BeginShot-1;i<=EndShot-1;i++){
	cout<<"Shot Num: "<<i+1<<" x: "<<sc[i][0]<<" z: "<<sc[i][1]<<"num of traces"<<ng[i]<<endl;
	for(int j=0;j<ng[i];j++)
	  {
	    cout<<"receiver Num: "<<j+1<<" x: "<<gc[i][0][j]<<" z: "<<gc[i][1][j]<<endl;
	  }
      }
    }
  if(outfile)
    {
      outfile->close();
      delete outfile;
    }
  if(cord !=0) delete cord;
  MyAlloc<int>::free(ng);
  MyAlloc<float>::free(gc);
  MyAlloc<float>::free(sc);
  
}
