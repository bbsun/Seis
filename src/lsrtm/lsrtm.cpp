#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "../util/memory.h"
#include "../util/arraymath.h"
#include "../util/coord.h"
#include "../util/parser.h"
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
}
