#ifndef INVERSION_H_H
#define INVERSION_H_H
#include "datastr.h"
#include <string>
class Inversion
{
 public:
  Inversion(int rank,int nprocs);
  ~Inversion();
  void readParamFile(std::string file);
  void getConfig();
  void test();
 private:
  float ** getVel(float &velfxsub, int &nzsub, int &nxsub, float ** vel,int nz, int nx, float velfx, int is);
  void swapModel(float ** m, float ** d, float mfx, float dfx, int nzm,int nxm,int nzd,int nxd, bool add);
  void swapReord(float **CSG, int is, int ng, int nt, float ** rec, float fx, int nx, bool adj);
 private:
  int rank;
  int nprocs;
  paramLSRTM param;
  int     *ng = 0;
  float  **sc = 0; // sc[ns][2]: [ns][0]-> x  [ns][1]-> z
  float ***gc = 0; // sc[ns][2][ngmax]
};
#endif
