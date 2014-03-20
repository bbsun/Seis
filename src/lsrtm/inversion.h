#ifndef INVERSION_H_H
#define INVERSION_H_H
#include "datastr.h"
#include <string>
class Inversion
{
 public:
  Inversion(int rank,int nprocs);
  void readParamFile(std::string file);
  void getConfig();
  ~Inversion();
 private:
  int rank;
  int nprocs;
  paramLSRTM param;
  int     *ng = 0;
  float  **sc = 0; // sc[ns][2]: [ns][0]-> x  [ns][1]-> z
  float ***gc = 0; // sc[ns][2][ngmax]
};
#endif
