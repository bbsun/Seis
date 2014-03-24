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
  float ** getVel (float &velfxsub, int &nzsub, int &nxsub, float ** vel,int nz, int nx, float velfx, int is);
  int   *  getIgz ( float velfx, int nx ,int is);
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
  void swapModel(float ** m, float ** d, float mfx, float dfx, int nzm,int nxm,int nzd,int nxd, bool add);
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
  void swapReord(float **CSG, int is, int ng, int nt, float ** rec, float fx, int nx, bool adj);
  
  std::string obtainCSGName  (int is);
  std::string obtainBornName (int is);
  std::string obtianImageName(int is);
  static std::string obtainNameDat( std::string dir, std::string name,int index);
  static std::string obtainNameSu ( std::string dir, std::string name,int index);
 private:
  int rank;
  int nprocs;
  paramLSRTM param;
  int     *ng = 0;
  float  **sc = 0; // sc[ns][2]: [ns][0]-> x  [ns][1]-> z
  float ***gc = 0; // sc[ns][2][ngmax]
};
#endif
