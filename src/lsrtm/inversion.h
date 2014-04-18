/**
 *@file  inversion.h 
 *@brief  seismic inversion.
 *@detail this class deals with seismic inversion problem
 * 1 least square reverse time inversion
 * 2 more to add.
 *@author Bingbing Sun
 *@version April-15-2014
 */
#ifndef INVERSION_H_H
#define INVERSION_H_H
#include "datastr.h"
#include <string>
#define RTM_IMG     1
#define LSRTM_STEP  2
class Inversion
{
 public:
  /**
   *Constructor.
   *@param rank   thread index
   *@param nprocs total number of threads
   */
  Inversion(int rank,int nprocs);
  /**
   *Destructor.
   */
  ~Inversion();
  /**
   * Read parameters for file.
   * All the parameters for inversion are saved in file,
   * we need to read from it to get work initialized. 
   *@param  file  the parameters file
   */
  void readParamFile(std::string file);
  /**
   * Read configurations from a file.
   * After this work done. we will have
   * ng numbers of gather for each shot
   * sc source coordinate for each shot
   * gc gather coordinate for each shot
   */
  void getConfig();
  /**
   * Testing function.
   * We use this function to debug.
   */
  void test();
 private:
  /**
   * synthetic common shot gather modeling.
   *@param [in] v              velocity model used for calculcation 
   */
  void modeling_MPI( float ** v );
  /**
   * illumation compensation for common shot migration
   *@param[out] illumation compensation matrix
   */ 
   void illum_MPI(float ** ss);
  /**
   * mpi implementation of born modeling. 
   *@param [in] dv             velocity perturbation for calculation
   */
  void  forward_MPI( float ** dv );
  /**
   * mpi implementation of RTM. 
   *@param [out] img           RTM image
   *@param migTag              tag for RTM_IMG or a step of LSRTM_STEP
   */
  void  adjoint_MPI( float ** img, int migTag);
  /**
   * mpi implementation of LSRTM operator.
   *@param [in]  in             input array
   *@param [out] out            output array
   */
  void   conope_MPI( float ** in, float**out);
  /**
   * mpi_implementation of LSRTM by conjugate gradient method.
   *@param maxIter              maxIteration number
   *@param reError              relative error for exit
   */
  float ** cgsolver_MPI( int maxIter, float reError); 
  /**
   * we use a dynamic mission pattern to finish our jobs in MPI.
   * this subroutine should be called in the master thread, for example, rank==0
   * we do not need any parameters here, because we have saved all the information in 
   * the private member. 
   */
  void masterRun();
  /**
   * Return the velocity model for finite difference calcuation.
   * returned velocity model velsub will be used in calculation
   * the returned parameters velfxsub and nxsub is useful for
   * insert the record into the calculation model.
   *@param [out] velfxsub       first sample of the output in x-axis
   *@param [out] nzsub          depth sample number of the output
   *@param [out] nxsub          distance sample number of the output
   *@param [in ] vel            velocity model
   *@param       nz             depth sample number of the input
   *@param       nx             distance sample number of the input
   *@param       velfx          first sample of the input in x-axis
   *@param       is             shot number index, in range [0, ns-1]
   *@return                     sub velocty model for calculation
   */
  float ** getVel (float &velfxsub, int &nzsub, int &nxsub, float ** vel,int nz, int nx, float velfx, int is);
  /**
   * Get the depth sample index for each gather of shot is
   *@param velfx                fisrt sample of the sub velocity in x-axis
   *@param nx                   distance smaple number of the sub velocity 
   *@param is                   shot index
   *@return                     a array containing the sample depth index iz[nx]
   */
  int   *  getIgz ( float velfx, int nx ,int is);
  /**
   * cast the grid for calculation to the grid for inversion.
   * swap the two model, the entrance of d will not be zero if its 
   * coordinate is within the range of m.
   *@param [in]  m                    grid for calculation
   *@param [out] d                    grid for inversion
   *@param       mfx                  first sample of m in x-axis
   *@param       dfx                  first sample of d in x-axis
   *@param       nzm                  depth sample number of m
   *@param       nxm                  distance sample number of m
   *@param       nzd                  depth sample number of d
   *@param       nxd                  distance sample number of d
   *@param       add                  clear d or not
   */
  void swapModel(float ** m, float ** d, float mfx, float dfx, int nzm,int nxm,int nzd,int nxd, bool add);
  /**
   * get the record for calcualtion or push back the calculated record to CSG.
   * adjoint = false, CSG --> rec
   * adojint = true , rec --> CSG
   *@param [in] CSG                   commmon shot gather 
   *@param      is                    shot index
   *@param      ng                    number of gathers of CSG
   *@param      nt                    number of sample 
   *@param [in] rec                   rec for calculation 
   *@param      fx                    first sample of rec in x-axis
   *@param      nx                    number of gather for rec
   *@param      adj                   adjoint or not 
   */
  void swapReord(float **CSG, int is, int ng, int nt, float ** rec, float fx, int nx, bool adj);
  /**
   * get the file name for a common shot gather.
   *@param     is                     shot index
   *@return                           file name of the common shot gather       
   */
  std::string obtainCSGName  (int is);
  /**
   * get the file name for a born record gather.
   *@param    is                     shot index
   *@return                          file name of the born recod gather
   */
  std::string obtainBornName (int is);
  /**
   * get the file name for the image.
   *@param   is                     shot index
   *@return                         file name of the image
   */
  std::string obtainImageName(int is);
  /**
   * get a file name for dat binary file. 
   *@param  dir                     directory of the file
   *@param name                     name of the file
   *@param index                    index number
   *@return                         filename
   */
  static std::string obtainNameDat( std::string dir, std::string name,int index);
  /**
   * get a file name for su file.
   *@param dir                      directory of the file
   *@param name                     name of the file
   *@param index                    index number
   *@return                         file name
   */
  static std::string obtainNameSu ( std::string dir, std::string name,int index);
  /**
   * remove direct arrival from the common shot gather.
   *@param csg                      common shot gater, it will be used as both input and output
   *@param v                        velocity model 
   *@param is                       shot index
   */
  void removeDirect(float ** csg, float ** v, int is);
 private:
  int rank;            /**< thread index                                                         */
  int nprocs;          /**< total number of threads                                              */
  paramLSRTM param;    /**< parameters                                                           */
  int     *ng = 0;     /**< number of gathers for each shot                                      */
  float  **sc = 0;     /**< source coordinate for each shot: sc[ns][2]: [ns][0]-> x  [ns][1]-> z */
  float ***gc = 0;     /**<  receiver coordinate for each shot :  sc[ns][2][ngmax]               */
};
#endif
