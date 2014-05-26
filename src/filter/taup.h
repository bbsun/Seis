/**
 *@file taup.h 
 *@brief taup transform and its adjoint.
 *@author Bingbing Sun
 *@version May-22-2014
 */
#ifndef TAUP_H_H
#define TAUP_H_H
#include <vector>
using std::vector;
/**
 *Taup transform.
 *This is taup transform applied in the time domain.
 *We have two version of taup transform:
 *one saves the tempoary value and need large memory while
 *the remains gives direct implemenation.Following show to use
 * this class. 
 * @code
 // large temporay memoray required 
 Taup taup(dt,nt,x0,dx,nx,p0,dp,np);
 taup.apply(xt,xp,false); // forward taup transform
 taup.apply(xt,xp,true ); // adjoint taup transform
 // on fly calculation, don't need large temporay memoray 
 Taup::apply(xt,xp,false,dt,nt,x0,dx,nx,p0,dp,np); // forward taup transform
 Taup::apply(xt,xp,true ,dt,nt,x0,dx,nx,p0,dp,np); // adjoint taup transform
 @endcode
 */
class Taup
{
 public:
  /**
   *Constructor.
   *Needs large memory: approx. nx*np*nt and it would be efficient 
   * when calculating local taup transform
   *@param dt time sample step
   *@param nt number of time sample
   *@param x0 first sample of x 
   *@param dx x sample step
   *@param nx number of sample of x
   *@param p0 first sample of p
   *@param dp p sample step
   *@param np number of sample of p 
   */
  Taup(float dt, int nt, float x0,float dx, int nx, float p0,float dp, int np);
  /**
   * Deconstructor.
   */
  ~Taup();
  /**
   *Apply taup transform.
   *@param tx data in t-x domain
   *@param tp data in t-p domain
   *@param adj false for taup transform, true for adjoint taup transform 
   */
  void apply(float ** tx,float ** tp,bool adj);
  /**
   *Apply taup transform . 
   *Implement directly without saving temporay variable: v_iit, v_it , v_ip, v_aa
   *@param tx data in t-x domain
   *@param tp data in t-p domain
   *@param adj false for taup transform, true for adjoint taup transform
   *@param dt time sample step
   *@param nt number of time sample
   *@param x0 first sample of x
   *@param dx x sample step
   *@param nx number of sample of x
   *@param p0 first sample of p
   *@param dp p sample step
   *@param np number of sample of p 
   */
  static void apply(float ** tx, float **tp, bool adj, float dt, int nt, float x0, float dx, int nx, float p0, float dp, int np);
 private:
  /**
   * Initialize the temporary variable : v_iit, v_it, v_ip, v_aa
   */
  void init();  
  vector<int>   * v_iit=0 ;  /**< temporay variable */ 
  vector<int>   * v_it=0 ;  /**< temporay variable */
  vector<int>   * v_ip=0 ;  /**< temporay variable */ 
  vector<float> * v_aa=0 ;  /**< temporay varibale */
  float dt;                /**< time sample step  */
  int   nt;                /**< number of time sample */
  float x0;                /**< first sample of x     */
  float dx;                /**< x sample step         */
  int   nx;                /**< number of sample of x */
  float p0;                /**< first sample of     p */
  float dp;                /**< p sample step         */
  int   np;                /**< number of sample of p */
};
#endif
