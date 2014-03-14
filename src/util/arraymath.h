#ifndef ARRAY_MATH_H_H
#define ARRAY_MATH_H_H
#define COPY 1
#define SCAL 2
#define ADD  3
#define SUB  4
#define MUL  5
#define DIV  6
#define ADDSCAL 7
#define SQRT 8
#define DOTPRODUCT 10001
#define VALUE      20001
#define RANDOM       20002
#define  ERROR_STRING "input or output is none for COPY"
#include <string.h>
#include "debug.h"
#include <stdlib.h>
#include <math.h>
template <typename T>
T min(T* in, int n)
{
  check(n>0,"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T minF = in[0];
  for(int i=1;i<n;i++)
    minF = minF>in[i]?in[i]:maxF;
  return minF;
}

template <typename T>
T min(T** in, int n1, int n2)
{
  int n = n1 *n2;
  check((n1>0) &&(n2>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T minF = in[0][0];
  for(int i=1;i<n;i++)
    minF = minF<in[0][i]?in[0][i]:maxF;
  return minF;
}

template <typename T>
T min(T*** in, int n1, int n2, int n3)
{
  int n = n1 *n2 * n3;
  check((n1>0) &&(n2>0) &&(n3>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T minF = in[0][0][0];
  for(int i=1;i<n;i++)
    minF = minF<in[0][0][i]?in[0][0][i]:maxF;
  return minF;
}

template <typename T>
T min(T*** in, int n1, int n2, int n3, int n4)
{
  int n = n1 *n2 * n3 * n4;
  check((n1>0) &&(n2>0) &&(n3>0) &&(n4>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T minF = in[0][0][0][0];
  for(int i=1;i<n;i++)
    minF = minF<in[0][0][0][i]?in[0][0][0][i]:maxF;
  return minF;
}

template <typename T>
T max(T* in, int n)
{
  check(n>0,"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T maxF = in[0];
  for(int i=1;i<n;i++)
    maxF = maxF<in[i]?in[i]:maxF;
  return maxF;
}

template <typename T>
T max(T** in, int n1, int n2)
{
  int n = n1 *n2;
  check((n1>0) &&(n2>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T maxF = in[0][0];
  for(int i=1;i<n;i++)
    maxF = maxF<in[0][i]?in[0][i]:maxF;
  return maxF;
}

template <typename T>
T max(T*** in, int n1, int n2, int n3)
{
  int n = n1 *n2 * n3;
  check((n1>0) &&(n2>0) &&(n3>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T maxF = in[0][0][0];
  for(int i=1;i<n;i++)
    maxF = maxF<in[0][0][i]?in[0][0][i]:maxF;
  return maxF;
}

template <typename T>
T max(T*** in, int n1, int n2, int n3, int n4)
{
  int n = n1 *n2 * n3 * n4;
  check((n1>0) &&(n2>0) &&(n3>0) &&(n4>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  T maxF = in[0][0][0][0];
  for(int i=1;i<n;i++)
    maxF = maxF<in[0][0][0][i]?in[0][0][0][i]:maxF;
  return maxF;
}


template <typename T>
void opern(T* in, int tag, int n, T val=0)
{
  check(n>0,"size of array incorrect");
  check(in!=0, ERROR_STRING);
  switch (tag){
  case VALUE:
    for(int i=0;i<n;i++)
      in[i] = val;
    break;
  case RANDOM:
    for(int i=0;i<n;i++)
      in[i] = (T) rand()/RAND_MAX;
    break;
  default:
    check(false,"operation is not defined");
  }
}

template <typename T>
void opern(T** in, int tag, int n1,int n2, T val=0)
{
  int n= n1*n2;
  check((n1>0) &&(n2>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  switch (tag){
  case VALUE:
    for(int i=0;i<n;i++)
      in[0][i] = val;
    break;
  case RANDOM:
    for(int i=0;i<n;i++)
      in[0][i] = (T) rand()/RAND_MAX;
    break;
  default:
    check(false,"operation is not defined");
  }
}

template <typename T>
void opern(T*** in, int tag, int n1,int n2, int n3,T val=0)
{
  int n= n1*n2*n3;
  check((n1>0) &&(n2>0) &&(n3>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  switch (tag){
  case VALUE:
    for(int i=0;i<n;i++)
      in[0][0][i] = val;
    break;
  case RANDOM:
    for(int i=0;i<n;i++)
      in[0][0][i] = (T) rand()/RAND_MAX;
    break;
  default:
    check(false,"operation is not defined");
  }
}

template <typename T>
void opern(T**** in, int tag, int n1,int n2, int n3, int n4,T val=0)
{
  int n= n1*n2*n3*n4;
  check((n1>0) &&(n2>0) &&(n3>0) &&(n4>0),"size of array incorrect");
  check(in!=0, ERROR_STRING);
  switch (tag){
  case VALUE:
    for(int i=0;i<n;i++)
      in[0][0][0][i] = val;
    break;
  case RANDOM:
    for(int i=0;i<n;i++)
      in[0][0][0][i] = (T) rand()/RAND_MAX;
    break;
  default:
    check(false,"operation is not defined");
  }
}

template <typename T>
double opern(T* in1, T* in2, int tag, int n)
{
  check(n>0,"size of array incorrect");
  double out=0;
  switch (tag){
  case DOTPRODUCT:
    check(in1!=0 && in2!=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      out +=(double)(in1[i]*in2[i]);
    break;
  default:
    check(false,"operation is not defined");
    }
  return out;
}

template <typename T>
double opern(T** in1, T** in2, int tag, int n1, int n2)
{
  int n= n1*n2;
  check((n1>0) && (n2>0),"size of array incorrect");
  double out=0;
  switch (tag){
  case DOTPRODUCT:
    check(in1!=0 && in2!=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      out +=(double)(in1[0][i]*in2[0][i]);
    break;
  default:
    check(false,"operation is not defined");
    }
  return out;
}

template <typename T>
double opern(T*** in1, T*** in2, int tag, int n1, int n2, int n3)
{
  int n= n1*n2*n3;
  check((n1>0) && (n2>0) && (n3>0),"size of array incorrect");
  double out=0;
  switch (tag){
  case DOTPRODUCT:
    check(in1!=0 && in2!=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      out +=(double)(in1[0][0][i]*in2[0][0][i]);
    break;
  default:
    check(false,"operation is not defined");
    }
  return out;
}

template <typename T>
double opern(T*** in1, T*** in2, int tag, int n1, int n2, int n3, int n4)
{
  int n= n1*n2*n3*n4;
  check((n1>0) && (n2>0) && (n3>0) && (n4>0), "size of array incorrect");
  double out=0;
  switch (tag){
  case DOTPRODUCT:
    check(in1!=0 && in2!=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      out +=(double)(in1[0][0][0][i]*in2[0][0][0][i]);
    break;
  default:
    check(false,"operation is not defined");
  }
  return out;
}

template <typename T>
void opern(T* des,T*src1, T* src2,int tag, int n, double a1=0, double a2=0)
{
  check(n>0,"size of array incorrect");
  switch (tag){
  case COPY:
    check(des!=0 && src1 !=0, ERROR_STRING);
      memcpy(des,src1,sizeof(T)*n);
    break;
  case SCAL:
    check( (des!=0) && (src1 !=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[i] =a1*src1[i];
    break;
  case ADD:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[i] =src1[i]+ src2[i];
    break;
  case SUB:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[i] =src1[i]-src2[i];
    break;
  case MUL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[i] =src1[i]*src2[i];
    break;
  case DIV:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      check(src2[i]>0,"divide zero");
      des[i] =src1[i]/src2[i];
    }
    break;
  case ADDSCAL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      des[i] = a1 * src1[i] + a2 * src2[i];
    }
    break;
  case SQRT:
    check(des!=0 && src1 !=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      des[i] = sqrt(src1[i]);
    break;
  default:
    check(false,"operation is not defined");
    
  }
}
template <typename T>
void opern(T** des,T**src1, T**src2,int tag, int n1,int n2, double a1=0, double a2=0)
{
  int n= n1*n2;
  check(n>0,"size of array incorrect");
  switch (tag){
  case COPY:
    check(des!=0 && src1 !=0, ERROR_STRING);
    memcpy(des[0],src1[0],sizeof(T)*n);
    break;
  case SCAL:
    check( (des!=0) && (src1 !=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][i] =a1*src1[0][i];
    break;
  case ADD:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][i] =src1[0][i]+ src2[0][i];
    break;
  case SUB:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][i] =src1[0][i]-src2[0][i];
    break;
  case MUL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][i] =src1[0][i]*src2[0][i];
    break;
  case DIV:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      check(src2[0][i]>0,"divide zero");
      des[0][i] =src1[0][i]/src2[0][i];
    }
    break;
  case ADDSCAL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      des[0][i] = a1 * src1[0][i] + a2 * src2[0][i];
    }
    break;
  case SQRT:
    check(des!=0 && src1 !=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][i] = sqrt(src1[0][i]);
    break;
  default:
    check(false,"operation is not defined");
  }
}
template <typename T>
void opern(T*** des,T***src1, T***src2,int tag, int n1,int n2,int n3,  double a1=0, double a2=0)
{
  int n= n1*n2*n3;
  check(n>0,"size of array incorrect");
  switch (tag){
  case COPY:
    check(des!=0 && src1 !=0, ERROR_STRING);
    memcpy(des[0][0],src1[0][0],sizeof(T)*n);
    break;
  case SCAL:
    check( (des!=0) && (src1 !=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][i] =a1*src1[0][0][i];
    break;
  case ADD:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][i] =src1[0][0][i]+ src2[0][0][i];
    break;
  case SUB:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][i] =src1[0][0][i]-src2[0][0][i];
    break;
  case MUL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][i] =src1[0][0][i]*src2[0][0][i];
    break;
  case DIV:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      check(src2[0][0][i]>0,"divide zero");
      des[0][0][i] =src1[0][0][i]/src2[0][0][i];
    }
    break;
  case ADDSCAL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      des[0][0][i] = a1 * src1[0][0][i] + a2 * src2[0][0][i];
    }
    break;
  case SQRT:
    check(des!=0 && src1 !=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][i] = sqrt(src1[0][0][i]);
    break;
  default:
    check(false,"operation is not defined");
  }
}
template <typename T>
void opern(T**** des,T****src1, T****src2,int tag, int n1,int n2,int n3,int n4,double a1=0, double a2=0)
{
  int n= n1*n2*n3*n4;
  check(n>0,"size of array incorrect");
  switch (tag){
  case COPY:
    check(des!=0 && src1 !=0, ERROR_STRING);
    memcpy(des[0][0][0],src1[0][0][0],sizeof(T)*n);
    break;
  case SCAL:
    check( (des!=0) && (src1 !=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][0][i] =a1*src1[0][0][0][i];
    break;
  case ADD:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][0][i] =src1[0][0][0][i]+ src2[0][0][0][i];
    break;
  case SUB:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][0][i] =src1[0][0][0][i]-src2[0][0][0][i];
    break;
  case MUL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][0][i] =src1[0][0][0][i]*src2[0][0][0][i];
    break;
  case DIV:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      check(src2[0][0][0][i]>0,"divide zero");
      des[0][0][0][i] =src1[0][0][0][i]/src2[0][0][0][i];
    }
    break;
  case ADDSCAL:
    check( (des!=0) && (src1 !=0) && (src2!=0),ERROR_STRING);
    for(int i=0;i<n;i++){
      des[0][0][0][i] = a1 * src1[0][0][0][i] + a2 * src2[0][0][0][i];
    }
    break;
  case SQRT:
    check(des!=0 && src1 !=0, ERROR_STRING);
    for(int i=0;i<n;i++)
      des[0][0][0][i] = sqrt(src1[0][0][0][i]);
    break;
  default:
    check(false,"operation is not defined");
  }
}
#endif
