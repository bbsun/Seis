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
#define  ERROR_STRING "input or output is none for COPY"
#include <string.h>
#include "debug.h"
#include <stdlib.h>
#include <math.h>
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
