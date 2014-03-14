/**
 *@file memory.h
 *@brief memory allocation for 1d, 2d, 3d and 4d array
 *@author Bingbing Sun
 *@version March-13-2014
 */
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <string.h>
using std::endl;
using std::cout;
#ifndef MEMORY_H_H
#define MEMORY_H_H 
template <typename T> class MyAlloc
{
 public:
  static T* alc(size_t n1)
  {
    if (n1 <=0) return NULL;
    T *buf = new T[n1];
    if (buf==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf!=NULL);
    memset(buf,'\0',n1*sizeof(T));
    return buf;
  }
  
  static T** alc(size_t n1, size_t n2)
  {
    if (n1 <=0 || n2 <=0) return NULL;
    T **buf = new T*[n2];
    if (buf==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf!=NULL);
    buf[0] = new T[n1*n2];
    if (buf[0]==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf[0]!=NULL);
    memset(buf[0],'\0',n1*n2*sizeof(T));
    for (size_t i=1;i<n2;i++)
      buf[i]=&buf[0][i*n1];
    return buf;
  }
  static T*** alc(size_t n1, size_t n2, size_t n3)
  {
    if (n1 <=0 || n2 <=0 || n3 <= 0) return NULL;
    T ***buf = new T**[n3];
    if (buf==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf!=NULL);
    buf[0] = new T*[n2*n3];
    if (buf[0]==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf[0]!=NULL);
    buf[0][0] = new T[n1*n2*n3];
    if (buf[0][0]==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf[0][0]!=NULL);
    memset(buf[0][0],'\0',n1*n2*n3*sizeof(T));
    for (size_t i=1;i<n3;i++)
      buf[i]=&buf[0][i*n2];
    for (size_t i=1;i<n2*n3;i++)
      buf[0][i]=&buf[0][0][i*n1];
    return buf;
  }
   
  static T**** alc(size_t n1, size_t n2, size_t n3, size_t n4)
  {
    if (n1 <=0 || n2 <=0 || n3 <= 0 || n4<=0 ) return NULL;
    T ****buf = new T***[n4];
    if (buf==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf!=NULL);
    buf[0] = new T**[n3*n4];
    if (buf[0]==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf[0]!=NULL);
    buf[0][0] = new T*[n2*n3*n4];
    if (buf[0][0]==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf[0][0]!=NULL);
    buf[0][0][0] = new T[n1*n2*n3*n4];
    if (buf[0][0][0]==NULL)
      cout << "Not Enough Memory.\n" <<endl;
    assert(buf[0][0][0]!=NULL);
    memset(buf[0][0][0],'\0',n1*n2*n3*n4*sizeof(T));
    for (size_t i=1;i<n4;i++)
      buf[i]=&buf[0][i*n3];
    for (size_t i=1;i<n3*n4;i++)
      buf[0][i]=&buf[0][0][i*n2];
    for (size_t i=1;i<n2*n3*n4;i++)
      buf[0][0][i]=&buf[0][0][0][i*n1];
    return buf;
  }
  static void free(T* buf)
  {
    if (buf==NULL)
      return;
    delete[] buf;
    buf=NULL;
  }
  static void free(T** buf)
  {
    if (buf==NULL)
      return;
    delete[] buf[0];
    delete[] buf;
    buf==NULL;
  }
  static void free(T*** buf)
  {
    if (buf==NULL)
      return;
    delete[] buf[0][0];
    delete[] buf[0];
    delete[] buf;
    buf==NULL;
  }
  static void free(T**** buf)
  {
    if (buf==NULL)
      return;
    delete[] buf[0][0][0];
    delete[] buf[0][0];
    delete[] buf[0];
    delete[] buf;
    buf==NULL;
  }
};
#endif
