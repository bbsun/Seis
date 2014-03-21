#ifndef INPUT_OUTPUT_H_H
#define INPUT_OUTPUT_H_H
#include <string>
#include <fstream>
#include <iostream>
#include "../util/debug.h"
#include <stdlib.h>
template <typename T>
void read(std::string filename, int n1, T* data)
{
  std::ifstream fr(filename, std::ios::binary );
  if(!fr){
    error("can not open file : ");
    std::cout<<filename<<std::endl;
    exit(0);
  }
  fr.read((char *) &data[0], sizeof(T)* n1);
  fr.close();
}

template <typename T>
void read(std::string filename, int n1, int n2, T**data)
{
  read(filename, n1*n2,&data[0][0]);
}

template <typename T>
void read(std::string filename, int n1, int n2, int n3, T***data)
{
  read(filename, n1*n2*n3,&data[0][0][0]);
}

template <typename T>
void read(std::string filename, int n1, int n2, int n3, int n4, T****data)
{
  read(filename, n1*n2*n3*n4, &data[0][0][0][0]);
}


template <typename T>
void write(std::string filename, int n1, T* data)
{
  std::ofstream fw(filename, std::ios::binary|std::ios::trunc);
  if(!fw){
    error("can not open file: ");
    std::cout<<filename<<std::endl;
    exit(0);
  }
  fw.write((char *) &data[0], sizeof(T)* n1);
  
  fw.close();
}

template <typename T>
void write(std::string filename, int n1, int n2, T**data)
{
  write(filename, n1*n2,&data[0][0]);
}

template <typename T>
void write(std::string filename, int n1, int n2, int n3, T***data)
{
  write(filename, n1*n2*n3,&data[0][0][0]);
}

template <typename T>
void write(std::string filename, int n1, int n2, int n3, int n4, T****data)
{
  write(filename, n1*n2*n3*n4, &data[0][0][0][0]);
}

void writesu(std::string filename, int n1,          float * data);
void writesu(std::string filename, int n1, int n2,  float ** data);
#endif
