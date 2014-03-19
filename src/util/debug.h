#ifndef DEGUG_H_H
#include <iostream>
#include <string>
void check(bool eval, char * str);
void error(std::string errinfo);
template <typename T>
void dump(char* name,T* a , int n)
{
  std::cout<<" "<<name;
  for(int i=0; i< n; i++)
    std::cout<<" "<<""<<" "<<""<<a[i]<<",";
  std::cout<<std::endl;
}
template <typename T>
void dump(char * name, T** a, int n1, int n2)
{
  for(int i=0;i<n2;i++){
    cout<<" line : "<<i;
    dump(name,a[i],n1);
  }
}
template <typename T>
void dump(char * name, T*** a, int n1, int n2, int n3)
{
  for(int i=0;i<n3;i++){
    cout<<"slice : "<<i<<endl;
    dump(name,a[i3],n1,n2);
  }
}
#endif
