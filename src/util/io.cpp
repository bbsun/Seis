#include "io.h"
#include "mysu.h"
#include <string>
#include <iostream>
#include <string.h>
using std::string  ;
using std::ofstream;
void writeSu(string filename,int n1,         float * data)
{
  suheader header;
  memset(&header,'\0',sizeof(header));
  header.d1 = 1;
  header.d2 = 1;
  header.f1 = 0;
  header.ns = n1;
  header.ep = 1;
  header.fldr = 1;
  header.dt = 1000000;
  header.tracl = 0;
  header.tracr = 0;
  header.f1    = 0;
  header.f2    = 1;
  std::ofstream fw(filename, std::ios::binary|std::ios::trunc);
  
  header.tracl ++;
  header.tracr ++;
  header.cdp   = 1;
  header.cdpt  = 1;
  fw.write((char*)&header,sizeof(header));
  fw.write((char*)&data[0],sizeof(float)*n1);
  
  fw.close();
}
void writeSu(string filename,int n1,int n2,  float ** data)
{
  suheader header;
  memset(&header,'\0',sizeof(header));
  header.d1 = 1;
  header.d2 = 1;
  header.ns = n1;
  header.ep = 1;
  header.fldr = 1;
  header.dt =   1000;
  header.tracl = 0;
  header.tracr = 0;
  header.f1    = 0;
  header.f2    = 1;
  std::ofstream fw(filename, std::ios::binary|std::ios::trunc);
  for(int i=0; i<n2;i++)
    {
      header.tracl = header.tracl + 1;
      header.tracr = header.tracr + 1;
      header.cdp   = i+1;
      header.cdpt  = i+1;
      fw.write((char*)&header,sizeof(header));
      fw.write((char*)&(data[i][0]),sizeof(float)*n1);
    }
}
