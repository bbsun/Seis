#include "io.h"
#include "mysu.h"
#include <string>
using std::string  ;
using std::ofstream;
void writesu(string filename,int n1,         float * data)
{
  suheader header;
  header.d1 = 1;
  header.f1 = 0;
  header.ns = n1;
  header.ep = 1;
  header.fldr = 1;
  header.dt = 1000000;
  header.tracl = 0;
  header.tracr = 0;
  std::ofstream fw(filename, std::ios::binary );
  
  header.tracl ++;
  header.tracr ++;
  header.cdp   = 1;
  header.cdpt  = 1;
  fw.write((char*)&header,sizeof(header));
  fw.write((char*)&data[0],sizeof(float)*n1);
  
  fw.close();
}
void writesu(string filename,int n1,int n2,  float ** data)
{
  suheader header;
  header.d1 = 1;
  header.f1 = 0;
  header.ns = n1;
  header.ep = 1;
  header.fldr = 1;
  header.dt = 1000000;
  header.tracl = 0;
  header.tracr = 0;
  std::ofstream fw(filename, std::ios::binary );
  for(int i=0; i<n2;i++)
    {
      header.tracl ++;
      header.tracr ++;
      header.cdp   = i+1;
      header.cdpt  = i+1;
      fw.write((char*)&header,sizeof(header));
      fw.write((char*)&data[i][0],sizeof(float)*n1);
    }
}
