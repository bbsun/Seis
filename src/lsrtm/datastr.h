#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H
#define DEFSTRING "**"
#define DEFINT    -1
#define DEFFLOAT  -9.9
#include <string>
#define DEFAULT_NPML 10
#define DEFAULT_MASK 0
#define DEFAULT_MAXITER 30

class paramString
{
 public:
  paramString(const paramString & other);
  paramString(std::string key, std::string val, std::string def, bool opt);
  void print();
  
 public:
  std::string key;
  std::string val;
  std::string def = DEFSTRING;
  bool        opt = false;
};

class  paramInt
{
 public:
  paramInt(const paramInt & other );
  paramInt(std::string key, int val, int def, bool opt);
  void print();
 public:
  std::string key;
  int    val;
  int    def = DEFINT;
  bool   opt = false;
};

class paramFloat
{
 public:
  paramFloat(const paramFloat & other);
  paramFloat(std::string key, float val, float def, bool opt);
  void print();
 public:
  std::string key;
  float  val;
  float def = DEFFLOAT;
  bool  opt = false;
};

class paramLSRTM
{
 public:
  void print();
 public:
  paramInt   planeTag{paramInt("planeTag",0,0,1)};
  paramInt   nthread{paramInt("nthread",0,1,1)};
  paramInt   ns{paramInt("ns",0,0,0)},   ngmax{paramInt("ngmax",0,0,0)};
  paramInt   nx{paramInt("nx",0,0,0)},   nz{paramInt("nz",0,0,0)},   nt{paramInt("nt",0,0,0)},   npml{paramInt("npml",0,DEFAULT_NPML,1)} ;
  paramInt   mask{paramInt("mask",0,DEFAULT_MASK,1)}, maxiter{paramInt("maxiter",0,DEFAULT_MAXITER,1)}, delay{paramInt("delay",0,0,1)}, delaycal{paramInt("delaycal",0,0,1)};
  paramInt   lpad{paramInt("lpad",0,0,1)}, rpad{paramInt("rpad",0,0,1)};
  paramFloat dx{paramFloat("dx",0,0,0)}, dz{paramFloat("dz",0,0,0)}, dt{paramFloat("dt",0,0,0)}, fr{paramFloat("fr",0,0,0)};
  paramFloat velfx{paramFloat("velfx",0,0,1)};
  paramFloat pmin{paramFloat("pmin",0,0,0)};
  paramFloat dp{paramFloat("dp",0,0,0)};
  paramInt   np{paramInt("np",0,0,0)};
  paramString precsg{paramString("precsg","","",0)};
  paramString wdir{paramString("wdir","","",0)};
  paramString v0file{paramString("v0file","","",0)};
  paramString vfile{paramString("vfile","","",1)};
  paramString maskfile{paramString("maskfile","","",1)};
  paramString coordfile{paramString("coordfile","","",0)};
};
#endif
