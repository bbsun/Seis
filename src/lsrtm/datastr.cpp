#include <iostream>
#include "datastr.h"
#include <stdio.h>
#include "../util/memory.h"
#include "../util/debug.h"
using std::cout;
using std::endl;
using std::string;
paramString::paramString(const paramString & other)
{
  this->key = other.key;
  this->opt = other.opt;
  this->val = other.val;
  this->def = other.def;
}
paramString::paramString(string key, string val, string def, bool opt)
{
  this->key = key;
  this->val = val;
  this->opt = opt;
  if(opt)
  this->def = def;
}
void paramString::print()
{
  cout<<key <<endl;
  cout<<"    val "<< val <<endl;
  cout<<"    def "<< def <<endl;
  cout<<"    opt "<< opt << endl;
}
paramInt::paramInt(const paramInt & other)
{
  this->key = other.key;
  this->opt = other.opt;
  this->val = other.val;
  this->def = other.def;
}
paramInt::paramInt(string key, int val, int def, bool opt)
{
  this->key = key;
  this->val = val;
  this->opt = opt;
  if(opt)
  this->def = def;
}
void paramInt::print()
{
  cout<< key <<endl;
  cout<<"    val "<< val <<endl;
  cout<<"    def "<< def <<endl;
  cout<<"    opt "<< opt << endl;
}
paramFloat::paramFloat(const paramFloat & other)
{
  this->key = other.key;
  this->opt = other.opt;
  this->val = other.val;
  this->def = other.def;
}
paramFloat::paramFloat(string key, float val, float def, bool opt)
{
  this->key = key;
  this->val = val;
  this->opt = opt;
  if(opt)
    this->def = def;
}

void paramFloat::print()
{
  cout<< key <<endl;
  cout<<"    val "<< val <<endl;
  cout<<"    def "<< def <<endl;
  cout<<"    opt "<< opt << endl;
}
void paramLSRTM::print()
{
  ns.print();
  ngmax.print();
  nx.print();
  nz.print();
  nt.print();
  npml.print();
  mask.print();
  maxiter.print();
  delay.print();
  dx.print();
  dz.print();
  dt.print();
  fr.print();
  velfx.print();
  precsg.print();
  wdir.print();
  v0file.print();
  maskfile.print();
  coordfile.print();
}
