#include <iostream>
#include "datastr.h"
#include <stdio.h>
#include <iomanip>
#include "../util/memory.h"
#include "../util/debug.h"
using std::cout;
using std::endl;
using std::string;
using std::setw;
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
  cout<<setw(10)  << key << ": ";
  cout<< val <<endl;
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
  cout<<setw(10)  << key << ": ";
  cout<< val <<endl;
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
  cout<<setw(10)  << key << ": ";
  cout<< val <<endl;
}
void paramLSRTM::print()
{
  nthread.print();
  ns.print();
  ngmax.print();
  nx.print();
  nz.print();
  nt.print();
  lpad.print();
  rpad.print();
  npml.print();
  mask.print();
  maxiter.print();
  delay.print();
  delaycal.print();
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
