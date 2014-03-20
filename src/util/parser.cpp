#include "parser.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "debug.h"
using std::ifstream;
using std::string;
ParserFromFile::ParserFromFile()
{
  std::cout<<"default input file :"<< _file <<std::endl;
}
ParserFromFile::ParserFromFile(string file )
{
  this->_file = file;
}
void ParserFromFile::getFloat(string key, float & val,float def,bool opt)
{
  bool findKey   = false;
  const string cmt="!";
  string file = _file;
  ifstream in(file);
  if(!in) {
    error("File can not open : ");
    std::cout << file << std::endl;
    exit(0);
  }
  string line;
  while(getline(in,line)){
    line.find(cmt)!=string::npos?line=line.erase(line.find(cmt)):line;
    if(line.find(key)!=string::npos){
      findKey = true;
      line = line.substr(line.find(key)+key.size());
      if(line.find_first_of("+-.0123456789") != string::npos){
	auto index = line.find_first_of("+-.0123456789");
	string tmp;
	while(!isspace(line[index]) && index <line.size()){
	  tmp += line[index];
	  index++;
	}
	val = stof(tmp);  // process the value (transformation may needed)
	in.close();
	return ;
      }
    }
  }
  in.close();
  if(opt){
    val = def;
    return;
  }
  if(!findKey){
    error("error: ");
    std::cout<<"key "<< key << " not found"<< std::endl;
    exit(0);
    return;
  }
  error("error: ");
  std::cout<<"key "<< key << " value missing"<< std::endl;
  exit(0);
  return;
  }
void ParserFromFile::getInt(string key, int & val,int def,bool opt)
{
  bool findKey   = false;
  const string cmt="!";
  string file = _file;
  ifstream in(file);
  if(!in) {
    error("File can not open : ");
    std::cout << file << std::endl;
    exit(0);
  }
  string line;
  while(getline(in,line)){
    line.find(cmt)!=string::npos?line=line.erase(line.find(cmt)):line;
    if(line.find(key)!=string::npos){
      findKey = true;
      line = line.substr(line.find(key)+key.size());
      if(line.find_first_of("+-0123456789") != string::npos){
	auto index = line.find_first_of("+-0123456789");
	string tmp;
	while(!isspace(line[index]) && index <line.size()){
	  tmp += line[index];
	  index++;
	}
	val = stoi(tmp);  // process the value (transformation may needed)
	in.close();
	return ;
      }
    }
  }
  in.close();
  if(opt){
    val = def;
    return;
  }
  if(!findKey){
    error("error: ");
    std::cout<<"key "<< key << " not found"<< std::endl;
    exit(0);
    return;
  }
  error("error: ");
  std::cout<<"key "<< key << " value missing"<< std::endl;
  exit(0);
  return;
  }
void ParserFromFile::getString(string key, string & val,string def,bool opt)
{
  bool findKey   = false;
  const string cmt="!";
  string file = _file;
  ifstream in(file);
  if(!in) {
    error("File can not open : ");
    std::cout << file << std::endl;
    exit(0);
  }
  string line;
  while(getline(in,line)){
    line.find(cmt)!=string::npos?line=line.erase(line.find(cmt)):line;
    if(line.find(key)!=string::npos){
      findKey = true;
      line = line.substr(line.find(key)+key.size());
      if(line.find_first_of("/qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM") != string::npos){
	auto index = line.find_first_of("/qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM");
	string tmp;
	while(!isspace(line[index]) && index <line.size()){
	  tmp += line[index];
	  index++;
	}
	val = tmp;  // process the value (transformation may needed)
	in.close();
	return ;
      }
    }
  }
  in.close();
  if(opt){
    val = def;
    return;
  }
  if(!findKey){
    error("error: ");
    std::cout<<"key "<< key << " not found"<< std::endl;
    exit(0);
    return;
  }
  error("error: ");
  std::cout<<"key "<< key << " value missing"<< std::endl;
  exit(0);
  return;
  }
void ParserFromFile::getInt(paramInt &a)
{
  getInt(a.key,a.val,a.def,a.opt);
}
void ParserFromFile::getFloat(paramFloat & a)
{
  getFloat(a.key,a.val,a.def,a.opt);
}
void ParserFromFile::getString(paramString & a)
{
  getString(a.key,a.val,a.def,a.opt);
}
