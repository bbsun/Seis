#ifndef PARSER_H_H
#define PARSER_H_H
#include <string>
#include "../lsrtm/datastr.h"
class ParserFromFile
{
 public:
  ParserFromFile();
  ParserFromFile(std::string file );
  void   getInt(std::string key,          int& val,          int def, bool opt);
  void  getFloat(std::string key,      float & val,        float def, bool opt);
  void getString(std::string key, std::string & val, std::string def, bool opt);
  void    getInt(   paramInt & a);
  void  getFloat( paramFloat & a);
  void getString(paramString & a);
 private:
  std::string _file="param.txt";
};
#endif
