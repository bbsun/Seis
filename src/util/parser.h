#ifndef PARSER_H_H
#define PARSER_H_H
#include <string>
class ParserFromFile
{
 public:
  ParserFromFile();
  ParserFromFile(std::string file );
  void getString(std::string key, std::string & val,std::string def,bool opt);
 private:
  std::string _file="param.txt";
};
#endif
