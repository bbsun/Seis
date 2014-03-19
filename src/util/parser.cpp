#include "parser.h"
#include <fstream>
#include <iostream>
using std::string;
using std::npos;
ParserFromFile::ParserFromFile(string file )
{
  this->_file = file
}

void ParserFromFile::getString(string key, string & val,string def,bool opt);
{
  bool findKey=false;
  bool findValue=false;
    const string cmt="!";
    string file=_file;
    ifstream in(file);
    string line;
    while(getline(in,line))
      {
	line.find(cmt)?line=line.erase(line.find(cmt)):;
	  if(line.find(key)!=npos)
	    {
	      findKey = true;
	      line = line.substr(line.find(key)+key.size());
	      string tmp=0;
		if((auto index=line.find_first_of("+-.0123456789"))!=npos)
		  {
		    findValue = true;
		    while(!isspace(line[index]) && index <line.size()){
		      tmp+=line[index];
		      index++;
		    }
		  }
		val = tmp;
		return;
	    }
      }
    // set defualt value
      
  }
