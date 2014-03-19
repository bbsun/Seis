#include "debug.h"
#include <iostream>
#include <stdlib.h>
using std::string;
void check(bool eval, char * str)
{
  if(!eval)
    {
      std::cout<<"ERROR: "<< str << std::endl;
      exit(0);
    }
}
void error(string errinfo)
{
printf("\033[1;4;%dm%s\033[0m",31,errinfo);
}
