#include "debug.h"
#include <iostream>
#include <stdlib.h>
void check(bool eval, char * str)
{
  if(!eval)
    {
      std::cout<<"ERROR: "<< str << std::endl;
      exit(0);
    }
}
