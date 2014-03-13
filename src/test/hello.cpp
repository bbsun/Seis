#include <stdio.h>
#include "add.h"
#include <iostream>
#include "hello.h"
#include "file2.h"
#include <string>
#include "../util/memory.h"
using std::string;
#define print std::cout 
#define ed    std::endl;
float * aaa;
extern int ix;
int ix = 1000;
int a_out;
// variable defined outside function has a default value of zero
int main(int argc, char *argv[], char *envp[])
{
  float * f=MyAlloc<float>::alc(10);
  MyAlloc<float>::free(f);
  const string stringxx("const string");
  print << stringxx <<ed;
  for(auto & c:stringxx)
    {
      print<<c;
    }
  string testpunc(" ni haoma,wo jiao sun bingbing,.. asdasd!!!");
  print<<testpunc<<ed;
  for(auto &c: testpunc)
    {
      c = ispunct(c)?' ':c;
    }
  print<<testpunc<<ed;
  string s;
  print << s[2] <<ed;
  string stringa("abscdef");
  string stringb = stringa;
  string stringc = stringa;
  print<<stringa<<ed;
  for(char &c:stringa)
    c='X';
  print<<stringa<<ed;
  for(decltype(stringb.size()) index=0;index<stringb.size();++index)
    stringb[index] = 'X';
  print<<stringb<<ed;
  decltype(stringc.size()) index=0;
  while(stringc[index])
    {
	  stringc[index] = 'X';
	  index++;
    }
  print<<stringc<<ed;
  return 0;
  string teststring("asdasdasddsadadasd  asdasd asdasd asds d ");
  for(decltype(teststring.size()) index=0;
      index != teststring.size() && !isspace(teststring[index]); ++ index)
    teststring[index] = toupper(teststring[index]);
  print << teststring<<ed;
  for(auto c : teststring)
    print<<c<<ed;
  for(auto &c :teststring)
    {
      c = toupper(c);
    }
    for(auto c : teststring)
    print<<c<<ed;
  string big;
  string word;
  string line;
  while(std::cin >> word)
    {
      big +=" "+word;
      print<<"now big is "<<big<<ed;
    }
  while( getline(std::cin,line))
    print<< line   << ed;
  while(std::cin >> word)
    print << word << ed;
  string sa("bbsunok");
  print<<"size of "<<sa <<sa.size()<<ed;
  return 0;
  print<<"You have inputed total "<< argc <<ed;
  
  for( int i=0; i<argc ; i++)
    {
      printf("arg%d : %s\n",i,argv[i]);
    }
  
  printf("The follow is envp : \n");
  for(int i=0; envp[i];i++)
   {
      printf("%s \n", envp[i]);
   }
  return 0;
  double dval = 3.14;
  double *fd = &dval;
  print<<*fd<<ed;
  (*fd)++;
  print<<*fd<<ed;
  print<<aaa<<" "<<ed;
  int vv=1;
  print<<vv<<" "<<ed;
  int & vvf1 = vv;
  print<<vvf1<<" "<<ed;
  int & vvf2 = vvf1;
  print<<vvf2<<" "<<ed;
  int ii=100,sum=0;
  for(int ii = 0; ii != 10; ++ ii)
    sum +=ii;
  print<<ii<<"  "<<sum<<ed;
  int _=1;
  print<<_<<ed;
  fun2();
  return 0;
  int a_in;
  // variable defined inside has no defualt value.
  print<<" a_out: "<<a_out<<ed;
  print<<" a_in:  "<<a_in <<ed;
  print<<" ix     "<<ix   <<ed;
  int input_value;
  int ix = 3.14; 
  print<<ix<<ed;
  int jx = {3.14};
  print<<jx<<ed;
  std::cin>>input_value;
  double wage{0};
  double salary = wage=999.99;
  std::cout<<wage<<std::endl;
  std::cout<<salary<<std::endl;
  long double ld = 3.1415926;
  int inta{ld};
  int intb = {ld};
  std::cout<<ld<<std::endl;
  std::cout<<inta<<std::endl;
  std::cout<<intb<<std::endl;
  std::cout<< "a very very long long "
    "long long string "
    "long loog string "<<std::endl;
  unsigned u1=20, u2 = 10;
  std::cout<< u1 - u2 <<std::endl;
  std::cout<< u2 - u1 <<std::endl;
  unsigned u = 10;
  int i = -42;
  std::cout<< i + i << std::endl;
  std::cout<< u + i << std::endl;
  char b = 'c';
  std::cout<<"size of char is "<< sizeof(char)<<std::endl;
  auto axx=10;
  std::cout<<axx<<std::endl;
  printf("hello world");
  std::cout << "Enter two numbers:"<<std::endl;
  float v1 = 0, v2 = 0;
  //std::cin>>v1>>v2;
  std::cout<<"The sum of "<< v1 << " and "<< v2 
	   <<" is " << v1 + v2 << std::endl;
  return 0;
}
