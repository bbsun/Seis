#include <stdio.h>
#include "add.h"
#include <iostream>
#include "hello.h"
#include "file2.h"
#include <string>
#include <fstream>
#include "../util/memory.h"
#include "../util/arraymath.h"
#include "../util/coord.h"
using std::string;
#define print std::cout 
#define ed    std::endl;
using std::ifstream;
float * aaa;
extern int ix;
int ix = 1000;
int a_out;
// variable defined outside function has a default value of zero
void cppstudy(int argc, char* argv[],char * envp[])
{
  string s2 = "pi = 3.1415926 ! parameters for the pi";
  const string arg = "pi";
  const string cmt = "!";
  std::cout<<s2<<std::endl;
  if(s2.find(cmt))
    s2 =s2.erase(s2.find(cmt));
  std::cout<<s2<<std::endl;
  string s3=s2.substr(s2.find(arg)+arg.size());
  string tmpxx;
  std::cout<<s3<<std::endl;
  auto indexX = s3.find_first_of("+-.0123456789");
  while(! isspace(s3[indexX]) && indexX <s2.size()){
    tmpxx+=s3[indexX];
    indexX++;
  }
  float tmpf = stof(tmpxx);
  std::cout<<tmpf<<std::endl;
  const string s="keep out";
  for(auto &c:s)
    std::cout<<c<<std::endl;
  print<<"argc "<<argc<<ed;
  print<<"argv"<<argv[0]<<ed;
  // string study
  string tmpX;

  getline(std::cin,tmpX);
  decltype(tmpX.size()) index=0;
  while(index < tmpX.size())
    {
       if(isalnum(tmpX[index]))
	tmpX[index]='X';
      index++;
    }
  std::cout<<tmpX<<std::endl;
  for(char &c:tmpX)
    {
      if(isalnum(c))
	c='X';
    }
  for(auto &c:tmpX)
    {
      if(isalnum(c))
	c='X';
    }
  std::cout<<tmpX<<std::endl;
  while(getline(std::cin,tmpX))
    std::cout<<tmpX<<" "<<tmpX.size()<<std::endl;
  
  
  // std::cout<<tmpX<<" its size is "<<tmpX.size()<<std::endl;
  std::cout<<"finish"<<std::endl;
  if(false){
  string s1,s2;
  std::cin >>s1>>s2;
  std::cout<< s1 <<s2<<endl;
  }
  string testfile("/home/sunbb/sh/param.txt");
  ifstream in(testfile);
  string tmp;
  std::cout<<"read again and again"<<std::endl;
  while(in>>tmp)
    std::cout<<tmp<<std::endl;
  std::cout<<"read each line "<<std::endl;
  in.close();
  in.open(argv[1]);
  while(getline(in,tmp))
    std::cout<<tmp<<std::endl;
  exit(0);
}
int main(int argc, char *argv[], char *envp[])
{
  if(true)
    cppstudy(argc,argv,envp);
  char coordsFile[1024] = "/home/sunbb/sh/coords_ns256.txt";
  char outFile[1024]="/home/sunbb/sh/coords_ns256_out.txt";
  Coords::printShotInfo(1,10,coordsFile,outFile);
  Coords::printShotInfo(1,10,coordsFile,0);
  return 0;
  float * f=MyAlloc<float>::alc(10);
  float * d=MyAlloc<float>::alc(10);
  float * k=MyAlloc<float>::alc(10);
  for(int i=0;i<10;i++){
    d[i]=i*1.0f;
    k[i]=11.0f-i*1.0f;
  }
  dump("d",d,10);
  dump("k",k,10);
  print<<"out = d.*d"<<ed;
  double a = opern(d,d,DOTPRODUCT,10);
  cout<<a<<endl;
  print<<"d=10"<<ed;
  opern(d,VALUE,10,10.0f); dump("d",d,10);
  print<<"d=random"<<ed;
  opern(d,RANDOM,10);        dump("d",d,10);
  print<<"f = d"<<ed;
  opern(f,d,(float *)0,COPY,10); dump("f",f,10);
  print<<"f = sqrt(d)"<<ed
  opern(f,d, d, SQRT,10);        dump("f",f,10);
  print<<"f = d*10"<<ed;
  opern(f,d, d, SCAL,10,10.0);   dump("f",f,10);
  print<<"f = d + k"<<ed;
  opern(f,d, k, ADD,10);   dump("f",f,10);
  print<<"f = d - k"<<ed;
  opern(f,d, k, SUB,10);   dump("f",f,10);
  print<<"f = d * k"<<ed;
  opern(f,d, k, MUL,10);   dump("f",f,10);
  print<<"f = d / k"<<ed;
  opern(f,d, k, DIV,10);   dump("f",f,10);dump("f",f,10);
  print<<"f = 2*d + 3*k"<<ed;
  opern(f,d, k, ADDSCAL,10,2.0,3.0);   dump("f",f,10);

  float ** f2=MyAlloc<float>::alc(3,3);
  float ** d2=MyAlloc<float>::alc(3,3);
  float ** k2=MyAlloc<float>::alc(3,3);
  for(int j=0;j<3;j++)
  for(int i=0;i<3;i++){
    int ii=i+1;
    int jj=j+1;
    d2[i][j]=jj+ii*1.0f;
    k2[i][j]=ii*jj*1.0f;
  }
  dump("d2",d2,3,3);
  dump("k2",k2,3,3);
  print<<"f2 = d2"<<ed;
  opern(f2,d2, d2,COPY,3,3); dump("f2",f2,3,3);
  print<<"f2 = sqrt(d2)"<<ed
    opern(f2,d2, d2, SQRT,3,3);        dump("f2",f2,3,3);
  print<<"f2 = d2*10"<<ed;
  opern(f2,d2, d2, SCAL,3,3,10.0);   dump("f2",f2,3,3);
  print<<"f2 = d2 + k2"<<ed;
  opern(f2,d2, k2, ADD,3,3);   dump("f2",f2,3,3);
  print<<"f2 = d2 - k2"<<ed;
  opern(f2,d2, k2, SUB,3,3);   dump("f2",f2,3,3);
  print<<"f2 = d2 * k2"<<ed;
  opern(f2,d2, k2, MUL,3,3);   dump("f2",f2,3,3);
  print<<"f2 = d2 / k2"<<ed;
  opern(f2,d2, k2, DIV,3,3);   dump("f2",f2,3,3);
  print<<"f2= 2*d2 + 3*k2"<<ed;
  opern(f2,d2, k2, ADDSCAL,3,3,2.0,3.0);   dump("f2",f2,3,3);
  return 0;
  MyAlloc<float>::free(f);
  double ***ff=MyAlloc<double>::alc(10,10,10);
  MyAlloc<double>::free(ff);
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
