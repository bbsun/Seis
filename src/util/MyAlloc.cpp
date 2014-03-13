#include "MyAlloc.h"
#include <iostream>

using namespace std;

MyAlloc::MyAlloc(void)
{
}


int MyAlloc::test()
{

for (int i=0;i<1;i++){
	float *f = myAlloc<float>(10);
	for (int i=0;i<10;i++)
		f[i] = 1;
	memset(f,'\0',10);
	for (int i=0;i<10;i++)
		cout << f[i] << endl;
	//delete(f);
	myFree(f);
}
	return 0;
}

