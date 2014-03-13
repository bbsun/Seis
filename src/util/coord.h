#ifndef _COORDS_H_H
#define _COORDS_H_H
class Coords{
public:
	Coords();
	Coords(const Coords & other);
	~Coords();
	static void getSources(char * CoordFileName, int ns, float ** &sc);
	static void getReceivers(char * CoordFileName, int ns, int * ng,float ***&gc);
	static void countReceivers(char * CoordFileName, int ns, int * &ng);
	static int  countShots(char * CoordFileName);
	static void printShotInfo(int BeginShot,int EndShot,char *CoordFileName,char* OutputFileName);
};
#endif