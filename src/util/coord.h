#ifndef _COORDS_H_H
#define _COORDS_H_H
#include <string>
class Coords{
public:
	Coords();
	Coords(const Coords & other);
	~Coords();
	static void getSources(std::string  CoordFileName, int ns, float ** &sc);
	static void getReceivers(std::string  CoordFileName, int ns, int * ng,float ***&gc);
	static void countReceivers(std::string  CoordFileName, int ns, int * &ng);
	static int  countShots(std::string  CoordFileName);
	static void printShotInfo(int BeginShot, int EndShot, int ns, int ngmax, int * ng, float **sc, float ***gc);
	static void printShotInfo(int BeginShot,int EndShot,std::string CoordFileName,std::string CoordFileNameOut);
};
#endif
