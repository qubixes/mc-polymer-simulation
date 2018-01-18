#ifndef GPUPOL_DEF_H_INC
#define GPUPOL_DEF_H_INC
#ifndef UINT_DEF
#define UINT_DEF
typedef unsigned int uint;
#endif
#include <stdlib.h>
typedef struct uint4_t {
	uint x, y, z, w;
}uint4_t;


typedef struct Coor{
	double x,y,z;
}Coor;

typedef struct TUVCoor{
	uint t,u,v;
}TUVCoor;



typedef struct Polymer{
	uint* bonds;
	Coor* polCoor;
	TUVCoor startTUV;
	int startLeftSL;
// 	Coor* oldModes;
	uint* label;
	uint labelStart;
	int length;
}Polymer;

typedef struct SimState{
	uint* lattice, *gpuLattice, *slLattice, *labLattice;
	uint* seeds;
	int nPol;
	Polymer* pol;
	int* label2Id;
	int* id2Label;
	int* polNumTranslate;
	int labelHead;
	int headBit;
	int labelBit;
	uint4_t rngState;
}SimState;


typedef struct SimProperties{
	size_t nwg;
	size_t ws;
	
	long tMax;
	long curT;
	uint R; //Number of 32 bit values that represent state of the RNG.
	long nSteps;
	long writeInterval;
	uint nwt;
	uint nwu;
	uint nwv;
	uint LT, LU, LV, L;
	uint latSize;
	uint slSize, labSize;
	
	int equilibrated;
	int double_step;
	uint seed;
	uint gpuLatSize;
	int polLength;
	long fastEq;
	char* dirOut;
	int nDevices;
	uint* trans;
	uint maxPolLength;
	double density;
}SimProperties;
#endif
