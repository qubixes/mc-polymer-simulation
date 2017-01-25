#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "rng.h"
#include "timer.h"

typedef unsigned int uint;

#define TRUE 1
#define FALSE 0
#define LAT_POW2 TRUE
#define POL_RING 1
#define POL_LIN 2

#ifndef POL_TYPE 
	#define POL_TYPE POL_LIN
#endif

#define MIN(X,Y) ((X<Y)?X:Y)
#define MAX(X,Y) ((X>Y)?X:Y)

// #define BLT 5 // 2^5=32
// #define BLU 5 // 2^6=64
// #define BLV 5 // 2^7=128

// #define LT (1<<BLT)
// #define LU (1<<BLU)
// #define LV (1<<BLV)

// #define ALLOC_T (2*LT)
// #define ALLOC_U (2*LU)
// #define ALLOC_V (LV)
#define NOPT 16

// #define LAT_SIZE (LT*LU*LV)
// #define LAT_ALLOC ((ALLOC_T*ALLOC_V*ALLOC_U)/32)

// #define LARGE_T_UNIT 1
// #define LARGE_U_UNIT ALLOC_T
// #define LARGE_V_UNIT (ALLOC_T*ALLOC_U)
// #define LARGE_W_UNIT (ALLOC_T*ALLOC_U*ALLOC_V*2)
// #define LARGE_TUV (LARGE_T_UNIT|LARGE_U_UNIT|LARGE_V_UNIT|LARGE_W_UNIT)

// #define ADD_BASE (LT+LU*ALLOC_T+LV*ALLOC_T*ALLOC_U)
// #define T_MASK (ALLOC_T-1)
// #define U_MASK ((ALLOC_U-1)*ALLOC_T)
// #define V_MASK ((ALLOC_V-1)*ALLOC_T*ALLOC_U)

// #define BOUND_MASK 0x007f7f7f
// #define BOUND_MASK ((LT-1)|((LU-1)*ALLOC_T)|((LV-1)*ALLOC_T*ALLOC_U))
// #define TUV_MASK 0x00ffffff
#define BLOCK_SELECT 8

typedef struct Constants{
	uint BL;
	uint L;
	uint ALLOC_T, ALLOC_U, ALLOC_V;
	uint LAT_SIZE;
	uint LAT_ALLOC;
	uint LARGE_T_UNIT, LARGE_U_UNIT, LARGE_V_UNIT, LARGE_W_UNIT;
	uint LARGE_TUVW;
	uint ADD_BASE;
	uint T_MASK, U_MASK, V_MASK;
	uint BOUND_MASK;
}Constants;

typedef struct LookupTables{
	unsigned int transTable[16*16*16];
	int curNTrans[16*16];
	int addOccup[16*16*16];
	int remOccup[16*16*16];
	double accRat[16*16*16];
}LookupTables;

typedef struct CurState{
	uint* lattice;
	
	uint* unit1Cache;
	uint* unit2Cache;
	uint* coorCache;
	
	uint* unitPol;
	uint* coorPol;
	
	int intsPerPol;
	int intsPerMono;
	
	int nPol;
	int polSize;
	
	int allocated;
	Constants con;
}CurState;

typedef struct Config{
	long tMax;
	long* time;
	int* nSample;
	int* sampled;
	int nTime;
	long totT;
	long tMaxSS;
	long interval;
	long eqInterval;
	long nEq;
	int nRun;
	unsigned int seed;
	char* dir;
	double density;
	double EField;
	int efPolId;
	uint rngState[4];
}Config;

typedef struct CoorI{
	int x[3];
}CoorI;

typedef struct CoorD{
	double x[3];
}CoorD;

typedef struct PolyConfig{
	CoorI* tuv;
	CoorI* xyz;
	int polSize;
// 	CoorD cms;
}PolyConfig;


typedef struct Result{
	PolyConfig pcfg;
	double* gen;
	int* genCounts;
}Result;


Config cfg;
CurState cs[10];
LookupTables tab;
Result res[10];

uint TCoor(unsigned int coor, Constants* con);
uint UCoor(unsigned int coor, Constants* con);
uint VCoor(unsigned int coor, Constants* con);
int X(int t, int u, int v);
int Y(int t, int u, int v);
int Z(int t, int u, int v);
int XUnit(int unit);
int YUnit(int unit);
int ZUnit(int unit);

int TUVtoCoor(int t, int u, int v, Constants* con);
uint SmallToLargeUnit(uint uSmall, Constants* con);
uint LargeToSmallUnit(uint uLarge, Constants* con);
int CharToHex(char c);
uint IsValid(uint a);
uint ValidateAddUnitVectors(uint a, uint b, uint* c);
uint AddUnitVectors(uint a, uint b);
uint AddUnitToCoor(uint addUnit, uint coor, Constants* con);
double TRelax(int polSize);
int GetNPol(double density, int polSize, int LAT_SIZE);
void SetLattice(int coor, CurState* cs);
void UnsetLattice(int coor, CurState* cs);
uint OccupLattice(int coor, CurState* cs);
void CopyState(CurState* srcCS, CurState* dstCS);
void ConfigInit();
void PCInit(CurState* cs, PolyConfig* pcfg);
int DoMCStep(CurState* cs, int nTime);
int DoStep(CurState* cs);

void SetTransTable(Constants* con);
void CreatePolymers(Config* cfg, CurState* cs);
void ResultInit(Config* cfg, CurState* cs, Result* res);
void SetConstants(Constants* con, int BL);
void CSInit(CurState* cs, int BL, int polSize, int nPol);
void DoubleLattice(CurState* cs, CurState* cs_next);
void CSClear(CurState* cs);

void WritePolymer(CurState* cs, int iPol, FILE* pFile);
int WriteLatticeFile(CurState* cs, char* file);
int CheckIntegrity(CurState* cs);
int ReadLatticeFile(CurState* cs, char* file);
int GetDtMax(char* dir, long* dt, long* tMax);
void WriteConfig(Config* cfg);
void PrintTrans();
void WriteGenom(Config* cfg, CurState* cs, Result* res);
void WriteResults(Config* cfg, CurState* cs, Result* res);
void PrintCoor(int coor, Constants* con);
void PrintPol(int iPol, CurState* cs);

void LoadPolymer(PolyConfig* pcfg, CurState* cs, int iPol);
void MeasVars(CurState* cs, Result* res);
void AddGenom(Result* res);

