#define _XOPEN_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "timer.h"
#include "file.h"

// #define MAX_LT 256 /* Maximum size of lattice in t-direction */
// #define MAX_LU 256 /* Maxumum size of lattice in u-direction */
// #define MAX_LV 256 /* Maximum size of lattice in v-direction */
// #define MAX_LSIZE (MAX_LT*MAX_LU*MAX_LV)
#define PI (3.14159265359)
// #define NMODES 10
#define SPAC_MODES 25
#define MIN(X,Y) ((X<Y)?X:Y)
#define MAX(X,Y) ((X>Y)?X:Y)
#define POL_TYPE_LIN 1
#define POL_TYPE_RING 2
#define POL_TYPE_MIXED 3

int tuvRelPos[12][3];

int curPolId;
int LT, LU, LV;
int *t, *u, *v;
char* strIn;
fpos_t* filePos;
char** allFiles;
double *genom, *unitCor;
double rGyr;
// long THERM;
char sampleDir[10000];
char* baseDir;

typedef struct LatPos{
	int firstMono;
	int nOcc;
}LatPoint;

// typedef struct MonoList{
// 	int* next;
// }MonoList;


typedef struct Coor{
	double x,y,z;
}Coor;

typedef struct DInt{
	int x,y,ig;
}DInt;

typedef struct IDouble{
	int id; double val;
}IDouble;

typedef struct Histogram{
	double dBin;
	int nBin;
	double* avgVal;
	long* count;
	long totCount;
}Histogram;

typedef struct SimProperties{
	long maxNMono;
	int nTime, nEqd, nTherm;
	long dT, tStart;
	int nPol;
	int nDev;
	int equalLengths;
	char* resDir;
	char* sampleDir;
	char* neFile;
	char exec[100];
	double density;
	double Ne;
	double tFac;
	int doubleStep;
	int equilibrated;
	int LT, LU, LV, LSIZE;
	int updRouseStat, updRouseDyn, updGenom, updUnit, updSPRouse, updSL, updDif, updGyr;
	int updSpacDif, updShearMod, updRee, updPC, updAvgPos, updMagDip;
	int polTypeMelt;
}SimProperties;


typedef struct PolyConfig{
	int *x, *y, *z;
	int *t, *u, *v;
	double* sxmode, *symode, *szmode;
	double* cxmode, *cymode, *czmode;
	double* unitCor;
	double* genom;
	double** sinfac, **cosfac;
	int* modeList;
	int nModes;
	long nMono, polSize, polType;
	double rGyr;
	double slRat;
	Coor cms;
	double tuvCMS[3];
	double stress[3];
	double ree[3];
	double magDip[3];
}PolyConfig;

typedef struct TDT{
	int t;
	int dt;
	int idt;
}TDT;

typedef struct TDTTable{
	int nTDT;
	TDT* tdt;
	int nDt;
}TDTTable;
	

typedef struct PCMonoList{
	int next;
	int img;
	int pos;
}PCMonoList;

typedef struct PolyTimeLapse{
	PolyConfig* polys;
	double** sinfac, **cosfac;
	int* modeList;
	LatPoint* lattice;
	PCMonoList* monoList;
	int L, LIMG; /// For internal contact probability calculation.
	DInt* genomList;
	IDouble* ssMonoList;
	int nGenom;
	int nIg;
	int nGenomBin;
	int* genomIdList;
	long* genomCount;
	
	Histogram* magDipHist;
	double** genomR;
	double* sqrtList;
	int* pointDensity;
	
	int nModes;
	int nPolAdded;
	long pcBins;
	double** avgPosition;
	double* rGyrT;
	double* magDipTime;
	double** avgModesStat;
	double** avgModesDyn;
	double* cmsDif;
	double* mmDif;
	double* smDif;
	double* emDif;
	double* avgUnitCor;
	double* avgGenom;
	double* avgSL;
	double* avgShearMod;
	double* avgRee;
	double** genomProb;
	double** pc;
	double** pcss;
	double* pcAvg;
	double* pcssAvg;
	double avgRGyr;
	double avgMagDip;
	double* magDipCor;
	double**** sAvgSpacMode;
	double**** cAvgSpacMode;
	double**** avgSpacDif;
	TDTTable tTable;
	int polId, devId;
	int nTherm, nEqd;
}PolyTimeLapse;



typedef struct SpacDifPoint{
// 	double radius;
	int t,u,v;
	int nCoorInside;
}SpacDifPoint;

typedef struct LList{
	struct LList* next;
	int sdId;
}LList;

typedef struct SpacDif{
	SpacDifPoint* sDPoints;
	int nSDPoints;
	
	LList** ballList; //Given a set of tuv coordinates, which SD balls correspond to them?
// 	int nBalls[LT*LU*LV]; //Given a set of tuv coordinates, how many SD balls does it lie within?
}SpacDif;


// PolyTimeLapse ;
SpacDif sd;

double RSQ(Coor c);
void TUV2XYZ(int tuv[3], int xyz[3]);
double ShortestDistanceSQ(int coor1[3], int coor2[3]);
int CompareIDouble(const void* id1, const void* id2);
void RemPeriod(int coor[3], int coorRet[3]);
int IsValid(int a);
int TUV2Coor(int tuv[3]);
void Coor2TUV(int coor, int tuv[3]);
int AddCoor2TUV(int coor1, int tuv[3]);
void DTUV2XYZ(double tuv[3], double xyz[3]);
double DRXYZ(double xyz[3], double xyz2[3]);
Histogram* HistogramSum(Histogram* hist, int nHist);


void PTLAllocate(SimProperties* sp, PolyTimeLapse* ptl);
int GetNUpdates(SimProperties* sp, char* sampleDir);
void InitArrays(SimProperties* sp, PolyTimeLapse* ptl);
void PTLInit(SimProperties* sp, PolyTimeLapse* ptl);
void PCInit(SimProperties* sp,  PolyConfig* pc, double** sinfac, double** cosfac, int* modeList, int nModes);
void DestrArrays(SimProperties* sp, PolyTimeLapse* ptl);
void PTLDestr(SimProperties* sp, PolyTimeLapse* ptl);
void PCDestr(SimProperties* sp, PolyConfig* pc);
void InitFilePos(SimProperties* sp, int devId);
TDTTable* TTableNew(SimProperties* sp, int tFirst);
void TTableDestr(SimProperties* sp, PolyTimeLapse* ptl);
void SpacDifInit(SpacDif* sd);
void SpacDifDestr(SpacDif* sd);
void LoadPTL(SimProperties* sp, PolyTimeLapse* ptl, int polId, int devId);
void InitRelPos();
int GetNClosestNeigh(int* occLat, int* retList, int tuv[3], int volume);
double TRelaxStretched(int polSize, int polType, double nTau);
double TRelax(SimProperties* sp);
void HistogramAlloc(Histogram* hist, double dBin, int nBin);
Histogram* NewHistogram(double dBin, int nBin);
void HistogramReset(Histogram* hist);


void AddAverages(SimProperties* sp, PolyTimeLapse* ptl);
void AddSpaceRouse(SimProperties* sp, PolyTimeLapse* ptl);
void AddSpacDif(SimProperties* sp, PolyTimeLapse* ptl, SpacDif* sd);
void ComputeObsPoly(SimProperties* sp, PolyConfig* pcfg);
void ComputeGenom(SimProperties* sp, PolyConfig* pcfg);
void ComputeModes(SimProperties* sp, PolyConfig* pcfg);
void ComputeUnitCor(SimProperties* sp, PolyConfig* pcfg);
void ComputeGyration(SimProperties* sp, PolyConfig* pcfg);
void ComputeCMS(SimProperties* sp, PolyConfig* pcfg);
void ComputeSL(SimProperties* sp, PolyConfig* pcfg);
void ComputeTUVCMS(SimProperties* sp, PolyConfig* pcfg);
void ComputeShearMod(SimProperties* sp, PolyConfig* pcfg);
void ComputeRee(SimProperties* sp, PolyConfig* pcfg);
void ComputeMagDip(SimProperties* sp, PolyConfig* pcfg);

void AddRGyr(SimProperties* sp, PolyTimeLapse* ptl);
void AddSL(SimProperties* sp, PolyTimeLapse* ptl);
void AddGenom(SimProperties* sp, PolyTimeLapse* ptl);
void AddUnit(SimProperties* sp, PolyTimeLapse* ptl);
void AddModesStat(SimProperties* sp, PolyTimeLapse* ptl);
void AddModesDyn(SimProperties* sp, PolyTimeLapse* ptl);
void AddMono(SimProperties* sp, PolyTimeLapse* ptl, int mono, double* res);
void AddCMS(SimProperties* sp, PolyTimeLapse* ptl);
void AddMM(SimProperties* sp, PolyTimeLapse* ptl);
void AddSM(SimProperties* sp, PolyTimeLapse* ptl);
void AddEM(SimProperties* sp, PolyTimeLapse* ptl);
void AddSpaceRouse(SimProperties* sp, PolyTimeLapse* ptl);
void AddShearMod(SimProperties* sp, PolyTimeLapse* ptl);
void AddRee(SimProperties* sp, PolyTimeLapse* ptl);
void AddContactProbability(SimProperties* sp, PolyTimeLapse* ptl);
void AddGenomNew(SimProperties* sp, PolyTimeLapse* ptl);
void AddGenomNew2(SimProperties* sp, PolyTimeLapse* ptl);
void AddAvgPos(SimProperties* sp, PolyTimeLapse* ptl);
void AddMagDip(SimProperties* sp, PolyTimeLapse* ptl);

void SetSimProps(SimProperties* sp, char* sampleDir);
void WriteGyration(SimProperties* sp, PolyTimeLapse* ptl);
void WriteGenom(SimProperties* sp, PolyTimeLapse* ptl);
void WriteUnitCor(SimProperties* sp, PolyTimeLapse* ptl);
void WriteAllFiles(SimProperties* sp, PolyTimeLapse *ptl);
void WriteModesStat(SimProperties* sp, PolyTimeLapse* ptl);
void WriteModesDyn(SimProperties* sp, PolyTimeLapse* ptl);
void WriteCMSDiff(SimProperties* sp, PolyTimeLapse* ptl);
void WriteSL(SimProperties* sp, PolyTimeLapse* ptl);
void WriteDiff(SimProperties* sp, PolyTimeLapse* ptl, double* dif, char* base, long nEqd);
void WriteSpaceRouse(SimProperties* sp, PolyTimeLapse* ptl);
void WriteSpacDif(SimProperties* sp, PolyTimeLapse* ptl);
void WriteShearMod(SimProperties* sp, PolyTimeLapse* ptl);
void WriteRee(SimProperties* sp, PolyTimeLapse* ptl);
void WriteMagDip(SimProperties* sp, PolyTimeLapse* ptl);
void WriteContactProbability(SimProperties* sp, PolyTimeLapse* ptl);
void WriteAvgContactProbability(SimProperties* sp, PolyTimeLapse* ptl);
int CharToHex(char c);
void WriteAvgPos(SimProperties* sp, PolyTimeLapse* ptl);
void LoadPol(SimProperties* sp, PolyTimeLapse* pt);
double** ReadAvgPos(SimProperties* sp, PolyTimeLapse* ptl);

