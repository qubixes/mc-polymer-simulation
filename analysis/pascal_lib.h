#include <stdlib.h>
#include <stdio.h>
#include "lowm_modes.h"

#define MAX(X,Y) ((X>Y)?X:Y)
#define MIN(X,Y) ((X<Y)?X:Y)

#define BOUNDARY_PERIOD 0
#define BOUNDARY_STATIC 1


typedef struct PasData {
	int nPol;
	int* polId;
	int* nMono;
	int nTotMono;
	int maxNMono;
	double*** pos;
	int* centromers;
} PasData;

typedef struct Data{
	int*** tuv;
	double*** xyz;
	int* nMono;
	int* polTypes;
	double (*TUV2Distance) (double*, double*, int);
	int nPol;
	int maxNMono;
	int nTotMono;
	int L;
}Data;

/// Histogram that starts at 0;
// typedef struct Histogram{
// 	int nBins;
// 	double dBin;
// 	int* counts;
// 	int nCounts;
// }Histogram;

Data* NewData(int maxNMono, int nPol, int L, int boundaryCond);
PasData* ReadInfoFile(char* infoFile);
PasData* PascalDataInit(char* infoFile, char* configFile);
void PrintPasData(PasData* pData);
Data* ReadData(char* file, int boundaryCond);
double TUV2DistancePeriodic(double tuv1[3], double tuv2[3], int L);
double TUV2DistanceStatic(double tuv1[3], double tuv2[3], int L);
void AddCharToTUV(char c, int tuv[3]);
void DblTUV2XYZ(double tuv[3], double xyz[3]);
double Squabs(double xyz[3]);
