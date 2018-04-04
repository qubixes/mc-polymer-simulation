#include <stdlib.h>
#include <stdio.h>

#define MAX(X,Y) ((X>Y)?X:Y)
#define MIN(X,Y) ((X<Y)?X:Y)

typedef struct PasData {
	int nPol;
	int* polId;
	int* nMono;
	int nTotMono;
	int maxNMono;
	double*** pos;
	int* centromers;
} PasData;



PasData* ReadInfoFile(char* infoFile);
PasData* PascalDataInit(char* infoFile, char* configFile);
void PrintPasData(PasData* pData);