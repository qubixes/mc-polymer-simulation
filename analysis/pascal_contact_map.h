#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rng.h"
#include "pascal_lib.h"	

typedef struct Node {
	struct Node* next;
	int iPol, iMono;
} Node;

typedef struct Contact {
	int iPol, iMono;
	int jPol, jMono;
	double strength;
} Contact;

typedef struct PCMatrix {
	Node** matrix;
	int nBins[3];
	int nTotBins;
	double dxyz;
	double xyzStart[3];
} PCMatrix;

void PrintPasData(PasData* pData);
PCMatrix* GeneratePCMatrix(PasData* pData);
Node* NewNode(int iPol, int iMono);
int Index(PCMatrix* pcm, double xyz[3]);
void PrintContactMatrix(PCMatrix* pcm, PasData* pData, char* file, unsigned int rng[4], long nSamples);
