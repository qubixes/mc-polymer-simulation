#include <stdio.h>
#include <stdlib.h>
#define IS_MAIN

#include "rng.h"
#include "denspol_lib.h"
#define LAT_SHAPE_SPHERE 0
#define LAT_SHAPE_EMPTY 1


typedef struct BoxState{
	int nPol;
	int maxPolSize;
	int L;
	
	int** startTUV;
	int** bonds;
	int** coor;
	int* nMono;
	int* polTypes;
	int* topo;
	int* bondOcc;
	unsigned int seed[4];
}BoxState;

void PrintCoor(int coor, int L);
int TUV2CoorX(int t, int u, int v, int L);
int AddUnitToCoor(int unit, int pos, int L);
int UpscaleCoor(int coor, int L);
BoxState* NewBoxState(int L, int nPol, int maxPolSize);
BoxState* LoadSystem(char* file, char* topoFile);
BoxState* UpscaleBox(BoxState* bs, char* lengthFile, int* topoStraight, int latShape, double density);
int* LoadStraightTopo(char* file);
void WriteState(BoxState* bs, char* file, char* topoFile);
void UpdateWithSphereBoundary(BoxState* box, char* file, double density);


int main(int argc, char** argv){
	char* file, *fileOut, *topoFile, *topoFileOut, *strTopoFile, *lengthFile=NULL;
	int latShape=LAT_SHAPE_EMPTY;
	double density=0.0;
	BoxState* bs, *newBs;
	int* topoStraight;
	if(argc<6){
		printf("Not enough arguments\n");
		return 192;
	}
	
	file        = argv[1];
	topoFile    = argv[2];
	fileOut     = argv[3];
	topoFileOut = argv[4];
	strTopoFile = argv[5];
	if(argc>=8){
		if(!strcmp(argv[6], "sphere"))
			latShape = LAT_SHAPE_SPHERE;
		lengthFile = argv[7];
		density = atof(argv[8]);
	}
	
	topoStraight = LoadStraightTopo(strTopoFile);
	bs = LoadSystem(file, topoFile);
	newBs = UpscaleBox(bs, lengthFile, topoStraight, latShape, density);
	WriteState(newBs, fileOut, topoFileOut);
	
	return 0;
}

int TUV2CoorX(int t, int u, int v, int L){
	return t+L*u+L*L*v;
}

int AddUnitToCoor(int unit, int pos, int L){
	int t = pos%L;
	int u = (pos/L)%L;
	int v = pos/(L*L);
	
	int dw =  (unit>>3)&0x1;
	int dt =  (unit&0x1)     - dw;
	int du = ((unit>>1)&0x1) - dw;
	int dv = ((unit>>2)&0x1) - dw;
	
	t += dt+L; u += du+L; v += dv+L;
	
	return TUV2CoorX(t%L,u%L,v%L,L);
}

int UpscaleCoor(int pos, int L){
	int t = pos%L;
	int u = (pos/L)%L;
	int v = pos/(L*L);
	
	t *= 2; u *= 2; v *= 2;
	return TUV2CoorX(t,u,v,2*L);
}

void ResizeBoxState(BoxState* box, int* nMono){
	int maxMono=0;
	for(int iPol=0; iPol<box->nPol; iPol++){
		box->nMono[iPol] = nMono[iPol];
		box->bonds[iPol] = realloc(box->bonds[iPol], sizeof(int)*nMono[iPol]);
		box->coor[iPol]  = realloc(box->coor[iPol], sizeof(int)*nMono[iPol]);
		maxMono = MAX(nMono[iPol], maxMono);
	}
	box->maxPolSize = maxMono;
}

BoxState* NewBoxState(int L, int nPol, int maxPolSize){
	BoxState* bs = malloc(sizeof(BoxState));
	bs->L          = L;
	bs->maxPolSize = maxPolSize;
	bs->nPol       = nPol;
	int LSize = L*L*L;
	
	bs->bonds   = malloc(sizeof(int*)*nPol);
	bs->coor    = malloc(sizeof(int*)*nPol);
	bs->startTUV= malloc(sizeof(int*)*nPol);
	bs->nMono   = malloc(sizeof(int) *nPol);
	bs->polTypes= malloc(sizeof(int) *nPol);
	bs->topo    = malloc(sizeof(int) *LSize);
	bs->bondOcc = malloc(sizeof(int) *LSize);
	
	
	for(int iPol=0; iPol<nPol; iPol++){
		bs->bonds[iPol] = malloc(sizeof(int)*maxPolSize);
		bs->coor[iPol] = malloc(sizeof(int)*maxPolSize);
		bs->startTUV[iPol] = malloc(sizeof(int)*3);
	}
	
	
	for(int i=0; i<LSize; i++){
		bs->bondOcc[i] = 0;
		bs->topo[i] = 0;
	}
	return bs;
}

BoxState* LoadSystem(char* file, char* topoFile){
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	int LT, maxPolSize, nPol;
	
	
	fscanf(pFile, "%*s %i", &LT);
	fscanf(pFile, "%*s %i", &LT);
	fscanf(pFile, "%*s %i", &LT);
	fscanf(pFile, "%*s %i", &nPol);
	fscanf(pFile, "%*s %i", &maxPolSize);
	
	BoxState* bs = NewBoxState(LT, nPol, maxPolSize);
	char* str = malloc(sizeof(char)*(maxPolSize+1));
// 	int t,u,v;
	
	for(int iPol=0; iPol<nPol; iPol++){
		fscanf(pFile, "%*s %i", &bs->nMono[iPol]);
		fscanf(pFile, "%i %i %i", &bs->startTUV[iPol][0], &bs->startTUV[iPol][1], &bs->startTUV[iPol][2]);
		fscanf(pFile, "%s", str);
		
		int coor = TUV2CoorX(bs->startTUV[iPol][0], bs->startTUV[iPol][1], bs->startTUV[iPol][2], LT);
		bs->coor[iPol][0] = coor;
		for(int iBond=0; iBond<bs->nMono[iPol]; iBond++){
			bs->bonds[iPol][iBond] = CharToHex(str[iBond]);
			coor = AddUnitToCoor(bs->bonds[iPol][iBond], coor, bs->L);
			if(iBond+1<bs->nMono[iPol])
				bs->coor[iPol][iBond+1] = coor;
		}
		if(str[bs->nMono[iPol]-1] == 'f') bs->polTypes[iPol] = POL_TYPE_LIN;
		else bs->polTypes[iPol] = POL_TYPE_RING;
	}
	
	fclose(pFile);
	
	pFile = fopen(topoFile, "r");
	if(!pFile){
		printf("Error opening file %s\n", topoFile);
		exit(192);
	}
	
	for(int i=0; i<4; i++){
		fscanf(pFile, "%u", bs->seed+i);
	}
	
	for(int i=0; i<LT*LT*LT; i++){
		fscanf(pFile, "%i %i", bs->topo+i, bs->bondOcc+i);
	}
	
	fclose(pFile);
	return bs;
}

BoxState* UpscaleBox(BoxState* bs, char* lengthFile, int* topoStraight, int latShape, double density){
	BoxState* newBs = NewBoxState(bs->L*2, bs->nPol, bs->maxPolSize*8);
	int nLatticeUsed, nBondUsed;
	if(latShape == LAT_SHAPE_SPHERE)
		nLatticeUsed = SetLatticeSphere(newBs->topo, newBs->bondOcc, newBs->L, &nBondUsed);
	else
		nLatticeUsed = SetLatticeEmpty(newBs->topo, newBs->bondOcc, newBs->L, &nBondUsed);
	
	int* newNMono;
	if(lengthFile)
		newNMono = ComputePolLengths(lengthFile, nBondUsed, density);
	else{
		newNMono = malloc(sizeof(int)*newBs->nPol);
		for(int iPol=0; iPol<newBs->nPol; iPol++)
			newNMono[iPol] = 8*bs->nMono[iPol];
	}
	ResizeBoxState(newBs, newNMono);

	int maxNMono=0;
	for(int iPol=0; iPol<bs->nPol; iPol++)
		maxNMono = MAX(maxNMono, newNMono[iPol]);
	
	for(int iPol=0; iPol<bs->nPol; iPol++){
		double avgSlAdded;
		if(bs->polTypes[iPol] == POL_TYPE_LIN)
			avgSlAdded = (newNMono[iPol]-(2*bs->nMono[iPol]-1))/(double)(2*bs->nMono[iPol]-1);
		else
			avgSlAdded = (newNMono[iPol]-2*bs->nMono[iPol])/(double)(2*bs->nMono[iPol]);
		
		avgSlAdded = MAX(avgSlAdded, 0);
		
		for(int i=0; i<3; i++)
			newBs->startTUV[iPol][i] = bs->startTUV[iPol][i]*2;
		
		int curCoor= TUV2CoorX(newBs->startTUV[iPol][0], newBs->startTUV[iPol][1], newBs->startTUV[iPol][2], newBs->L);
		
		double leftover=0;
		int jBondNew=0;
		for(int iBond=0; iBond<bs->nMono[iPol]; iBond++){
			int bond = bs->bonds[iPol][iBond];
			if(bond == 0xf){
				for(int i=0; i<(int)(leftover+avgSlAdded+0.5); i++)
					newBs->bonds[iPol][jBondNew++] = 0x0;
				newBs->bonds[iPol][jBondNew++] = 0xf;
				break;
			}
			
			double dSl = leftover+avgSlAdded;
			int nSl = (int)(dSl+0.5);
			leftover = dSl-nSl;
			for(int i=0; i<nSl; i++)
				newBs->bonds[iPol][jBondNew++] = 0x0;
			
			newBs->bonds[iPol][jBondNew++] = bond;
			
			if(bond != 0x0){
				curCoor = AddUnitToCoor(bond, curCoor, newBs->L);
				newBs->topo[curCoor]    = topoStraight[bond];
				newBs->bondOcc[curCoor] = (1<<bond)|(1<<(bond^0xf));
			}
			
			dSl = leftover+avgSlAdded;
			nSl = (int)(dSl+0.5);
			leftover = dSl-nSl;
			for(int i=0; i<nSl; i++)
				newBs->bonds[iPol][jBondNew++] = 0x0;
			
			newBs->bonds[iPol][jBondNew++] = bond;
			curCoor = AddUnitToCoor(bs->bonds[iPol][iBond], curCoor, newBs->L);
		}
		newBs->nMono[iPol] = newNMono[iPol];
	}
	
	for(int coor=0; coor<bs->L*bs->L*bs->L; coor++){
		int newCoor = UpscaleCoor(coor, bs->L);
		newBs->topo[newCoor]    = bs->topo[coor];
		newBs->bondOcc[newCoor] = bs->bondOcc[coor];
	}
	
	for(int i=0; i<4; i++) newBs->seed[i] = bs->seed[i];
	
	
	return newBs;
}

int* LoadStraightTopo(char* file){
	FILE* pFile = fopen(file, "r");
	if(!pFile){ printf("Error opening file %s\n", file); exit(192); }
	int* straightTopo = malloc(sizeof(int)*16);
	for(int i=0; i<16; i++) straightTopo[i] = -1;
	
	int bond, topo;
	while(!feof(pFile) && fscanf(pFile, "%x %i", &bond, &topo)){
		straightTopo[bond] = topo;
	}
	
// 	for(int i=0; i<16; i++) {
// 		if(straightTopo[i] != -1)
// 			printf("%x %i\n", i, straightTopo[i]);
// 	}
	return straightTopo;
}

void WriteState(BoxState* bs, char* file, char* topoFile){
	FILE* pFile = fopen(file, "w");
	if(!pFile){ printf("Error opening file %s\n", file); exit(192); }

	fprintf(pFile, "LT= %i\n", bs->L);
	fprintf(pFile, "LU= %i\n", bs->L);
	fprintf(pFile, "LV= %i\n", bs->L);
	fprintf(pFile, "np= %i\n", bs->nPol);
	fprintf(pFile, "maxPolLength= %i\n", bs->maxPolSize);
	
	for(int iPol=0; iPol<bs->nPol; iPol++){
		fprintf(pFile, "len= %i\n", bs->nMono[iPol]);
		fprintf(pFile, "%i %i %i\t", bs->startTUV[iPol][0], bs->startTUV[iPol][1], bs->startTUV[iPol][2]);
		for(int iBond=0; iBond<bs->nMono[iPol]; iBond++){
			fprintf(pFile, "%x", bs->bonds[iPol][iBond]);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	
	pFile = fopen(topoFile, "w");
	if(!pFile){ printf("Error opening file %s\n", topoFile); exit(192); }

	for(int i=0; i<4; i++) fprintf(pFile, "%u ", bs->seed[i]);
	fprintf(pFile, "\n");
	for(int coor=0; coor<bs->L*bs->L*bs->L; coor++)
		fprintf(pFile, "%i %i\n", bs->topo[coor], bs->bondOcc[coor]);
	fclose(pFile);
}

