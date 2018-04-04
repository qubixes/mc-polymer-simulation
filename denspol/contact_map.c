#include <stdlib.h>
#include <stdio.h>
#define IS_MAIN
#include "denspol_lib.h"

typedef struct DInt{
	int iPol, iMono;
}DInt;

typedef struct Contact{
	DInt monomers[2];
}Contact;

typedef struct Data{
	DInt** lattice;
	int* nIdList;
	char* str;
	int* polTypes;
	int* polSizes;
	int nPol;
	int maxPolSize;
	int L;
}Data;

typedef struct RunProperties{
	char* fileIn;
	char* fileOut;
	unsigned int seed;
	long nSamples;
}RunProperties;

Data* ReadData(char* file);
void PrintContactMatrix(Data* data, RunProperties* rp);

int main(int argc, char** argv){
	if(argc<3){
		printf("Need two arguments!\n");
		return 192;
	}
	
	RunProperties rp;
	rp.nSamples = (long)1e8;
	rp.seed     = 1209384;
	
	rp.fileIn = argv[1];
	rp.fileOut = argv[2];
	if(argc > 3)
		rp.seed = atoi(argv[3]);
	if(argc > 4)
		rp.nSamples = (long)atof(argv[4]);
	
	Data* data = ReadData(rp.fileIn);
	PrintContactMatrix(data, &rp);
}

Data* NewData(int polSize, int nPol, int L){
	Data* data    = malloc(sizeof(Data));
	data->str     = malloc(sizeof(char)*(polSize+1));
	data->nIdList = malloc(sizeof(double)*(L*L*L));
	data->lattice = malloc(sizeof(DInt*)*(L*L*L));
	data->polTypes= malloc(sizeof(int)*nPol);
	data->polSizes= malloc(sizeof(int)*nPol);
	for(int coor=0; coor<L*L*L; coor++){
		data->lattice[coor] = malloc(sizeof(DInt)*(100));
		data->nIdList[coor] = 0;
	}
	data->L = L;
	data->nPol = nPol;
	return data;
}

Data* ReadData(char* file){
	int nPol, maxLen, L;
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int i=0; i<3; i++)
		fscanf(pFile, "%*s %i", &L);
	
	fscanf(pFile, "%*s %i", &nPol);
	
	fscanf(pFile, "%*s %i", &maxLen);
// 	printf("maxPolSize = %i\n", maxPolSize);
	Data* data = NewData(maxLen, nPol, L);
	data->maxPolSize=0;
	for(int iPol=0; iPol<nPol; iPol++){
		int tuvStart[3];
		fscanf(pFile, "%*s %i %i %i %i %s", &data->polSizes[iPol], tuvStart, tuvStart+1, tuvStart+2, data->str);
		int polSize = data->polSizes[iPol];
		data->maxPolSize = MAX(polSize,data->maxPolSize);
		data->polTypes[iPol] = (data->str[polSize-1] == 'f')?POL_TYPE_LIN:POL_TYPE_RING;
		int coor = TUV2Coor(tuvStart[0], tuvStart[1], tuvStart[2], L);
		int bond;
		int zeroBefore=0;
		int zero=0;
		while( (bond = CharToHex(data->str[polSize-1-zeroBefore])) == 0) zeroBefore++;
		for(int iMono=0; iMono<polSize-zeroBefore; iMono+=zero+1){
			zero=0;
			while( (bond = CharToHex(data->str[iMono+zero])) == 0){
				zero++;
			}
			
			int monoId;
			if(iMono == 0){
				monoId = (iMono+(zero-(zeroBefore+zero)/2+polSize))%polSize;
			}
			else
				monoId = iMono+zero/2;
			
			data->lattice[coor][data->nIdList[coor]  ].iMono = monoId;
			data->lattice[coor][data->nIdList[coor]++].iPol  = iPol;
			
			int OOB; ///NOTE: please fix this, check out of bounds set boundary conditions.
			coor = AddUnitToCoorPeriod(bond, coor, L, &OOB);
			if(bond == 0xf) break;
		}
	}
	fclose(pFile);
	return data;
}

void PrintContactMatrix(Data* data, RunProperties* rp){
	int L=data->L;
	
	int dAlloc = 1000;
	int cAlloc = dAlloc;
	Contact* contacts = malloc(sizeof(Contact)*cAlloc);
	int nContacts=0;
	unsigned int rng[4];
	Seed(rng, rp->seed);
	
	for(int coor=0; coor<L*L*L; coor++){
		for(int indexMono=0; indexMono<data->nIdList[coor]; indexMono++){
			for(int jIndexMono=indexMono+1; jIndexMono<data->nIdList[coor]; jIndexMono++){
				contacts[nContacts].monomers[0].iPol  = data->lattice[coor][indexMono].iPol;
				contacts[nContacts].monomers[0].iMono = data->lattice[coor][indexMono].iMono;
				contacts[nContacts].monomers[1].iPol  = data->lattice[coor][jIndexMono].iPol;
				contacts[nContacts].monomers[1].iMono = data->lattice[coor][jIndexMono].iMono;
				nContacts++;
				if(cAlloc <= nContacts){
					cAlloc += dAlloc;
					contacts = realloc(contacts, sizeof(Contact)*cAlloc);
				}
			}
		}
	}
	
	FILE* pFile = fopen(rp->fileOut, "w");
	if(!pFile){
		printf("Contact_map: error opening file %s for writing\n", rp->fileOut);
		exit(192);
	}
	fprintf(pFile, "#nPol= %i\n", data->nPol);
	fprintf(pFile, "#maxlen= %i\n", data->maxPolSize);
	fprintf(pFile, "#nContacts= %li\n", MIN(rp->nSamples, nContacts));
	for(int i=0; i<data->nPol; i++){
		if(data->polTypes[i] == POL_TYPE_LIN)
			fprintf(pFile, "lin");
		else if(data->polTypes[i] == POL_TYPE_RING)
			fprintf(pFile, "ring");
		else
			fprintf(pFile, "???");
		fprintf(pFile, " %i\n", data->polSizes[i]);
	}
		
	for(long i=0; i<rp->nSamples && i<nContacts; i++){
		long j = DRng(rng)*(nContacts-i);
		fprintf(pFile, "%i %i %i %i %lf\n", contacts[j].monomers[0].iPol, contacts[j].monomers[0].iMono, contacts[j].monomers[1].iPol, contacts[j].monomers[1].iMono, 1.0);
		contacts[j] = contacts[nContacts-i-1];
	}
	fclose(pFile);
}
