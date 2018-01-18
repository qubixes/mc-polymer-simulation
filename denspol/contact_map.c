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
	int nPol;
	int N;
	int L;
}Data;

typedef struct RunProperties{
	char* fileIn;
	char* fileOut;
	long nSamples;
}RunProperties;

Data* ReadData(char* file);
void PrintContactMatrix(Data* data, char* file, long nSamples);

int main(int argc, char** argv){
	if(argc<3){
		printf("Need two arguments!\n");
		return 192;
	}
	
	RunProperties rp;
	
	rp.fileIn = argv[1];
	rp.fileOut = argv[2];
	if(argc > 3)
		rp.nSamples = (long)atof(argv[3]);
	else
		rp.nSamples = (long)1e8;
	
// 	printf("nSamples = %li\n", rp.nSamples);
// 	printf("argc=%i, %s\n", argc, argv[3]);
	Data* data = ReadData(rp.fileIn);
	PrintContactMatrix(data, rp.fileOut, rp.nSamples);
}

Data* NewData(int polSize, int nPol, int L){
	Data* data    = malloc(sizeof(Data));
	data->str     = malloc(sizeof(char)*(polSize+1));
	data->nIdList = malloc(sizeof(double)*(L*L*L));
	data->lattice = malloc(sizeof(DInt*)*(L*L*L));
	for(int coor=0; coor<L*L*L; coor++){
		data->lattice[coor] = malloc(sizeof(DInt)*(100));
		data->nIdList[coor] = 0;
	}
	data->L = L;
	data->nPol = nPol;
	return data;
}

Data* ReadData(char* file){
	int nPol, maxPolSize, L;
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int i=0; i<3; i++)
		fscanf(pFile, "%*s %i", &L);
	
	fscanf(pFile, "%*s %i", &nPol);
	
	fscanf(pFile, "%*s %i", &maxPolSize);
// 	printf("maxPolSize = %i\n", maxPolSize);
	Data* data = NewData(maxPolSize, nPol, L);
	
	for(int iPol=0; iPol<nPol; iPol++){
		int polSize;
		int tuvStart[3];
		fscanf(pFile, "%*s %i %i %i %i %s", &polSize, tuvStart, tuvStart+1, tuvStart+2, data->str);
		data->N = polSize;
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
			
			coor = AddUnitToCoor(bond, coor, L);
			if(bond == 0xf) break;
		}
	}
	fclose(pFile);
	return data;
}

void PrintContactMatrix(Data* data, char* file, long nSamples){
	int L=data->L;
	
	int dAlloc = 1000;
	int cAlloc = dAlloc;
	Contact* contacts = malloc(sizeof(Contact)*cAlloc);
	int nContacts=0;
	unsigned int rng[4];
	Seed(rng, 1209384);
	
	for(int coor=0; coor<L*L*L; coor++){
		for(int indexMono=0; indexMono<data->nIdList[coor]; indexMono++){
			for(int jIndexMono=indexMono+1; jIndexMono<data->nIdList[coor]; jIndexMono++){
				contacts[nContacts].monomers[0].iPol = data->lattice[coor][indexMono].iPol;
				contacts[nContacts].monomers[0].iMono= data->lattice[coor][indexMono].iMono;
				contacts[nContacts].monomers[1].iPol = data->lattice[coor][jIndexMono].iPol;
				contacts[nContacts].monomers[1].iMono= data->lattice[coor][jIndexMono].iMono;
				nContacts++;
				if(cAlloc <= nContacts){
					cAlloc += dAlloc;
					contacts = realloc(contacts, sizeof(Contact)*cAlloc);
				}
			}
		}
	}
	
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "#len= %i\n", data->N);
	for(long i=0; i<nSamples && i<nContacts; i++){
		long j = DRng(rng)*(nContacts-i);
		fprintf(pFile, "%i %i %i %i %lf\n", contacts[j].monomers[0].iPol, contacts[j].monomers[0].iMono, contacts[j].monomers[1].iPol, contacts[j].monomers[1].iMono, 1.0);
		contacts[j] = contacts[nContacts-i-1];
	}
	fclose(pFile);
}
