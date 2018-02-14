#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define IS_MAIN
#include "lowm_modes.h"
#include "rng.h"

typedef struct Data{
	int*** tuv;
	int nPol;
	int polSize;
	int polType;
	int L;
}Data;

typedef struct ContactData{
	long nContact;
	int polSize;
	int nPol;
	int** contacts;
	double* contactStrengths;
}ContactData;

typedef struct RunProperties{
	char* simFile;
	char* contactFile;
}RunProperties;

Data* ReadData(char* file);
ContactData* ReadContactData(char* file);
double* SimilarityToContact(Data* simData, ContactData* conData);

int main(int argc, char** argv){
	if(argc<3){
		printf("Need two arguments!\n");
		return 192;
	}
	RunProperties rp;
	
	rp.simFile = argv[1];
	rp.contactFile = argv[2];
	
	unsigned int rng[4];
	Seed(rng, 129846412);
	
	Data* simData = ReadData(rp.simFile);
	ContactData* contactData = ReadContactData(rp.contactFile);
	
	if(simData->nPol != contactData->nPol){
		printf("Comparison incompatible between contact matrix and simulation data.\nNeed same polymer size and number of polymers\n");
		exit(192);
	}	
	
	double* similarity = SimilarityToContact(simData, contactData);
	printf("%lf %lf\n", similarity[0], similarity[1]);
}

Data* NewData(int polSize, int nPol, int L){
	Data* data    = malloc(sizeof(Data));
	data->nPol    = nPol;
	data->polSize = polSize;
	data->L       = L;
	
	data->tuv = malloc(sizeof(int**)*nPol);
	for(int i=0; i<nPol; i++){
		data->tuv[i] = malloc(sizeof(int*)*(polSize+1));
		for(int j=0; j<=polSize; j++)
			data->tuv[i][j] = malloc(sizeof(int)*3);
	}
	
	return data;
}

ContactData* NewContactData(int polSize, int nPol, long nContact){
	ContactData* cData = malloc(sizeof(ContactData));
	cData->nPol        = nPol;
	cData->polSize     = polSize;
	cData->nContact    = nContact;
	
	cData->contactStrengths = malloc(sizeof(double)*nContact);
	cData->contacts = malloc(sizeof(int*)*nContact);
	for(long i=0; i<nContact; i++){
		cData->contacts[i] = malloc(sizeof(int)*4);
	}
	return cData;
}

void AddCharToTUV(char c, int tuv[3]){
	int bond = CharToHex(c);
	
	for(int i=0; i<3; i++){
		tuv[i] += (bond>>i)&0x1;
		tuv[i] -= (bond>>3)&0x1;
	}
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
	Data* data = NewData(maxPolSize, nPol, L);
	char* str  = malloc(sizeof(char)*(maxPolSize+1));
	
	for(int iPol=0; iPol<nPol; iPol++){
		int tuv[3];
		fscanf(pFile, "%*s %i %i %i %i %s", &data->polSize, tuv, tuv+1, tuv+2, str);
		for(int iMono=0; iMono<data->polSize; iMono++){
			for(int k=0; k<3; k++)
				data->tuv[iPol][iMono][k] = tuv[k];
			AddCharToTUV(str[iMono], tuv);
		}
		for(int k=0; k<3; k++)
			data->tuv[iPol][data->polSize][k] = tuv[k];
	}
	if(str[data->polSize-1] == 'f') data->polType = POL_TYPE_LIN;
	else data->polType = POL_TYPE_RING;
	fclose(pFile);
	return data;
}

ContactData* ReadContactData(char* file){
	int nPol, maxPolSize;
	long nContact;
	
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	fscanf(pFile, "%*s %i" , &nPol);
	fscanf(pFile, "%*s %i" , &maxPolSize);
	fscanf(pFile, "%*s %li", &nContact);
	ContactData* cData = NewContactData(maxPolSize, nPol, nContact);
	
	for(long ic=0; ic<nContact; ic++){
		fscanf(pFile, "%i %i %i %i %lf", cData->contacts[ic], cData->contacts[ic]+1, cData->contacts[ic]+2, cData->contacts[ic]+3, cData->contactStrengths+ic);
	}
	fclose(pFile);
	return cData;
}


void DblTUV2XYZ(double tuv[3], double xyz[3]){
	xyz[0] = tuv[0]+tuv[1]-tuv[2];
	xyz[1] = tuv[0]-tuv[1]       ;
	xyz[2] =               tuv[2];
}

double Squabs(double xyz[3]){
	int tot=0;
	for(int k=0; k<3; k++) tot += xyz[k]*xyz[k];
	return tot;
}

/// Calculate the distance between two points on the lattice with periodic 
/// boundary conditions.
/// 
/// It's looks a bit weird, but dtuv starts at the (-L,-L,-L) - (0,0,0) periodic box 
/// for slightly easier indexing.
/// 
/// One important note: it is destructive, so don't put original data in!

double TUV2Distance(double tuv1[3], double tuv2[3], int L){
	double dxyz[3];
	
	double dtuv[3], newtuv[3];
	
	
	
	for(int k=0; k<3; k++){
		while(tuv1[k] <  0) tuv1[k] += L;
		while(tuv1[k] >= L) tuv1[k] -= L;
		while(tuv2[k] <  0) tuv2[k] += L;
		while(tuv2[k] >= L) tuv2[k] -= L;
		dtuv[k] = tuv1[k]-tuv2[k]-L;
	}
	
	double minDist=L*L*L*L;
	
	newtuv[0]=dtuv[0];
	for(int dt=-L; dt <= L; dt += L, newtuv[0] += L){
		newtuv[1]=dtuv[1];
		for(int du=-L; du <= L; du += L, newtuv[1] += L){
			newtuv[2]=dtuv[2];
			for(int dv=-L; dv <= L; dv += L, newtuv[2] += L){
				DblTUV2XYZ(newtuv, dxyz);
				double newDist = Squabs(dxyz);
				minDist = MIN(minDist, newDist);
			}
		}
	}
	
	return minDist;
}

double* SimilarityToContact(Data* simData, ContactData* conData){
	double* similarity = malloc(sizeof(double)*2);
	double* nSample = malloc(sizeof(long)*2);
	double ratio = simData->polSize/(double)conData->polSize;
	double ituv[3], jtuv[3];
	
	for(int i=0; i<2; i++){
		similarity[i]=0; nSample[i]=0;
	}
	
	for(int iCon=0; iCon<conData->nContact; iCon++){
		int iPol  = conData->contacts[iCon][0];
		int iMono = (int)(conData->contacts[iCon][1]*ratio+0.5);
		int jPol  = conData->contacts[iCon][2];
		int jMono = (int)(conData->contacts[iCon][3]*ratio+0.5);
		
		for(int k=0; k<3; k++){
			ituv[k] = simData->tuv[iPol][iMono][k];
			jtuv[k] = simData->tuv[jPol][jMono][k];
		}
		
		double distNorm = sqrt(TUV2Distance(ituv, jtuv, simData->L))/(double)(simData->L);
		int k=((iPol == jPol)?0:1);
		similarity[k] += distNorm;
		nSample[k]++;
	}
	
	for(int k=0; k<2; k++) similarity[k] /= nSample[k];
	free(nSample);
	return similarity;
}

