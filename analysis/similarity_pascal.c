#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define IS_MAIN
#include "pascal_lib.h"
#include "lowm_modes.h"
#include "rng.h"

#define BOUNDARY_PERIOD 0
#define BOUNDARY_STATIC 1

typedef struct Data{
	int*** tuv;
	int* nMono;
	int* polTypes;
	double (*TUV2Distance) (double*, double*, int);
	int nPol;
	int maxNMono;
	int nTotMono;
	int L;
}Data;

typedef struct RunProperties{
	char* pasConfigFile;
	char* pasInfoFile;
	char* recFile;
}RunProperties;

Data* ReadData(char* file, int boundaryCond);
double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1);
double* SystemSimilarity(PasData* pData, Data* simData, long maxSamples, unsigned int rng[4]);

double TUV2DistancePeriodic(double tuv1[3], double tuv2[3], int L);
double TUV2DistanceStatic(double tuv1[3], double tuv2[3], int L);

/** The internal error is in the first and 3rd column.
  * The external error is in the second and 4th column.
  **/

int main(int argc, char** argv){
	if(argc<3){
		printf("Need two arguments!\n");
		return 192;
	}
	RunProperties rp;
	
	rp.pasInfoFile = argv[1];
	rp.pasConfigFile = argv[2];
	rp.recFile = argv[3];
	
	unsigned int rng[4];
	Seed(rng, 129846412);
	
	PasData* pData = PascalDataInit(rp.pasInfoFile, rp.pasConfigFile);
	Data* simData = ReadData(rp.recFile, BOUNDARY_STATIC);
	long maxSamples = 1e6;
	
	double* simForw = SystemSimilarity(pData, simData, maxSamples, rng);
	
	printf("%lf %lf %lf %lf\n", simForw[0], simForw[1], simForw[0], simForw[1]);
}

Data* NewData(int maxNMono, int nPol, int L, int boundaryCond){
	Data* data    = malloc(sizeof(Data));
	data->nPol    = nPol;
	data->maxNMono= maxNMono;
	data->L       = L;
	
	data->tuv = malloc(sizeof(int**)*nPol);
	for(int i=0; i<nPol; i++){
		data->tuv[i] = malloc(sizeof(int*)*maxNMono);
		for(int j=0; j<maxNMono; j++)
			data->tuv[i][j] = malloc(sizeof(int)*3);
	}
	
	data->nMono    = malloc(sizeof(int)*nPol);
	data->polTypes = malloc(sizeof(int)*nPol);
	
	if(boundaryCond == BOUNDARY_PERIOD)
		data->TUV2Distance = &TUV2DistancePeriodic;
	else
		data->TUV2Distance = &TUV2DistanceStatic;
	
	return data;
}

void AddCharToTUV(char c, int tuv[3]){
	int bond = CharToHex(c);
	
	for(int i=0; i<3; i++){
		tuv[i] += (bond>>i)&0x1;
		tuv[i] -= (bond>>3)&0x1;
	}
}

Data* ReadData(char* file, int boundaryCond){
	int nPol, maxNMono, L;
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int i=0; i<3; i++)
		fscanf(pFile, "%*s %i", &L);
	
	fscanf(pFile, "%*s %i", &nPol);
	
	fscanf(pFile, "%*s %i", &maxNMono);
	Data* data = NewData(maxNMono, nPol, L, boundaryCond);
	char* str  = malloc(sizeof(char)*(maxNMono+1));
	
	data->nTotMono=0;
	for(int iPol=0; iPol<nPol; iPol++){
		int tuv[3];
		fscanf(pFile, "%*s %i %i %i %i %s", &data->nMono[iPol], tuv, tuv+1, tuv+2, str);
		data->nTotMono += data->nMono[iPol];
		for(int iMono=0; iMono<data->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++)
				data->tuv[iPol][iMono][k] = tuv[k];
			AddCharToTUV(str[iMono], tuv);
		}
		if(str[data->nMono[iPol]] == 'f') data->polTypes[iPol] = POL_TYPE_LIN;
		else data->polTypes[iPol] = POL_TYPE_RING;
	}
	fclose(pFile);
	return data;
}

void DblTUV2XYZ(double tuv[3], double xyz[3]){
	xyz[0] = (tuv[0]+tuv[1]-tuv[2])/sqrt(2);
	xyz[1] = (tuv[0]-tuv[1]       )/sqrt(2);
	xyz[2] = (              tuv[2])/sqrt(2);
}

double Squabs(double xyz[3]){
	int tot=0;
	for(int k=0; k<3; k++) tot += xyz[k]*xyz[k];
	return tot;
}

/// Calculate the distance between two points on the lattice with periodic boundary conditions.
/// 
/// It's looks a bit weird, but dtuv starts at the (-L,-L,-L) - (0,0,0) periodic box 
/// for slightly easier indexing.

double TUV2DistancePeriodic(double tuv1[3], double tuv2[3], int L){
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

double TUV2DistanceStatic(double tuv1[3], double tuv2[3], int L){
	double dtuv[3], dxyz[3];
	
	for(int k=0; k<3; k++)
		dtuv[k] = tuv1[k]-tuv2[k];
	DblTUV2XYZ(dtuv, dxyz);
	return Squabs(dxyz);
}

double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1){
// 	double ratio = simData->nMono[iPol]/(double)pData->nMono[jPol];
	double ratioI, ratioJ;
	
	if(simData->polTypes[iPol] == POL_TYPE_LIN)
		ratioI = (simData->nMono[iPol]-1)/(double)(pData->nMono[iPol]-1);
	else
		ratioI = simData->nMono[iPol]/(double)pData->nMono[iPol];
	
	if(simData->polTypes[jPol] == POL_TYPE_LIN)
		ratioJ = (simData->nMono[jPol]-1)/(double)(pData->nMono[jPol]-1);
	else
		ratioJ = simData->nMono[jPol]/(double)pData->nMono[jPol];
	
	double ituv2[3], jtuv2[3];
	
	double pDxyz[3];
	
	for(int k=0; k<3; k++){
		pDxyz[k] += pData->pos[iPol][iMono1][k] - pData->pos[jPol][jMono1][k];
	}
	
	double dIMono2 = iMono1*ratioI;
	double dJMono2 = jMono1*ratioJ;
	
	int iMono2 = (int)dIMono2;
	int jMono2 = (int)dJMono2;
	
	double dI = dIMono2-(int)dIMono2;
	double dJ = dJMono2-(int)dJMono2;
	
	for(int k=0; k<3; k++){
		ituv2[k]  = (1-dI) * simData->tuv[iPol][iMono2  ][k];
		ituv2[k] +=    dI  * simData->tuv[iPol][(iMono2+1)%simData->nMono[iPol]][k];
		jtuv2[k]  = (1-dJ) * simData->tuv[jPol][jMono2  ][k];
		jtuv2[k] +=    dJ  * simData->tuv[jPol][(jMono2+1)%simData->nMono[jPol]][k];
	}
	
	double pDist = Squabs(pDxyz);
	double simDist = simData->TUV2Distance(ituv2, jtuv2, simData->L);
	
	pDist   /= pow(pData->nTotMono/4.87, 1./3.);
	simDist /= pow(simData->nTotMono/10.18, 1./3.);
	
	return fabs(pDist-simDist);
}

double* SystemSimilarity(PasData* pData, Data* simData, long maxSamples, unsigned int rng[4]){
	if(pData->nPol != simData->nPol){
		printf("Error: cannot compare data files with a different number of polymers.\n");
		exit(192);
	}
	
	double* similarity=malloc(sizeof(double)*2);
	long* nSample = malloc(sizeof(long)*2);
	for(int k=0; k<2; k++){
		similarity[k]=0;
		nSample[k]=0;
	}
	
	long maxPosSamples = (long)pData->nPol*pData->nPol*pData->maxNMono*pData->maxNMono;
	long maxInterSamples = (long)pData->maxNMono*pData->maxNMono*pData->nPol;
	
	if(maxPosSamples < maxSamples){
// 		printf("Full comparison: %li < %li\n", maxPosSamples, maxSamples);
		for(int iPol=0; iPol<pData->nPol; iPol++){
			for(int iMono1=0; iMono1<pData->nMono[iPol]; iMono1++){
				for(int jPol=0; jPol<pData->nPol; jPol++){
					for(int jMono1=0; jMono1<pData->nMono[jPol]; jMono1++){
						if(iPol == jPol && iMono1 == jMono1) continue;
						int k= (iPol==jPol)?0:1;
						similarity[k] += Similarity(pData, simData, iPol, jPol, iMono1, jMono1);
						nSample[k]++;
					}
				}
			}
		}
	}
	else if(pData->nPol == 1){
// 		printf("One polymer, only sampling\n");
		int iPol=0, jPol=0;
		while(nSample[0] < maxSamples){
			int iMono1 = pData->nMono[0]*DRng(rng);
			int jMono1 = pData->nMono[0]*DRng(rng);
			if(iMono1 == jMono1) continue;
			similarity[0] += Similarity(pData, simData, iPol, jPol, iMono1, jMono1);
			nSample[0]++;
		}
	}
	else if(maxInterSamples<=maxSamples/2){
// 		printf("Full comparison intra polymers, sampling between\n");
		for(int iPol=0; iPol<pData->nPol; iPol++){
			int jPol = iPol;
			for(int iMono1=0; iMono1<pData->nMono[iPol]; iMono1++){
				for(int jMono1=0; jMono1<pData->nMono[jPol]; jMono1++){
					if(iMono1==jMono1) continue;
					similarity[0] += Similarity(pData, simData, iPol, jPol, iMono1, jMono1);
					nSample[0]++;
				}
			}
		}
		
		long samplesPerCombi= MAX(1, maxSamples/(pData->nPol*(pData->nPol-1)));
		
		for(int iPol=0; iPol<pData->nPol; iPol++){
			for(int jPol=iPol+1; jPol<pData->nPol; jPol++){
				long curSamples=0;
				while(curSamples<samplesPerCombi){
					int iMono1 = pData->nMono[iPol]*DRng(rng);
					int jMono1 = pData->nMono[jPol]*DRng(rng);
					similarity[1] += Similarity(pData, simData, iPol, jPol, iMono1, jMono1);
					nSample[1]++;
					curSamples++;
				}
			}
		}
	}
	else{
// 		printf("Only sampling\n");
		long intraSamplesPerPol = MAX(1, maxSamples/(2*pData->nPol));
		long interSamplesPerCombi = MAX(1, maxSamples/(pData->nPol*(pData->nPol-1)));
		
		for(int iPol=0; iPol<pData->nPol; iPol++){
			int jPol=iPol;
			long curSamples=0;
			while(curSamples<intraSamplesPerPol){
				int iMono1 = pData->nMono[iPol]*DRng(rng);
				int jMono1 = pData->nMono[jPol]*DRng(rng);
				if(iMono1 == jMono1) continue;
				similarity[0] += Similarity(pData, simData, iPol, jPol, iMono1, jMono1);
				nSample[0]++;
				curSamples++;
			}
		}
		
		for(int iPol=0; iPol<pData->nPol; iPol++){
			for(int jPol=iPol+1; jPol<simData->nPol; jPol++){
				long curSamples=0; 
				while(curSamples<interSamplesPerCombi){
				int iMono1 = pData->nMono[iPol]*DRng(rng);
				int jMono1 = pData->nMono[jPol]*DRng(rng);
				similarity[1] += Similarity(pData, simData, iPol, jPol, iMono1, jMono1);
				nSample[1]++;
				curSamples++;
				}
			}
		}
	}
	
	for(int k=0; k<2; k++){
		similarity[k] /= nSample[k];
	}
	
	return similarity;
}

