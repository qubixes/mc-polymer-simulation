#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define IS_MAIN
#include "pascal_lib.h"
#include "rng.h"

#define BOUNDARY_PERIOD 0
#define BOUNDARY_STATIC 1

typedef struct RunProperties{
	char* pasConfigFile;
	char* pasInfoFile;
	char* recFile;
}RunProperties;

double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1);
double* SystemSimilarity(PasData* pData, Data* simData, long maxSamples, unsigned int rng[4]);

/** The internal error is in the first and 3rd column.
  * The external error is in the second and 4th column.
  **/

int main(int argc, char** argv){
	if(argc<4){
		printf("Need three arguments!\n");
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


double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1){
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
		pDxyz[k] = pData->pos[iPol][iMono1][k] - pData->pos[jPol][jMono1][k];
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
	
	double pDist = sqrt(Squabs(pDxyz));
	double simDist = sqrt(simData->TUV2Distance(ituv2, jtuv2, simData->L));
	
	pDist   /= pow(pData->nTotMono/4.87, 1./3.);
	simDist /= pow(simData->nTotMono/10.18, 1./3.);
	
// 	printf("%le %le %le %le %le\n", pDist, Squabs(pDxyz), simDist, pow(pData->nTotMono/4.87, 1./3.), pow(simData->nTotMono/10.18, 1./3.));
	
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

