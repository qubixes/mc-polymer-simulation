#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define IS_MAIN
#include "pascal_lib.h"
#include "rng.h"

typedef struct RunProperties{
	char* pasConfigFile;
	char* pasInfoFile;
	char* tuvFile;
	char* outFile;
}RunProperties;

double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1);
void SystemCompare(PasData* pData, Data* simData, double** pasDist, double** pasDistSq, int** nPasDist);
void WriteCompare(char* file, PasData* pData, Data* simData, double* pasDist, double* pasDistSq, int* nPasDist);

/** The internal error is in the first and 3rd column.
  * The external error is in the second and 4th column.
  **/

int main(int argc, char** argv){
	if(argc<5){
		printf("Need four arguments!\n");
		return 192;
	}
	RunProperties rp;
	
	rp.pasInfoFile = argv[1];
	rp.pasConfigFile = argv[2];
	rp.tuvFile = argv[3];
	rp.outFile = argv[4];
	
	PasData* pData = PascalDataInit(rp.pasInfoFile, rp.pasConfigFile);
	Data* simData = ReadData(rp.tuvFile, BOUNDARY_STATIC);
	
	double *pasDist, *pasDistSq; 
	int *nPasDist;
	SystemCompare(pData, simData, &pasDist, &pasDistSq, &nPasDist);
	WriteCompare(rp.outFile, pData, simData, pasDist, pasDistSq, nPasDist);
}

double MagnificationRatio(int nMono, int nMonoOrig, int polType){
	if(polType == POL_TYPE_LIN)
		return ((nMono-1)/(double)(nMonoOrig-1));
	else
		return (nMono/(double)nMonoOrig);
}

void GetDistancesFromPas(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1, double* pasDist, double* simDist){
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
	
	*pasDist = sqrt(Squabs(pDxyz));
	*simDist = sqrt(simData->TUV2Distance(ituv2, jtuv2, simData->L));
}

void GetDistancesFromSim(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1, double* pasDist, double* simDist){
	double ratioI, ratioJ;
	
	ratioI = MagnificationRatio(pData->nMono[iPol], simData->nMono[iPol], simData->polTypes[iPol]);
	ratioJ = MagnificationRatio(pData->nMono[jPol], simData->nMono[jPol], simData->polTypes[jPol]);
	
	double ixyz2[3], jxyz2[3];
	double tuvI[3], tuvJ[3];
	
	double dIMono2 = iMono1*ratioI;
	double dJMono2 = jMono1*ratioJ;
	
	int iMono2 = (int)dIMono2;
	int jMono2 = (int)dJMono2;
	
	double dI = dIMono2-(int)dIMono2;
	double dJ = dJMono2-(int)dJMono2;
	
// 	if(jMono2 >= pData->nMono[jPol]){
// // 		printf("jMono1 = %i, jMono2 = %i, nMono[old] = %i, nMono[new] = %i, ratio = %lf\n", jMono1, jMono2, simData->nMono[jPol], pData->nMono[jPol], ratioJ);
// 		printf("eeeeeeeeeeeeeuuuuuuuuuuuuuwwwwwwwwwwwww\n");
// 	}
	
	for(int k=0; k<3; k++){
		ixyz2[k]  = (1-dI) * pData->pos[iPol][iMono2  ][k];
		ixyz2[k] +=    dI  * pData->pos[iPol][(iMono2+1)%pData->nMono[iPol]][k];
		jxyz2[k]  = (1-dJ) * pData->pos[jPol][jMono2  ][k];
		jxyz2[k] +=    dJ  * pData->pos[jPol][(jMono2+1)%pData->nMono[jPol]][k];
		
		tuvI[k] = simData->tuv[iPol][iMono1][k];
		tuvJ[k] = simData->tuv[jPol][jMono1][k];
	}
	
	double dxyz[3];
	for(int k=0; k<3; k++)
		dxyz[k] = ixyz2[k]-jxyz2[k];
	
	*pasDist = sqrt(Squabs(dxyz));
	*simDist = sqrt(simData->TUV2Distance(tuvI, tuvJ, simData->L));
}

// double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1){
// 	double pasDist, simDist;
// 	
// 	GetDistances(pData, simData, iPol, jPol, iMono1, jMono1, &pasDist, &simDist);
// 	
// 	pasDist /= pow(pData->nTotMono/4.87, 1./3.);
// 	simDist /= pow(simData->nTotMono/10.18, 1./3.);
// 	
// 	return fabs(pDist-simDist);
// }


void SystemCompare(PasData* pData, Data* simData, double** pasDist, double** pasDistSq, int** nPasDist){
	if(pData->nPol != simData->nPol){
		printf("Error: cannot compare data files with a different number of polymers.\n");
		exit(192);
	}
	int LS = simData->L*simData->L*simData->L;
	
	*pasDist   = malloc(sizeof(double)*LS);
	*pasDistSq = malloc(sizeof(double)*LS);
	*nPasDist  = malloc(sizeof(int)   *LS);
	
	for(int i=0; i<LS; i++) (*nPasDist)[i]=0;
	
	for(int iPol=0; iPol<simData->nPol; iPol++){
		for(int iMono1=0; iMono1<simData->nMono[iPol]; iMono1++){
			for(int jPol=0; jPol<simData->nPol; jPol++){
				for(int jMono1=0; jMono1<simData->nMono[jPol]; jMono1++){
					if(iPol == jPol && iMono1 == jMono1) continue;
					double pDist, sDist;
					
					GetDistancesFromSim(pData, simData, iPol, jPol, iMono1, jMono1, &pDist, &sDist);
					int sSqDistInt          = (int)(sDist*sDist*2+0.5);
					(*pasDist)[sSqDistInt]   += pDist;
					(*pasDistSq)[sSqDistInt] += pDist*pDist;
					(*nPasDist)[sSqDistInt]  ++;
				}
			}
		}
	}
}

void WriteCompare(char* file, PasData* pData, Data* simData, double* pasDist, double* pasDistSq, int* nPasDist){
	
	int LS = simData->L*simData->L*simData->L;
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int simSq=0; simSq<LS; simSq++){
		if(!nPasDist[simSq]) continue;
		double avgPasDist = pasDist[simSq]/(double)nPasDist[simSq];
		double avgPasDistSq = pasDistSq[simSq]/(double)nPasDist[simSq];
		
		double avgPasDistNorm = avgPasDist/pow(pData->nTotMono/4.87, 1./3.);
		double avgPasDistSqNorm = avgPasDistSq/pow(pData->nTotMono/4.87, 2./3.);
		
		double realSimDist = sqrt(simSq/2.0);
		double realSimDistNorm = realSimDist/pow(simData->nTotMono/10.18, 1./3.);
		
		double pasSD = sqrt(avgPasDistSq-avgPasDist*avgPasDist);
		double pasSDNorm = sqrt(avgPasDistSqNorm-avgPasDistNorm*avgPasDistNorm);
		
		fprintf(pFile, "%lf %lf %lf %lf %lf %lf\n", realSimDist, avgPasDist, pasSD, realSimDistNorm, avgPasDistNorm, pasSDNorm);
	}
	fclose(pFile);
}

