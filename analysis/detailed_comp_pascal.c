#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define IS_MAIN
#include "pascal_lib.h"
#include "rng.h"

// #define MAX(X,Y) ((X>Y)?X:Y)

typedef struct RunProperties{
	char* pasConfigFile;
	char* pasInfoFile;
	char* tuvFile;
	char* outFile;
	char* pasHistFile;
	char* simHistFile;
}RunProperties;

typedef struct Results{
	double* pasDist, *pasDistSq;
	int* nPasDist;
	double* radDist;
	double* pasRadDist;
}Results;

double Similarity(PasData* pData, Data* simData, int iPol, int jPol, int iMono1, int jMono1);
void SystemCompare(PasData* pData, Data* simData, long maxSamples, unsigned int rng[4], Results* res);
void WriteCompare(char* file, PasData* pData, Data* simData, Results* res);
void GetRadialDistances(Results* res, PasData* pData, Data* simData, Histogram* pasHist, Histogram* simHist);
void WriteHistogram(Histogram* hist, char* file);

/** The internal error is in the first and 3rd column.
  * The external error is in the second and 4th column.
  **/

int main(int argc, char** argv){
	if(argc<7){
		printf("Need six arguments!\n");
		return 192;
	}
	RunProperties rp;
	long maxSamples=(long)1e8;
	unsigned int seed = 198740511;
	unsigned int rng[4];
	
	Seed(rng, seed);
	
	rp.pasInfoFile = argv[1];
	rp.pasConfigFile = argv[2];
	rp.tuvFile = argv[3];
	rp.outFile = argv[4];
	rp.pasHistFile = argv[5];
	rp.simHistFile = argv[6];
	
	PasData* pData = PascalDataInit(rp.pasInfoFile, rp.pasConfigFile);
	Data* simData = ReadData(rp.tuvFile, BOUNDARY_STATIC);
	
	Results res;
	SystemCompare(pData, simData, maxSamples, rng, &res);
	WriteCompare(rp.outFile, pData, simData, &res);
	
	Histogram pasHist, simHist;
	GetRadialDistances(&res, pData, simData, &pasHist, &simHist);
	WriteHistogram(&pasHist, rp.pasHistFile);
	WriteHistogram(&simHist, rp.simHistFile);
}

double MagnificationRatio(int nMono, int nMonoOrig, int polType){
	if(polType == POL_TYPE_LIN)
		return ((nMono-1)/(double)(nMonoOrig-1));
	else
		return (nMono/(double)nMonoOrig);
}

double XYZDist(double xyz1[3], double xyz2[3]){
	double xyz[3];
	for(int k=0; k<3; k++) xyz[k] = xyz1[k]-xyz2[k];
	return Squabs(xyz);
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


void HistogramInit(Histogram* hist, int nBins, double dBin){
	hist->nBin = nBins;
	hist->dBin = dBin;
	hist->count = malloc(sizeof(int)*nBins);
	for(int i=0; i<nBins; i++)
		hist->count[i]=0;
	hist->totCount=0;
}

void WriteHistogram(Histogram* hist, char* file){
	
	int maxBin=0;
	for(int i=hist->nBin-1; i>0; i--){
		if(hist->count[i]){ maxBin=i+1; break; }
	}
	if(maxBin == 0){
		printf("Error: maximum bin == 0\n");
		exit(912);
	}
	if(maxBin<hist->nBin) maxBin++;
	
// 	printf("maxBin = %i, nBin=%i\n", maxBin, hist->nBin);
	FILE* pFile = fopen(file, "w");
	
	for(int i=0; i<maxBin; i++){
		fprintf(pFile, "%lf %lf\n", (i+0.5)*hist->dBin, hist->count[i]/(double)hist->totCount);
	}
	fclose(pFile);
}
	

void GetRadialDistances(Results* res, PasData* pData, Data* simData, Histogram* pasHist, Histogram* simHist){
	double xyzMin[3] = { 3e33, 3e33, 3e33};
	double xyzMax[3] = {-3e33,-3e33,-3e33};
	double dxyz[3];
	double cmsSim[3] = {0,0,0};
	/// First find the radial distances for the simulation data. 
	for(int iPol=0; iPol<simData->nPol; iPol++){
		for(int iMono=0; iMono<simData->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++){
				xyzMin[k] = MIN(simData->xyz[iPol][iMono][k], xyzMin[k]);
				xyzMax[k] = MAX(simData->xyz[iPol][iMono][k], xyzMax[k]);
				cmsSim[k] += simData->xyz[iPol][iMono][k];
			}
		}
	}
	
	for(int k=0; k<3; k++){
		cmsSim[k] /= simData->nTotMono;
		dxyz[k] = xyzMax[k]-xyzMin[k];
	}
	
	double maxSimDist = Squabs(dxyz)/2.0;
	int simNBins = MAX(MIN(300, simData->nTotMono/10), 10);
	double simDBin = maxSimDist/(simNBins-1);
	HistogramInit(simHist, simNBins, simDBin);
// 	printf("maxDist = %lf\n", maxSimDist);
	
	for(int iPol=0; iPol<simData->nPol; iPol++){
		for(int iMono=0; iMono<simData->nMono[iPol]; iMono++){
			double dist= XYZDist(simData->xyz[iPol][iMono], cmsSim);
			int bin = (int)(dist/maxSimDist*simNBins);
// 			printf("bin = %i/%i\n", bin, simNBins);
			simHist->count[bin]++;
			simHist->totCount++;
		}
	}
	
// 	for(int i=0; i<simHist->nBins; i++){
// 		if(simHist->counts[i]) printf("%i %i\n", i, simHist->counts[i]);
// 	}
	
	/// Find radial distances for Pascal data.
	double cmsPas[3] = {0,0,0};
	for(int k=0; k<3; k++){
		xyzMin[k]= 3e33;
		xyzMax[k]=-3e33;
		cmsPas[k]=0;
	}
	
	for(int iPol=0; iPol<pData->nPol; iPol++){
		for(int iMono=0; iMono<pData->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++){
				xyzMin[k] = MIN(pData->pos[iPol][iMono][k], xyzMin[k]);
				xyzMax[k] = MAX(pData->pos[iPol][iMono][k], xyzMax[k]);
				cmsPas[k] += pData->pos[iPol][iMono][k];
			}
		}
	}
	
	for(int k=0; k<3; k++){
		cmsPas[k] /= pData->nTotMono;
		dxyz[k] = xyzMax[k]-xyzMin[k];
	}
	
	double maxPasDist = Squabs(dxyz)/2.0;
	int pasNBins = MAX(MIN(300, pData->nTotMono/10),   10);
	double pasDBin = maxPasDist/(pasNBins-1);
	HistogramInit(pasHist, pasNBins, pasDBin);
	
	for(int iPol=0; iPol<pData->nPol; iPol++){
		for(int iMono=0; iMono<pData->nMono[iPol]; iMono++){
			double dist= XYZDist(pData->pos[iPol][iMono], cmsPas);
			int bin = (int)(dist/maxPasDist*pasNBins);
			pasHist->count[bin]++;
			pasHist->totCount++;
		}
	}
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

void MonoStoreRes(Results* res, PasData* pData, Data* simData, int iPol, int iMono, int jPol, int jMono){
	double pDist, sDist;
	
	GetDistancesFromSim(pData, simData, iPol, jPol, iMono, jMono, &pDist, &sDist);
	int sSqDistInt              = (int)(sDist*sDist*2+0.5);
	res->pasDist[sSqDistInt]   += pDist;
	res->pasDistSq[sSqDistInt] += pDist*pDist;
	res->nPasDist[sSqDistInt]  ++;
}
	

void SystemCompare(PasData* pData, Data* simData, long maxSamples, unsigned int rng[4], Results* res){
	if(pData->nPol != simData->nPol){
		printf("Error: cannot compare data files with a different number of polymers.\n");
		exit(192);
	}
	
	long nTotalDist = (long)simData->nTotMono*simData->nTotMono;
	int LS = simData->L*simData->L*simData->L;
	
	res->pasDist   = malloc(sizeof(double)*LS);
	res->pasDistSq = malloc(sizeof(double)*LS);
	res->nPasDist  = malloc(sizeof(int)   *LS);
	
	for(int i=0; i<LS; i++) res->nPasDist[i]=0;
	
	/** First do intra polymer computation **/
	
	for(int iPol=0; iPol<simData->nPol; iPol++){
		if(simData->nMono[iPol]*simData->nMono[iPol] < maxSamples/(double)simData->nPol){
			for(int iMono=0; iMono<simData->nMono[iPol]; iMono++){
				for(int jMono=0; jMono<simData->nMono[iPol]; jMono++){
					if(iMono == jMono) continue;
					MonoStoreRes(res, pData, simData, iPol, iMono, iPol, jMono);
				}
			}
		}
		else{
			int iMono = (int)(DRng(rng)*simData->nMono[iPol]);
			int jMono = (int)(DRng(rng)*simData->nMono[iPol]);
			MonoStoreRes(res, pData, simData, iPol, iMono, iPol, jMono);
		}
	}
	
	/** Then inter polymer measurement **/
	
	double sampleFac = maxSamples/(double)(nTotalDist);
	for(int iPol=0; iPol<simData->nPol; iPol++){
		for(int jPol=0; jPol<simData->nPol; jPol++){
			if(iPol == jPol) continue;
			if(nTotalDist<maxSamples){
				for(int iMono=0; iMono<simData->nMono[iPol]; iMono++){
					for(int jMono=0; jMono<simData->nMono[jPol]; jMono++){
						MonoStoreRes(res, pData, simData, iPol, iMono, jPol, jMono);
					}
				}
			}
			else{
				long curSamples = MAX(1,(long)(sampleFac*simData->nMono[iPol]*simData->nMono[jPol]));
				for(int iSample=0; iSample<curSamples; iSample++){
					int iMono = (int)(DRng(rng)*simData->nMono[iPol]);
					int jMono = (int)(DRng(rng)*simData->nMono[jPol]);
					MonoStoreRes(res, pData, simData, iPol, iMono, jPol, jMono);
				}
			}
		}
	}
}

void WriteCompare(char* file, PasData* pData, Data* simData, Results* res){
	
	int LS = simData->L*simData->L*simData->L;
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int simSq=0; simSq<LS; simSq++){
		if(!res->nPasDist[simSq]) continue;
		double avgPasDist = res->pasDist[simSq]/(double)res->nPasDist[simSq];
		double avgPasDistSq = res->pasDistSq[simSq]/(double)res->nPasDist[simSq];
		
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

