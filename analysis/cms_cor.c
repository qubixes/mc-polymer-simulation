#include "lowm_modes.h"
#include "file.h"
#include "timer.h"

/** WARNING This program is broken for files with unequal lengths for polymers. **/

typedef struct PartitionTable{
	int** polList;
	int* nPol;
	int nBox;
	int nMaxPol;
	int boxL;
	int xBox;
// 	CMSData* cms;
}PartitionTable;

typedef struct CMSData{
	double** TUV;
	double** XYZ;
	double** TUVRemPer;
	double** XYZRemPer;
}CMSData;

typedef struct Histogram{
	double** avgCor;
	int** counts;
	int* dt;
	int nTime;
	int nBins;
	double binSize;
}Histogram;

typedef struct AvgPosition{
	double **avgPos;
}AvgPosition;

void ComputeCorrelations(SimProperties* sp);
void LoadCMS(long t, CMSData* cms, SimProperties* sp);
CMSData* CMSDataNew(SimProperties* sp);
Histogram* HistogramNew(SimProperties* sp, int nBins, double binSize, TDTTable* tTable);
PartitionTable* PartitionNew(SimProperties* sp, int boxL, int nMaxPol);
double** LoadAvgPos(SimProperties* sp);
void CreatePartition(PartitionTable* pt, CMSData* cms, SimProperties* sp);
void ComputeCorrelations(SimProperties* sp);
void CalcDiff(CMSData* cms, CMSData* cmsDest, long t1, long t2, double** avgPos, double** diff, SimProperties* sp);
void AddToHistogram(Histogram* hst, double drInit, long iDt, double dtDiff);
void ProcessPartition(PartitionTable* pt, CMSData* cms, double** diff, SimProperties* sp, Histogram* hst, int iDt);
void WriteHistogram(Histogram* hst, SimProperties* sp);


int main(int argc, char** argv){
	SimProperties sp;
	char* sampleDir;
	char exec[10000];
	
	
	if(argc<2){ 
		printf("Need basedir\n");
		exit(0);
	}
	sampleDir = argv[1];
	
	sprintf(exec, "mkdir -p %s/cms", sampleDir); system(exec);
	int update =NeedsUpdatePath("cms/t=0_dev=0.res", "cms_cor.dat", sampleDir);
	    update*=NeedsUpdatePath("avg_pos.dat", "cms_cor.dat", sampleDir);
	if(!update){return 0;}
	SetSimProps(&sp, sampleDir);
	for(int iDev=0; iDev<sp.nDev; iDev++){
		ComputeCorrelations(&sp);
// 		ConvertData(&sp, iDev);
	}
	return 0;
}

void LoadCMS(long t, CMSData* cms, SimProperties* sp){
	char file[10000];
	sprintf(file, "%s/cms/t=%li_dev=0.res", sp->sampleDir, t*sp->dT);
	FILE* pFile = fopen(file, "r");
	if(!pFile){ printf("Error opening file %s\n", file); exit(192); }
	
	for(int i=0; i<sp->nPol; i++){
		fscanf(pFile, "%lf %lf %lf\n", cms->TUV[i], cms->TUV[i]+1, cms->TUV[i]+2);
		DTUV2XYZ(cms->TUV[i], cms->XYZ[i]);
		for(int k=0; k<3; k++){
			cms->TUVRemPer[i][k] = fmod(cms->TUV[i][k]+1000*LT, LT);
		}
		DTUV2XYZ(cms->TUVRemPer[i], cms->XYZRemPer[i]);
	}
	fclose(pFile);
}

CMSData* CMSDataNew(SimProperties* sp){
	CMSData* cms = malloc(sizeof(CMSData));
	cms->TUV = malloc(sizeof(double*)*sp->nPol);
	cms->XYZ = malloc(sizeof(double*)*sp->nPol);
	cms->TUVRemPer = malloc(sizeof(double*)*sp->nPol);
	cms->XYZRemPer = malloc(sizeof(double*)*sp->nPol);
	for(int i=0; i<sp->nPol; i++){
		cms->TUV[i] = malloc(sizeof(double)*3);
		cms->XYZ[i] = malloc(sizeof(double)*3);
		cms->TUVRemPer[i] = malloc(sizeof(double)*3);
		cms->XYZRemPer[i] = malloc(sizeof(double)*3);
	}
	return cms;
}

Histogram* HistogramNew(SimProperties* sp, int nBins, double binSize, TDTTable* tTable){
	Histogram* hst = malloc(sizeof(Histogram));
	int nDt = tTable->nDt;
	hst->nTime = nDt;
	hst->nBins = nBins;
	hst->binSize = binSize;
	
	hst->counts = malloc(sizeof(int*)*nDt);
	hst->avgCor = malloc(sizeof(double*)*nDt);
	hst->dt = malloc(sizeof(int)*nDt);
	
	for(int i=0; i<tTable->nTDT; i++){
		hst->dt[tTable->tdt[i].idt] = tTable->tdt[i].dt;
	}
	
	for(int i=0; i<nDt; i++){
		hst->counts[i] = malloc(sizeof(int)*nBins);
		hst->avgCor[i] = malloc(sizeof(double)*nBins);
		for(int bin=0; bin<nBins; bin++){
			hst->counts[i][bin]=0;
			hst->avgCor[i][bin]=0;
		}
	}
	return hst;
}

PartitionTable* PartitionNew(SimProperties* sp, int boxL, int nMaxPol){
	PartitionTable* pt = malloc(sizeof(PartitionTable));
	pt->xBox = LT/boxL;
	pt->boxL = boxL;
	pt->nMaxPol = nMaxPol;
	pt->nBox = pt->xBox*pt->xBox*pt->xBox;
	
	pt->nPol = malloc(sizeof(int)*pt->nBox);
	pt->polList = malloc(sizeof(int*)*pt->nBox);
	for(int i=0; i<pt->nBox; i++) pt->polList[i] = malloc(sizeof(int)*pt->nMaxPol);
	return pt;
}

double** LoadAvgPos(SimProperties* sp){
	char file[10000];
	double** avgPos;
	avgPos = malloc(sizeof(double*)*sp->nTime);
// 	for(int i=0; i<sp->nTime; i++) avgPos = malloc(sizeof
	
	sprintf(file, "%s/avg_pos.dat", sp->sampleDir);
	FILE* pFile = fopen(file, "r");
	for(int i=0; i<sp->nTime; i++){
		avgPos[i] = malloc(sizeof(double)*3);
		fscanf(pFile, "%lf %lf %lf\n", avgPos[i], avgPos[i]+1, avgPos[i]+2);
	}
	fclose(pFile);
	return avgPos;
}

void CreatePartition(PartitionTable* pt, CMSData* cms, SimProperties* sp){
	double t,u,v;
	int tb, ub, vb;
	
	for(int iBox=0; iBox<pt->nBox; iBox++){
		pt->nPol[iBox]=0;
	}
	
	for(int iPol=0; iPol<sp->nPol; iPol++){
		t=cms->TUVRemPer[iPol][0];
		u=cms->TUVRemPer[iPol][1];
		v=cms->TUVRemPer[iPol][2];
		
		tb = (int) (t/pt->boxL);
		ub = (int) (u/pt->boxL);
		vb = (int) (v/pt->boxL);
		
		int box = tb + ub*pt->xBox + vb*pt->xBox*pt->xBox;
		if(pt->nPol[box]<pt->nMaxPol)
			pt->polList[box][pt->nPol[box]++] = iPol;
	}
}

void ComputeCorrelations(SimProperties* sp){
	int boxL = MIN(LT/2, (int)pow(150*sp->maxNMono, 1./3.));
	while((LT/boxL)*boxL != LT) boxL++;
	printf("Using boxSize of %i\n", boxL);
	
	
	CMSData *cms1 = CMSDataNew(sp);
	CMSData *cms2 = CMSDataNew(sp);
	PartitionTable *pt = PartitionNew(sp, boxL, 1000);
	double** avgPos = LoadAvgPos(sp);
	TDTTable* tTable = TTableNew(sp, 1);
	Histogram* hst = HistogramNew(sp, boxL*boxL, 1, tTable);
	
	double** diff = malloc(sizeof(double*)*sp->nPol);
	for(int i=0; i<sp->nPol; i++) diff[i] = malloc(sizeof(double)*3);
	
	int curT=-1, curT2=-1;
	
	for(int i=0; i<tTable->nTDT; i++){
		if(tTable->tdt[i].t != curT){
			LoadCMS(tTable->tdt[i].t, cms1, sp);
			CreatePartition(pt, cms1, sp);
		}
		curT = tTable->tdt[i].t;
		curT2 = curT+tTable->tdt[i].dt;
		LoadCMS(curT2, cms2, sp);
		CalcDiff(cms1, cms2, curT, curT2, avgPos, diff, sp);
// 		printf("t=(%i,%i)\n", curT, curT2);
		ProcessPartition(pt, cms1, diff, sp, hst, tTable->tdt[i].idt);
		printf("\b \b\b \b\b \rMeasure: %i/%i [%.1lf%%]", i+1, tTable->nTDT, 100*(i+1)/(double)tTable->nTDT);
		fflush(NULL);
	}
	printf("\n");
	WriteHistogram(hst, sp);
}

void CalcDiff(CMSData* cms, CMSData* cmsDest, long t1, long t2, double** avgPos, double** diff, SimProperties* sp){
	for(int i=0; i<sp->nPol; i++){
		for(int k=0; k<3; k++){
			diff[i][k] = cmsDest->XYZ[i][k]-cms->XYZ[i][k] - (avgPos[t2][k]-avgPos[t1][k]);
		}
	}
}

void AddToHistogram(Histogram* hst, double drInit, long iDt, double drDiff){
	int bin = drInit/hst->binSize;
	if(bin >= hst->nBins){
// 		printf("Warning: not enough bins\n");
		bin = hst->nBins-1;
	}
	hst->avgCor[iDt][bin] += drDiff;
	hst->counts[iDt][bin]++;
}

void ProcessPartition(PartitionTable* pt, CMSData* cms, double** diff, SimProperties* sp, Histogram* hst, int iDt){
	for(int iBox=0; iBox<pt->nBox; iBox++){
		for(int iPol=0; iPol<pt->nPol[iBox]; iPol++){
			int polIdI = pt->polList[iBox][iPol];
			for(int jPol=iPol+1; jPol<pt->nPol[iBox]; jPol++){
				int polIdJ = pt->polList[iBox][jPol];
				double drInit = DRXYZ(cms->XYZRemPer[polIdI], cms->XYZRemPer[polIdJ]);
				
				double drDiff=0;
				for(int k=0; k<3; k++){
					drDiff += diff[iPol][k]*diff[jPol][k];
				}
// 				if(iDt == 0 && drInit < 1){
// 					printf("\n%lf %lf %lf  ", diff[iPol][0], diff[iPol][1], diff[iPol][2]);
// 					printf("%lf %lf %lf [drInit=%lf, cor=%lf]\n", diff[jPol][0], diff[jPol][1], diff[jPol][2], drInit, drDiff);
// // 					printf("(%lf %lf %lf) vs (%lf %lf %lf)\n", cms->TUV[polIdI][0], cms->TUV[polIdI][1], cms->TUV[polIdI][2], cms->TUV[polIdJ][0],cms->TUV[polIdJ][1], cms->TUV[polIdJ][2]);
// // 					printf("(%lf %lf %lf) vs (%lf %lf %lf)\n", cms->XYZ[polIdI][0], cms->XYZ[polIdI][1], cms->XYZ[polIdI][2], cms->XYZ[polIdJ][0],cms->XYZ[polIdJ][1], cms->XYZ[polIdJ][2]);
// 					
// 				}
				AddToHistogram(hst, drInit, iDt, drDiff);
			}
		}
	}
}

void WriteHistogram(Histogram* hst, SimProperties* sp){
	char file[10000];
	sprintf(file, "%s/cms_cor.dat", sp->sampleDir);
	int maxBin, exit;
	for(maxBin=hst->nBins-1, exit=0; maxBin>0 && !exit; maxBin--){
		for(int iDt=0; iDt<hst->nTime; iDt++){
			if(hst->counts[iDt][maxBin]){
				exit=1;
				break;
			}
		}
	}
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "# ");
	for(int iDt=0; iDt<hst->nTime; iDt++)
		fprintf(pFile, "%li ", hst->dt[iDt]*sp->dT);
	fprintf(pFile, "\n");
	
	for(int bin=0; bin<maxBin; bin++){
		fprintf(pFile, "%lf ", (bin+0.5)*hst->binSize);
		for(int iDt=0; iDt<hst->nTime; iDt++){
			if(hst->counts[iDt][bin])
				fprintf(pFile, "%lf ", hst->avgCor[iDt][bin]/hst->counts[iDt][bin]);
			else
				fprintf(pFile, "nan ");
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}
