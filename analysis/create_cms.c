#include "lowm_modes.h"
#include "file.h"
#include "timer.h"
#include <sys/resource.h>

void ConvertData(SimProperties* sp, int devId);

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
	int update=NeedsUpdatePath("simulation_settings.txt", "cms/t=0_dev=0.res", sampleDir);
	if(!update){return 0;}
// 	update=!update;
	SetSimProps(&sp, sampleDir);
	for(int iDev=0; iDev<sp.nDev; iDev++){
		ConvertData(&sp, iDev);
	}
	printf("\n\n");
	return 0;
}



void ConvertData(SimProperties* sp, int devId){
	char file[10000], fileOut[10000];
	char exec[10000];
	int curCfg=0;
	
	Timer timer;
	
// 	PolyConfig* pcfg[2];
	char* str = malloc(sizeof(char)*(sp->polSize+2));
	double* t[2], *u[2], *v[2];
	for(int i=0; i<2; i++){
		t[i] = malloc(sizeof(double)*sp->nPol);
		u[i] = malloc(sizeof(double)*sp->nPol);
		v[i] = malloc(sizeof(double)*sp->nPol);
	}
	
	
	
	sprintf(exec, "mkdir -p %s/cms/", sp->sampleDir);
	
	for(int iTime=0; iTime<sp->nTime; iTime++){
		TimerStart(&timer);
		sprintf(file, "%s/t=%li_dev=%i.res", sp->sampleDir, iTime*sp->dT, devId);
		sprintf(fileOut, "%s/cms/t=%li_dev=%i.res", sp->sampleDir, iTime*sp->dT, devId);
		FILE* pFile = fopen(file, "r");
		FILE* pFileOut = fopen(fileOut, "w");
		if(!pFile) printf("woops: not opening file %s\n", file);
		if(!pFileOut) printf("Not opening file for writing: %s\n", fileOut);
		
		for(int i=0; i<5; i++) fscanf(pFile, "%*s %*s");
		
		for(int iPol=0; iPol<sp->nPol; iPol++){
			fscanf(pFile, "%*s %*s");
			fscanf(pFile, "%lf %lf %lf", &t[curCfg][iPol], &u[curCfg][iPol], &v[curCfg][iPol]);
			fscanf(pFile, "%s", str);
			
			if(iTime != 0){
				while(t[curCfg][iPol] - t[curCfg^1][iPol] >  LT/2) t[curCfg][iPol] -= LT;
				while(t[curCfg][iPol] - t[curCfg^1][iPol] < -LT/2) t[curCfg][iPol] += LT;
				while(u[curCfg][iPol] - u[curCfg^1][iPol] >  LU/2) u[curCfg][iPol] -= LU;
				while(u[curCfg][iPol] - u[curCfg^1][iPol] < -LU/2) u[curCfg][iPol] += LU;
				while(v[curCfg][iPol] - v[curCfg^1][iPol] >  LV/2) v[curCfg][iPol] -= LV;
				while(v[curCfg][iPol] - v[curCfg^1][iPol] < -LV/2) v[curCfg][iPol] += LV;
			}
			double curT=t[curCfg][iPol];
			double curU=u[curCfg][iPol];
			double curV=v[curCfg][iPol];
			
			Coor cms={0,0,0};
			
			for(int iMono=0; iMono<sp->polSize; iMono++){
// 				cms.x += curT+curU-curV;
// 				cms.y += curT-curU;
// 				cms.z += curV;
				cms.x += curT;
				cms.y += curU;
				cms.z += curV;
				int step = CharToHex(str[iMono]);
				
				curT += ((step>>0)&0x1) - ((step>>3)&0x1);
				curU += ((step>>1)&0x1) - ((step>>3)&0x1);
				curV += ((step>>2)&0x1) - ((step>>3)&0x1);
			}
			cms.x /= sp->polSize;
			cms.y /= sp->polSize;
			cms.z /= sp->polSize;
			
			fprintf(pFileOut, "%lf %lf %lf\n", cms.x, cms.y, cms.z);
		}
		printf("\b \b\b \b\b \b\rTime %i/%i [%.1lf%%, %.1lf ms]", iTime, sp->nTime, (iTime+1)/(double)sp->nTime*100, TimerElapsed(&timer)*1e3);
		fflush(NULL);
		curCfg ^= 1;
		fclose(pFile);
		fclose(pFileOut);
	}	
}
