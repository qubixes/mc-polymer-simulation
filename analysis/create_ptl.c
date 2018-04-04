#include "lowm_modes.h"
#include "file.h"
#include <sys/resource.h>

typedef struct FileHandles{
	fpos_t* timeFPos;
	char** timeFNames;
	FILE** pFilePol;
}FileHandles;
void ConvertData(SimProperties* sp, int devId, int nPolOpen);

int main(int argc, char** argv){
	SimProperties sp;
	char* sampleDir;
	char exec[10000];
	int nPolOpen;
	struct rlimit limit;
	
	
	for(nPolOpen=1000; nPolOpen>=1; nPolOpen /=2){
		limit.rlim_cur = nPolOpen+50;
		limit.rlim_max = nPolOpen+50;
		if(setrlimit(RLIMIT_NOFILE, &limit) == 0) break;
	}
	
	
	if(argc<2){ 
		printf("Need basedir\n");
		exit(0);
	}
	sampleDir = argv[1];
	sp.neFile=NULL;
	
	sprintf(exec, "mkdir -p %s/ptl", sampleDir); system(exec);
	int update=NeedsUpdatePath("simulation_settings.txt", "ptl/pol=0_dev=0.res", sampleDir);
	if(!update){return 0;}
// 	update=!update;
	SetSimProps(&sp, sampleDir);
	for(int iDev=0; iDev<sp.nDev; iDev++){
		ConvertData(&sp, iDev, nPolOpen);
	}
	printf("\n\n");
	return 0;
}

void ConvertData(SimProperties* sp, int devId, int nPolOpen){
	int startPol=0;
	
	char file[10000];
	char polFile[1000];
	char* tString = malloc(sizeof(char)*(sp->maxNMono+1));
	FileHandles fh;
	FILE* timePFile;
	Timer timer;
	long totOps=sp->nTime*sp->nPol*sp->nDev, curOp;
	
	fh.timeFPos = malloc(sizeof(fpos_t)*sp->nTime);
	fh.timeFNames = malloc(sizeof(char*)*sp->nTime);
	fh.pFilePol = malloc(sizeof(FILE*)*nPolOpen);
	int t,u,v;

	for(int i=0; i<sp->nTime; i++){
		sprintf(file, "%s/t=%li_dev=%i.res", sp->sampleDir, i*sp->dT, devId);
		fh.timeFNames[i] = malloc(sizeof(char)*strlen(file)+1);
		strcpy(fh.timeFNames[i], file);
		FILE* pFile = fopen(file, "r");
		if(!pFile){
			printf("woops: not opening file %s (%i)\n", file, sp->nTime);
			exit(192);
		}
		
		for(int i=0; i<5; i++) fscanf(pFile, "%*s %*s");
		
		fgetpos(pFile, &(fh.timeFPos[i]));
		fclose(pFile);
	}
	while(startPol<sp->nPol){
		int nPolRead = MIN(nPolOpen, sp->nPol-startPol);
		for(int iPol=startPol; iPol<sp->nPol && iPol-startPol<nPolOpen; iPol++){
			sprintf(polFile, "%s/ptl/pol=%i_dev=%i.res", sp->sampleDir, iPol, devId);
			fh.pFilePol[iPol-startPol] = fopen(polFile, "w");
			if(!fh.pFilePol[iPol-startPol]) printf("Error opening file %s!\n", polFile);
		}
		
		for(int iTime=0; iTime<sp->nTime; iTime++){
			TimerStart(&timer);
			timePFile = fopen(fh.timeFNames[iTime], "r");
			fsetpos(timePFile, &(fh.timeFPos[iTime]));
			for(int jPol=0; jPol<nPolRead; jPol++){
				int polSize;
				fscanf(timePFile, "%*s %i", &polSize);
				if(fscanf(timePFile, "%i %i %i %s", &t, &u, &v, tString)<=0){
					printf("\nError reading file at file id %i\n", iTime);
					exit(192);
				}
				if(iTime == 0)
					fprintf(fh.pFilePol[jPol], "polId= %i\nnTime= %i\nlen= %i\n", jPol, sp->nTime, polSize);
				fprintf(fh.pFilePol[jPol], "%i %i %i %s\n", t,u,v,tString);
			}
			fgetpos(timePFile, &(fh.timeFPos[iTime]));
			fclose(timePFile);
			if(iTime%100 == 0 || iTime == sp->nTime-1){
				curOp = (long)devId*sp->nTime*sp->nPol+(long)startPol*sp->nTime+(long)iTime*nPolRead+nPolRead;
				printf("\b \b \b \b \rIO step: %.2f ms [%li/%li ->  %.1f%%]", TimerElapsed(&timer), curOp, totOps, 100*curOp/(double)totOps);
			}
			fflush(NULL);
		}
		
		for(int iPol=0; iPol<nPolRead; iPol++){
			fclose(fh.pFilePol[iPol]);
		}
		startPol += nPolOpen;
	}
// 	printf("\n\n");
}

long GetDT(char* sampleDir, char** firstFile, long* firstT){
	char exec[10000];
	char* file = malloc(10000*sizeof(char));
	char* retFile = malloc(10000*sizeof(char));
	FILE* pPipe;
	int i;
	long t;
	long minT=(long)1e13, secMinT=(long)1e14;
	
	*firstFile = retFile;
	sprintf(exec, "ls %s | grep 't=' |grep 'dev=0'", sampleDir);
	pPipe = popen(exec, "r");
	printf("Command: %s\n", exec);
	if(!pPipe){
		printf("error opening pipe\n");
		exit(0);
	}
	while(fscanf(pPipe, "%s", file)>0){
		i=0; 
		while(file[i] != '_' && file[i] != '\0') i++;
		if(file[i] != '_') continue;
		file[i]='\0';
		t = atoi(file+2);
		file[i]='_';
		if(t<minT){
			secMinT = minT;
			minT = t;
			strcpy(retFile, file);
			*firstT = t;
		}
		else if(t<secMinT){
			secMinT = t;
		}
	}
	pclose(pPipe);
	free(file);
	return secMinT-minT;
}


void SetSimProps(SimProperties* sp, char* sampleDir){
	FILE *pFile;
	char filename[10000];
	char* firstFile;
	char polType[100];
	char exec[1000];
	int nFiles;
	
	sp->sampleDir = sampleDir;
	printf("Sample directory: %s\n", sampleDir);
	sp->dT = GetDT(sampleDir, &firstFile, &(sp->tStart));
	sprintf(filename, "%s/%s", sampleDir, firstFile);
	printf("file = %s\n", filename);
	printf("dT = %li\n", sp->dT);
	pFile=fopen(filename,"r");
	fscanf(pFile, "%*s %i", &sp->LT);
	fscanf(pFile, "%*s %i", &sp->LU);
	fscanf(pFile, "%*s %i", &sp->LV);
	
// 	if(sp->LT > MAX_LT || sp->LU > MAX_LU || sp->LV > MAX_LV){
// 		printf("Box is too large: (%i, %i, %i) vs (%i,%i,%i)\n", sp->LT, sp->LU, sp->LV, MAX_LT, MAX_LU, MAX_LV);
// 		exit(0);
// 	}
	
	fscanf(pFile, "%*s %i", &sp->nPol);
	fscanf(pFile, "%*s %li", &sp->maxNMono);
	printf("nPol = %i\n", sp->nPol);
	printf("polSize=%li\n", sp->maxNMono);
	fclose(pFile);
	
	sprintf(exec, "ls %s | grep 't=' | grep 'dev=0' | wc -w", sampleDir);
	pFile = popen(exec, "r");
	fscanf(pFile, "%i", &sp->nTime);
	printf("Number of files detected: %i\n", sp->nTime);
	pclose(pFile);
	
	sprintf(exec, "ls %s | grep 't=' | wc -w", sampleDir);
	pFile = popen(exec, "r");
	fscanf(pFile, "%i", &nFiles);
	if(nFiles%sp->nTime != 0){
		printf("Error: uneven number of files for the different devices\n");
		exit(123);
	}
	sp->nDev = nFiles/sp->nTime;
	printf("Number of devices detected: %i\n", sp->nDev);
	pclose(pFile);
	
	sprintf(filename, "%s/simulation_settings.txt", sampleDir);
	sprintf(exec, "grep 'Polytype' %s", filename);
	pFile = popen(exec, "r");
	fscanf(pFile, "%*s %*s %s", polType);
	pclose(pFile);
	if(!strcmp(polType, "ring")){
		sp->polTypeMelt = POL_TYPE_RING;
	}
	else if (!strcmp(polType, "lin")){
		sp->polTypeMelt = POL_TYPE_LIN;
	}
	else if (!strcmp(polType, "mixed"))
		sp->polTypeMelt = POL_TYPE_MIXED;
	else{
		printf("Error: unknown polymer type %s\n", polType);
		exit(192);
	}
}


