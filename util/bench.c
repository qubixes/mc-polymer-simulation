#include "bench.h"

int compare (const void * a, const void * b)
{
	double da = *(double*)a; double db = *(double*)b;
	if(da>db) return 1;
	if(da==db) return 0;
	if(da<db) return -1;
	return 0;
}

int main(int argc, char** argv){
	char* simType;
	char exec[10000];
	char baseExec[10000];
	int nThreads=8;
	long baseTime=10;
	int nMeas=10;
	int nRemove=2;
	long interMax=10;
	
	char dorunDir[]="../";
	char dataDir[]="../testdata";
	Timer t;
	
	if(argc <=1) {
		printf("Need at least one argument: sim type!\n");
		return 193;
	}
	
	simType = argv[1];
	if(argc >= 3){
		nThreads=atoi(argv[2]);
	}
	
	if(!strcmp(simType, "--opencl") || !strcmp(simType, "--cuda")) baseTime *=20;
	
	//First increase MC steps to find time per MC.
	
	sprintf(baseExec, "%s/do_run.sh %s --dir %s --nmono 1000 --density 1.0 --threads %i", dorunDir, simType, dataDir, nThreads);
	long interval, nTime;
	
	double **timings;
	
	timings = malloc(sizeof(double*)*2);
	for(int i=0; i<2; i++){
		timings[i] = malloc(sizeof(double)*nMeas);
	}
	
	
	CleanDir(dataDir);
	double avgDif;
	
	for(int i=0; i<nMeas; i++){
		for(int j=1; j<=2; j++){
			interval = j*baseTime;
			nTime = interval;
			sprintf(exec, "%s --time %li --interval %li > /dev/null", baseExec, nTime, interval);
			TimerStart(&t);
			system(exec);
			timings[j-1][i] = TimerElapsed(&t);
			CleanDir(dataDir);
		}
	}
	
	avgDif = ComputeAverageDif(timings, nRemove, nMeas);
	
	double MCTime=1e9*avgDif/(double)(LAT_SIZE*baseTime);
	double MCSpeed=1e-6*(LAT_SIZE*baseTime)/avgDif;
	
	avgDif=0;
	for(int i=0; i<nMeas; i++){
		for(int j=0; j<2; j++){
			interval = ((j==1)?1:interMax)*(baseTime/interMax);
			nTime = interMax*(baseTime/interMax);
			sprintf(exec, "%s --time %li --interval %li > /dev/null", baseExec, nTime, interval);
			TimerStart(&t);
			system(exec);
			timings[j][i] = TimerElapsed(&t);
			CleanDir(dataDir);
		}
	}
	avgDif = ComputeAverageDif(timings, nRemove, nMeas);
	
	double IOTime=1e3*avgDif/(double)(interMax-1);
	double IOSpeed=(interMax-1)/avgDif;
	
	printf("%s %i %lf %lf %lf %lf\n", simType, nThreads, MCTime, MCSpeed, IOTime, IOSpeed);
}

void CleanDir(char* dir){
	char exec[10000];
	
	sprintf(exec, "X=`ls %s/long`; for dir in $X; do if [ $dir != \"t=0_dev=0.res\" ]; then rm %s/long/$dir; fi; done", dir, dir);
	system(exec);
}

double ComputeAverageDif(double** timings, int nRemove, int nMeas){
	double avgTimings[2]={0,0};
	
	for(int j=0; j<2; j++){
		qsort(timings[j], nMeas, sizeof(double), compare);
		for(int i=0; i< nMeas-nRemove; i++){
			avgTimings[j] += timings[j][i];
		}
		avgTimings[j] /= nMeas-nRemove;
	}
	
	return (avgTimings[1]-avgTimings[0]);
}
	