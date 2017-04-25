#include <stdlib.h>
#include <stdio.h>

typedef struct Data{
	int nTime;
	int nAlloc;
	double* t;
	double* val;
}Data;

typedef struct SimProperties{
	Data datIn;
	Data datOut;
	int onlyPositive;
	double tFac;
	char* file;
}SimProperties;

SimProperties sp;

void ReadData(Data* dat, char* file);
void SmoothData(SimProperties* sp);
void PrintData(Data* dat, int onlyPositive);

int main(int argc, char** argv){
	sp.onlyPositive = 1;
	if(argc<2){
		printf("Need file!\n");
		exit(0);
	}
	if(argc>=3){
		sp.onlyPositive = atoi(argv[2]);
	}
	
	sp.file = argv[1];
	sp.tFac = 1.2;
	
	ReadData(&(sp.datIn), sp.file);
	SmoothData(&sp);
	PrintData(&(sp.datOut), sp.onlyPositive);
	return 0;
}

void ReadData(Data* dat, char* file){
	int allocStep=1000;
	dat->nAlloc = allocStep;
	dat->t = malloc(sizeof(long)*dat->nAlloc);
	dat->val = malloc(sizeof(double)*dat->nAlloc);
	dat->nTime=0;
	
	FILE* pFile = fopen(file, "r");
	
	long t; double val;
	while(fscanf(pFile, "%li %lf\n", &t, &val) > 0){
		if(dat->nTime == dat->nAlloc){
			dat->nAlloc += allocStep;
			dat->t = realloc(dat->t, dat->nAlloc*sizeof(double));
			dat->val = realloc(dat->val, dat->nAlloc*sizeof(double));
		}
		dat->t[dat->nTime] = (double)t;
		dat->val[dat->nTime++] = val;
	}
	fclose(pFile);
}

void PrintData(Data* dat, int onlyPositive){
	for(int i=0; i<dat->nTime; i++){
		if(dat->val[i] < 0 && onlyPositive) return;
		printf("%li %lf\n", (long) dat->t[i], dat->val[i]);
	}
}

void SmoothData(SimProperties* sp){
	Data* dat = &sp->datIn;
	Data* datOut = &sp->datOut;
	
	datOut->nTime = dat->nTime;
	datOut->t = malloc(dat->nTime*sizeof(long));
	datOut->val = malloc(dat->nTime*sizeof(double));
	
	for(int iTime=0; iTime<dat->nTime; iTime++){
		double t1 = dat->t[iTime];
		double avgT=0, avgVal=0; int num=0;
		for(int jTime= iTime; jTime <dat->nTime; jTime++){
			double t2 = dat->t[jTime];
			
			if(t1*sp->tFac < t2) break;
			avgT += t2; avgVal += dat->val[jTime]; num++;
		}
		datOut->t[iTime] = avgT/num;
		datOut->val[iTime] = avgVal/num;
	}
}
