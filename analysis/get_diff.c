#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


typedef struct Data{
	long* t;
	double* rSq;
	double* D;
	int nTime;
	int nAlloc;
}Data;

void ReadData(Data* dat, char* file);
void GetDiff(Data* dat);

int main(int argc, char** argv){
	Data dat;
	char* file;
	if(argc<2) {
		printf("Error: need a file as input\n");
		return 192;
	}
	
	file = argv[1];
	
	ReadData(&dat, file);
	GetDiff(&dat);
	return 0;
}


void ReadData(Data* dat, char* file){
	int allocStep=1000;
	dat->nAlloc = allocStep;
	dat->t = malloc(sizeof(long)*dat->nAlloc);
	dat->rSq = malloc(sizeof(double)*dat->nAlloc);
	dat->D = malloc(sizeof(double)*dat->nAlloc);
	dat->nTime=0;
	
	FILE* pFile = fopen(file, "r");
	
	long t; double val;
	while(fscanf(pFile, "%li %lf\n", &t, &val) > 0){
		if(dat->nTime == dat->nAlloc){
			dat->nAlloc += allocStep;
			dat->t = realloc(dat->t, dat->nAlloc*sizeof(double));
			dat->rSq = realloc(dat->rSq, dat->nAlloc*sizeof(double));
			dat->D = realloc(dat->D, dat->nAlloc*sizeof(double));
		}
		dat->t[dat->nTime] = (double)t;
		dat->rSq[dat->nTime] = val;
		dat->D[dat->nTime++] = val/t;
	}
// 	printf("nTime = %i\n", dat->nTime);
	fclose(pFile);
}

void GetDiff(Data* dat){
	double DMin=1e301;
	double DMax=1032;
	int iStart=-1;
	
	for(int i=0; i<dat->nTime; i++){
		double D = dat->D[i];
		if(D<DMin) DMin = D;
	}
	
	for(int i=0; i<dat->nTime; i++){
		double D = dat->D[i];
		if(D<DMin*1.15){
			DMax = D;
			iStart = i;
			break;
		}
	}
// 	printf("DMax, DMin: %le, %le\n", DMax, DMin);
// 	double curD=0;
	int nThrow = (int)((dat->nTime-iStart-1)*0.1);
	nThrow = (nThrow>0)?nThrow:1;
	
	double curMinD;
	int iMinD;
	
	for(int j=0; j<nThrow; j++){
		curMinD = 1e301;
		iMinD = -1;
		for(int i=iStart; i<dat->nTime; i++){
			if(dat->D[i] < curMinD){
				curMinD = dat->D[i];
				iMinD = i;
			}
		}
		dat->D[iMinD] = 1e301;
	}
	int j=1;
	double aBest=1, bBest=0; 
	
	for(int i=0; i<iStart; i++){
		if(dat->t[i]*50 < dat->t[iStart]) continue;
		while(dat->t[j] < dat->t[i]*2) j++;
		if(j>= iStart) break;
		double a = log(dat->D[j]/dat->D[i])/log(dat->t[j]/(double)dat->t[i]);
// 		printf("%li %lf\n", dat->t[i], a);
		if(a < aBest){
			aBest = a;
			bBest = log(dat->D[i])-(double)(aBest*log(dat->t[i]));
		}
	}
	
	double tD = exp( (log(curMinD)-bBest)/aBest);
	double rSqD = tD*curMinD;
	
	printf("%le %le %le\n", curMinD/6.0, tD, rSqD);
}
