#include <stdlib.h>
#include <stdio.h>
#include "lowm_modes.h"

typedef struct Data{
	double *x, *y, *z;
	char* str;
	int N;
}Data;

typedef struct RunProperties{
	char* dir;
	char* binDir;
	char fileIn[1000];
	char dirOut[1000];
	int polId;
	long t;
}RunProperties;

Data* ReadData(char* file, int polId);
void RunSecAnalysis(Data* data, RunProperties* rp);

int main(int argc, char** argv){
	if(argc<3){
		printf("Need two arguments!\n");
		return 192;
	}
	RunProperties rp;
	
	rp.dir = argv[1];
	rp.binDir = argv[2];
	rp.t = atol(argv[3]);
	rp.polId = atoi(argv[4]);
	
	sprintf(rp.fileIn, "%s/t=%li_dev=0.res", rp.dir, rp.t);
	sprintf(rp.dirOut, "%s/sec_struct", rp.dir);
	
	Data* data = ReadData(rp.fileIn, rp.polId);
	RunSecAnalysis(data, &rp);
}

Data* NewData(int N){
	Data* data = malloc(sizeof(Data));
	data->x = malloc(sizeof(double)*N);
	data->y = malloc(sizeof(double)*N);
	data->z = malloc(sizeof(double)*N);
	data->str = malloc(sizeof(char)*(N+1));
	data->N = N;
	return data;
}

Data* ReadData(char* file, int polId){
	int nPol, maxPolSize;
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int i=0; i<3; i++)
		fscanf(pFile, "%*s %*s");
	
	fscanf(pFile, "%*s %i", &nPol);
	if(polId >= nPol){
		printf("Error: polId does not exist (too large)\n");
		exit(193);
	}
	
	fscanf(pFile, "%*s %i", &maxPolSize);
// 	printf("maxPolSize = %i\n", maxPolSize);
	Data* data = NewData(maxPolSize);
	
	for(int i=0; i<polId; i++)
		fscanf(pFile, "%*s %*s %*s %*s %*s %*s");
	
	fscanf(pFile, "%*s %i %*s %*s %*s %s", &data->N, data->str);
	fclose(pFile);
	
// 	printf("%s\n", data->str);
	
	int t=0,u=0,v=0;
	for(int i=0; i<data->N; i++){
		data->x[i] = (t+u-v)/sqrt(2);
		data->y[i] = (t-u)/sqrt(2);
		data->z[i] = v/sqrt(2);
// 		printf("%lf %lf %lf\n", data->x[i], data->y[i], data->z[i]);
		
		int bond = CharToHex(data->str[i]);
		int w = bond>>3;
		t += ((bond>>0)&0x1) - w;
		u += ((bond>>1)&0x1) - w;
		v += ((bond>>2)&0x1) - w;
	}
	return data;
}

Data* CoarseGrain(Data* data, int maxN){
	if(data->N <= maxN) return data;
	
	Data* newData = NewData(maxN);
	
	double frac = data->N/(double)newData->N;
	
	double leftOver=1;
	int curMono=0;
	for(int i=0; i<newData->N; i++){
		newData->x[i]=0; newData->y[i]=0; newData->z[i]=0;
		
		newData->x[i] += leftOver*data->x[curMono];
		newData->y[i] += leftOver*data->y[curMono];
		newData->z[i] += leftOver*data->z[curMono];
		
		int nWhole = (int)(frac-leftOver);
		double newFrac = frac-nWhole-leftOver;
		curMono++;
		for(int j=0; j<nWhole; j++, curMono++){
			newData->x[i] += data->x[curMono];
			newData->y[i] += data->y[curMono];
			newData->z[i] += data->z[curMono];
		}
		
		if(newFrac > 1e-3){
			newData->x[i] += newFrac*data->x[curMono];
			newData->y[i] += newFrac*data->y[curMono];
			newData->z[i] += newFrac*data->z[curMono];
		}
		
		newData->x[i] /= frac;
		newData->y[i] /= frac;
		newData->z[i] /= frac;
		
		leftOver = 1-newFrac;
	}
	return newData;
}

void RunSecAnalysis(Data* data, RunProperties* rp){
	char exec[1000];
// 	sprintf(exec, "sed 's/RING_LENGTH_SUB/%i/' %s/input_template.dat > %s/input.dat", data->N);
// 	system(exec);
	
	data = CoarseGrain(data, 10000);
	
	sprintf(exec, "%s/dJost.exe %s %i\n", rp->binDir, rp->dirOut, rp->polId);
	FILE* pFile = popen(exec, "w");
	for(int i=0; i<data->N; i++){
		fprintf(pFile, "%lf %lf %lf\n", data->x[i], data->y[i], data->z[i]);
	}
	pclose(pFile);
}
