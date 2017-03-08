#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define PI 3.14159265359

typedef struct Data{
	int N;
	double** pc;
}Data;

Data* ReadData(char* file);
void PrintData(Data* pData);
void Print2dData(Data* pData);
Data* ConvertPCtoRSQ(Data* pcData);
Data* NewData(int N);
void PrintAgregate(Data* data);

int main(int argc, char** argv){
	if(argc<2){
		printf("Need file to process\n");
		return 192;
	}
	
	Data* pData = ReadData(argv[1]);
	Data* rData = ConvertPCtoRSQ(pData);
	PrintAgregate(rData);
	return 0;
}

Data* NewData(int N){
	Data* data = malloc(sizeof(Data));
	data->N = N;
	data->pc = malloc(sizeof(double*)*data->N);
	for(int i=0; i<data->N; i++)
		data->pc[i] = malloc(sizeof(double)*data->N);
	return data;
}


Data* ReadData(char* file){
	char exec[1000];
	int N;
	int i,j;
	double p;

	
	FILE* pFile = fopen(file, "r");
	if(!pFile) return NULL;
	sprintf(exec, "wc -l %s", file);
	FILE* pPipe = popen(exec, "r");
	fscanf(pPipe, "%i", &N);
	N = (int)(sqrt(N));
	Data* pData = NewData(N);
	while(fscanf(pFile, "%i %i %lf", &i, &j, &p)>0){
		pData->pc[i][j] = p;
	}
	fclose(pPipe);
	fclose(pFile);
	return pData;
}

void PrintData(Data* pData){
	for(int i=0; i<pData->N; i++){
		for(int j=0; j<pData->N; j++){
			printf("%le ", pData->pc[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void Print2dData(Data* pData){
	for(int i=0; i<pData->N; i++){
		for(int j=0; j<pData->N; j++){
			printf("%i %i %le\n", i, j, pData->pc[i][j]);
		}
	}
	printf("\n");
}

void PrintAgregate(Data* data){
	double* agregate = malloc(sizeof(double)*data->N);
	memset(agregate, 0, sizeof(double)*data->N);
	
	
	for(int i=0; i<data->N; i++){
		for(int j=i; j<data->N; j++){
			agregate[j-i] += data->pc[i][j];
		}
	}
	
	for(int i=0; i<data->N; i++)
		printf("%i %lf\n", i, agregate[i]/(double)(data->N-i));
	printf("\n");
}


Data* ConvertPCtoRSQ(Data* pcData){
	Data* rData = NewData(pcData->N);
	for(int i=0; i<pcData->N; i++){
		for(int j=0; j<pcData->N; j++){
			if(pcData->pc[i][j])
				rData->pc[i][j] = 1/(2*PI*pow(pcData->pc[i][j], 2./3.));
			else
				rData->pc[i][j] = -1;
		}
	}
// 	Print2dData(rData);
	return rData;
}