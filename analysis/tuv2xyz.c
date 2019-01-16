#include <stdlib.h>
#include <stdio.h>
#include "lowm_modes.h"

typedef struct Data{
	double *x, *y, *z;
	char* str;
	int N;
}Data;

typedef struct RunProperties{
	char* fileIn;
	char* fileOut;
	int polId;
}RunProperties;

Data* ReadData(char* file, int polId);
void WriteXYZData(Data* data, char* file);

int main(int argc, char** argv){
	if(argc<4){
		printf("Need two arguments!\n");
		return 192;
	}
	RunProperties rp;
	
	rp.fileIn = argv[1];
	rp.fileOut = argv[2];
	rp.polId = atoi(argv[3]);
	
	Data* data = ReadData(rp.fileIn, rp.polId);
	WriteXYZData(data, rp.fileOut);
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



void WriteXYZData(Data* data, char* file){
	FILE* pFile = fopen(file, "w");
	for(int i=0; i<data->N; i++){
		fprintf(pFile, "%lf %lf %lf\n", data->x[i], data->y[i], data->z[i]);
	}
	pclose(pFile);
}
