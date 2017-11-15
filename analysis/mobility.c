#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(X,Y) ((X>Y)?X:Y)

void AnalyzeFile(char* file);

int main(int argc, char** argv){
	
	if(argc<2){
		printf("Need file to analyze\n");
		return 192;
	}
	
	AnalyzeFile(argv[1]);
	
	return 0;
}

double CompareString(char* str1, char* str2, int len){
	int maxSame=0;
	
	char* nzStr1 = malloc(sizeof(char)*(len+1));
	char* nzStr2 = malloc(sizeof(char)*(len+1));
	
	int nzLen1=0;
	int nzLen2=0;
	
	for(int iStr=0; iStr<len; iStr++){
		if(str1[iStr] != '0') nzStr1[nzLen1++] = str1[iStr];
		if(str2[iStr] != '0') nzStr2[nzLen2++] = str2[iStr];
	}
	
	for(int start=0; start<nzLen1; start++){
		int same=0;
		for(int i=0; i<nzLen2; i++){
			if(nzStr2[i] == nzStr1[(start+i)%nzLen1])
				same++;
		}
		maxSame = MAX(maxSame, same);
	}
	return maxSame/(double)MAX(nzLen1, nzLen2);
}

void AnalyzeFile(char* file){
	FILE* pFile =fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(190);
	}
	
	for(int i=0; i<5; i++) fscanf(pFile, "%*s");
	int length;
	fscanf(pFile, "%i", &length);
	
	char* str[2];
	for(int i=0; i<2; i++)
		str[i] = malloc(sizeof(char)*(length+1));
	
	int slot=0;
	int initial=0;
	while(fscanf(pFile, "%*i %*i %*i %s", str[slot])==1){
		if(initial){
			printf("%i %lf\n", initial, CompareString(str[slot], str[slot^0x1], length));
		}
		initial++;
		slot^=0x1;
	}
	
	
}
