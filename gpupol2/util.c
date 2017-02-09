#include "util.h"


// void PrintPolymer(int polSize, uint* unitPol, uint* coorPol){
//   int i;
//   
//   for(i=0; i<polSize; i++){
// 	 printf("l[%i]\t", i);
// 	 PrintUnitVector(unitPol[i]);
// 	 printf("\t");
// 	 PrintCoor(coorPol[i]);
// 	 printf("\n");
//   }
//   printf("l[%i]\t   -\t", polSize);
//   PrintCoor(coorPol[polSize]);
//   printf("\n");
// }
/*
void PrintCoor(uint coor){
  printf("(%i, %i, %i)", TC(coor), UC(coor), VC(coor));
}


void PrintUnitVector(uint vec){
	int i;
	if(vec==0) printf("   0");
	else{
		for(i=0; i<4-BitCount(vec); i++) printf(" ");
			if(TU(vec)) printf("t");
			if(UU(vec)) printf("u");
			if(VU(vec)) printf("v");
			if(WU(vec)) printf("w");
// 		for(i=3; i>=0; i--){
// 			if(vec&((uint)1<<(i*8))) printf("%c", 116+3-i);
// 		}
	}
// 	printf("(0x%x)", vec);
}*/

void PrintBonds(uint* bonds, int nBonds){
	for(int i=0; i<nBonds; i++){
		printf("%x", bonds[i]);
	}
	printf("\n");
}

int BitCount(uint vec){
	int i;
	int bc=0;
	for(i=0; i<32; i++)
		bc += (vec>>i)&(uint)1;
	return bc;
}


double TimeElapsed(struct timeval tStart, struct timeval tEnd){
	double dif;
	
	dif = tEnd.tv_sec-tStart.tv_sec + (tEnd.tv_usec-tStart.tv_usec)/1000000.0;
	return dif;
}

// void PrintCoor(Coor c){
// 	printf("(%lf %lf %lf) ", c.x,c.y,c.z);
// }

long GetTFromName(char* file){
	size_t len = strlen(file);
	char* tString = (char*) malloc(len+1);
	strcpy(tString, file);
	char *p=tString;
	long t=0;
	
	p = tString+len;
	while(*p != '/' && p>tString){
		p--;
		if(*p == '_')
			*p = '\0';
	}
	while(*p != '\0' && !(*p == 't' && *(p+1)=='=') ) p++;
	sscanf(p+2, "%li", &t);
	free(tString);
	return t;
}

unsigned char* BinaryFromFile(char* filename, size_t* n){
	FILE* pFile;
	unsigned char* binary;
	(*n)=0;
	
	pFile = fopen(filename, "rb");
	if(!pFile){
		*n = 0;
		return NULL;
	}
	while( fgetc(pFile) != EOF ) (*n)++;
	rewind(pFile);
	
	binary = (unsigned char*) malloc((*n)*sizeof(unsigned char));
	fread(binary, sizeof(unsigned char), *n, pFile);
	fclose(pFile);
	return binary;
}

float ElapsedTime(struct timeval sTime, struct timeval eTime){
	return (eTime.tv_sec-sTime.tv_sec + (eTime.tv_usec-sTime.tv_usec)/1000000.0);
}

char* CStringFromFile(char* filename){
	FILE* pFile;
	int n=0;
	char* cstring;
	char* include;
	char* incStart;
	char* lastEnd;
	char* includeEnd;
	char* tString;
	char* newString;
	char incStr[1000];
	int curSize=0;
	
	pFile = fopen(filename, "r");
	if(!pFile){
		printf("Error reading file %s\n", filename);
		return NULL;
	}
	while( fgetc(pFile) != EOF ) n++;
	rewind(pFile);
	
	newString=NULL;
	cstring = (char*) malloc((n+1)*sizeof(char));
	fread(cstring, sizeof(char), n, pFile);
	lastEnd = cstring;
	while((include = strstr(lastEnd, "#include ")) != NULL){
		incStart = include;
		while (*include != '"') include++;
		include++;
		includeEnd = include;
		while (*includeEnd != '"') includeEnd++;
		memcpy(incStr, include, includeEnd-include);
		incStr[includeEnd-include+1]='\0';
// 		puts(incStr);
		tString = CStringFromFile(incStr);
		if(!tString){
			printf("Error reading file in #include \"%s\" directive\n", incStr);
			return NULL;
		}
		newString = (char*) realloc(newString, curSize+incStart-lastEnd+strlen(tString)+1);
		memcpy(newString+curSize, lastEnd, incStart-lastEnd);
		memcpy(newString+(int)(curSize+incStart-lastEnd), tString, strlen(tString)+1);
		if(!tString) printf("???\n");
// 		puts(tString);
		curSize += incStart-lastEnd + strlen(tString);
		lastEnd = includeEnd+1;
		free(tString);
	}
	
	newString = (char*) realloc(newString, curSize+n-(lastEnd-cstring)+1);
	if(newString==NULL) return NULL;
	memcpy(newString + curSize, lastEnd, n-(lastEnd-cstring));
	newString[curSize+n-(lastEnd-cstring)]='\0';
// 	if(newString){
// 		printf("++++++++++++++++++++++++++++++++++++++++++++=\n");
// 		puts(newString);
// 		printf("---------------------------------------------\n");
// 	}
	fclose(pFile);
	free(cstring);
	return newString;
}

int CountBits(uint* ar, int nElem){
	int iBit, iWord, sl=0;
	
	for(iWord=0; iWord<nElem; iWord++){
		for(iBit=0; iBit<32; iBit++){
			sl += ((ar[iWord]>>iBit)&0x1);
		}
	}
	return sl;
}

