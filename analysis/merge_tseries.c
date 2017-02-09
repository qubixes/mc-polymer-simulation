#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct Data{
	long* t;
	double** val;
	int nTime;
	int nCol;
	int nAlloc;
	int nComment;
	char* comments;
}Data;

void ReadData(Data* dat, char* file);
void PrintData(Data* dat);
void MergeData(Data* sDat, Data* lDat, Data* mDat, int rouse);

int main(int argc, char** argv){
	Data datShort, datLong, datMerge;
	int rouseFile=0;
	
	if(argc < 3) {
		printf("Error: need two files to merge\n");
		return 192;
	}
	
	if(argc >3) rouseFile=1;
	
	char* shortFile = argv[1];
	char* longFile = argv[2];
	
	ReadData(&datShort, shortFile);
	ReadData(&datLong, longFile);
	MergeData(&datShort, &datLong, &datMerge, rouseFile);
	PrintData(&datMerge);
	
	return 0;
}

int GetNCols(FILE* pFile){
	rewind(pFile);
	char buf[20000];
	fgets(buf, 20000, pFile);
	long linePos = ftell(pFile);
	int nCol=0;
	rewind(pFile);
	for(nCol=0; ftell(pFile)<linePos; nCol++){
		fscanf(pFile, "%*s");
	}
	rewind(pFile);
	return nCol;
}


void ReadData(Data* dat, char* file){
	int allocStep=1000;
	char buf[20000];
	
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(182);
	}
		
	dat->nAlloc = allocStep;
	dat->nTime=0;
	dat->nComment=0;
	dat->nCol = GetNCols(pFile)-2;
// 	printf("#Number of columns: %i\n", dat->nCol);
	dat->val = malloc(sizeof(double*)*dat->nCol);
	dat->t = malloc(sizeof(long)*dat->nAlloc);
	for(int i=0; i<dat->nCol; i++){
		dat->val[i] = malloc(sizeof(double)*dat->nAlloc);
	}
	
	long t; char* str;
	while(1){
		fgets(buf, 20000, pFile);
		if(feof(pFile)) break;
		if(buf[0] == '#'){
			dat->comments = malloc(sizeof(char)*(strlen(buf)+1));
			dat->nComment=1;
			strcpy(dat->comments, buf);
// 			printf("%s", buf);
			continue;
		}
		
		str = buf;
		int offSet=0;
		if(sscanf(buf, "%li%n", &t, &offSet)<=0) break;
		str += offSet;
		
// 		if(fscanf(pFile, "%li", &t)<=0){
// 			if(feof(pFile)) break;
// 			else {
// 				if(fgets(buf, 20000, pFile)){
// 					fgets(buf, 20000, pFile);
// 					printf("%s", buf); 
// 					printf("At: %li\n", ftell(pFile));
// 				}
// 				continue;
// 			}
// 		}
		if(dat->nAlloc == dat->nTime){
			dat->nAlloc += allocStep;
			dat->t = realloc(dat->t, dat->nAlloc*sizeof(long));
			for(int i=0; i<dat->nCol; i++)
				dat->val[i] = realloc(dat->val[i], dat->nAlloc*sizeof(double));
		}
		for(int i=0; i<dat->nCol; i++){
			sscanf(str, "%lf%n", &(dat->val[i][dat->nTime]), &offSet);
			str += offSet;
// 			fscanf(pFile, "%lf", &(dat->val[i][dat->nTime]));
		}
		dat->t[dat->nTime++] = t;
	}
// 	printf("Read %i rows\n", dat->nTime);
	fclose(pFile);
}

void PrintData(Data* dat){
	if(dat->nComment)
		printf("%s", dat->comments);
	for(int i=0; i<dat->nTime; i++){
		printf("%li ", dat->t[i]);
		for(int j=0; j<dat->nCol; j++){
			printf("%le ", dat->val[j][i]);
		}
		printf("\n");
	}
}

double RouseInterpolate(double yt1, double yt2, long t1, long t2, double fNext){
	double a = - log(yt1/yt2)/(pow(t1,0.8)-pow(t2,0.8));
	double logy0 = log(yt1)+a*pow(t1, 0.8);
	double tInter = (1-fNext)*t1+fNext*t2;
	double yInter = exp(logy0-a*pow(tInter, 0.8));
	return yInter;
}

void MergeData(Data* sDat, Data* lDat, Data* mDat, int rouseFile){
	int nCol = sDat->nCol;
	if(rouseFile){
		for(int iCol=0; iCol<nCol; iCol++){
			double normL = lDat->val[iCol][0];
			double normS = sDat->val[iCol][0];
			double mult = normL/normS;
			for(int iTime=0; iTime<sDat->nTime; iTime++){
				sDat->val[iCol][iTime] *= mult;
			}
		}
	}
	
	mDat->nTime=0;
	mDat->nAlloc = sDat->nTime+lDat->nTime;
	mDat->nCol = nCol;
	mDat->t = malloc(mDat->nAlloc*sizeof(long));
	mDat->val = malloc(mDat->nCol*sizeof(double*));
	for(int i=0; i<nCol; i++)
		mDat->val[i] = malloc(mDat->nAlloc*sizeof(double));
	if(lDat->nComment){
		mDat->comments = malloc((strlen(lDat->comments)+1)*sizeof(char));
		mDat->nComment=1;
		strcpy(mDat->comments, lDat->comments);
	}
	else mDat->nComment=0;
	
	int sTime=0, lTime=0;
	long mInter = sDat->t[sDat->nTime-1]-lDat->t[0];
	int upd=1;
	int* useLong = malloc(sizeof(int)*nCol);
	for(int i=0; i<nCol; i++){ useLong[i]=1; }
	
	for(sTime=0; sTime < sDat->nTime; sTime++){
		while(sDat->t[sTime] > lDat->t[lTime+1]){ lTime++; upd=1;}
		
		mDat->t[mDat->nTime] = sDat->t[sTime];
		if(sDat->t[sTime] <= lDat->t[lTime]){
			for(int i=0; i<nCol; i++){
				mDat->val[i][mDat->nTime] = sDat->val[i][sTime];
			}
			mDat->nTime++;
			continue;
		}
		
		if(rouseFile && upd){
			for(int i=0; i<nCol; i++){
				int belowZero=0;
				for(int s2Time=sTime; s2Time < sDat->nTime && sDat->t[s2Time]<lDat->t[lTime+1]; s2Time++){
					if(sDat->val[i][s2Time]<=0 || lDat->t[lTime] == 0){
						belowZero=1;
						break;
					}
				}
				if(belowZero) useLong[i]=0;
				else useLong[i]=1;
			}
			upd=0;
		}

		
		double fNext = (sDat->t[sTime]-lDat->t[lTime])/(double)(lDat->t[lTime+1]-lDat->t[lTime]);
// 		printf("%li %lf %lf\n", sDat->t[sTime], fNext, fLong);
		for(int i=0; i<nCol; i++){
			double fLong = sqrt((sDat->t[sTime]-lDat->t[0])/(double)mInter);
			double longVal;
			if(rouseFile && lDat->val[i][lTime]>0 && lDat->val[i][lTime+1] > 0 && useLong[i]){
// 				printf("Hi\n");
				longVal = RouseInterpolate(lDat->val[i][lTime], lDat->val[i][lTime+1], lDat->t[lTime], lDat->t[lTime+1], fNext);
// 				if(i == 66)
// 					printf("%li %lf %lf %lf %lf %lf\n", sDat->t[sTime], fNext, fLong, lDat->val[i][lTime], lDat->val[i][lTime+1], longVal);
			}
			else if(rouseFile){
				longVal=0;
				fLong=0;
			}
			else
				longVal = (1-fNext)*lDat->val[i][lTime] + fNext*lDat->val[i][lTime+1];
			double shortVal = sDat->val[i][sTime];
			
			mDat->val[i][mDat->nTime] = (1-fLong)*shortVal + fLong*longVal;
		}
		mDat->nTime++;
	}
	while(lDat->t[lTime] <= sDat->t[sDat->nTime-1]) lTime++;
	
	for(;lTime<lDat->nTime; lTime++){
		mDat->t[mDat->nTime] = lDat->t[lTime];
		for(int i=0; i<nCol; i++){
			mDat->val[i][mDat->nTime] = lDat->val[i][lTime];
		}
		mDat->nTime++;
	}
}
