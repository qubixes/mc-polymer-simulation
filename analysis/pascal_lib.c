#include "pascal_lib.h"

PasData* ReadInfoFile(char* infoFile){
	PasData* pData = malloc(sizeof(PasData));
	
	FILE* pFile = fopen(infoFile, "r");
	
	if(!pFile){
		printf("Error reading file %s\n", infoFile);
		exit(192);
	}
	
	int iMono=0, nPol=0, lastPol=1234567890, curPol;
	
	while( fscanf(pFile, "%i %*s %i %*s %*s %*s %*s %*s %*s", &iMono, &curPol) == 2){
		if(abs(curPol) != lastPol){
			nPol++;
			lastPol = abs(curPol);
		}
	}
	nPol *= 2;
	rewind(pFile);
	pData->nPol = nPol;
	pData->nMono = malloc(sizeof(int)*nPol);
	pData->centromers = malloc(sizeof(int)*nPol);
	pData->polId = malloc(sizeof(int)*nPol);
	pData->pos = malloc(sizeof(double**)*nPol);
	
	int iPol=0; lastPol=396782797;
	int nCurBonds=0; 
	
	for(int iPol=0; iPol<nPol; iPol++)
		pData->centromers[iPol] = 0;
	
	while( fscanf(pFile, "%*i %*s %i %*s %*s %*s %*s %*s %*s", &curPol) == 1 ){
		if(lastPol == 396782797){
			nCurBonds=1;
			lastPol=curPol;
		}
		else if( curPol && curPol == -lastPol){
			pData->centromers[iPol] = nCurBonds;
			lastPol = curPol;
			nCurBonds++;
		}
		else if( abs(curPol) != abs(lastPol)){
			pData->nMono[iPol] = nCurBonds+1;
			pData->polId[iPol] = abs(lastPol);
			nCurBonds=1;
			iPol++;
			lastPol = curPol;
		}
		else
			nCurBonds++;
	}
	pData->nMono[iPol] = nCurBonds+1;
	pData->polId[iPol] = lastPol;
	
	for(int iPol=0; iPol<nPol/2; iPol++){
		pData->nMono[iPol+nPol/2] = pData->nMono[iPol];
		pData->polId[iPol+nPol/2] = pData->polId[iPol];
		pData->centromers[iPol+nPol/2] = pData->centromers[iPol];
	}
	
	fclose(pFile);
	
	pData->nTotMono = 0;
	pData->maxNMono = 0;
	for(int iPol=0; iPol<nPol; iPol++){
		pData->nTotMono += pData->nMono[iPol];
		pData->maxNMono = MAX(pData->maxNMono, pData->nMono[iPol]);
		pData->pos[iPol] = malloc(sizeof(double*)*pData->nMono[iPol]);
		for(int iMono=0; iMono<pData->nMono[iPol]; iMono++)
			pData->pos[iPol][iMono] = malloc(sizeof(double)*3);
	}
	
	return pData;
}



PasData* PascalDataInit(char* infoFile, char* configFile){
	PasData* pData = ReadInfoFile(infoFile);
	
	FILE* pFile = fopen(configFile, "r");
	if(!pFile){
		printf("Errof opening file %s\n", configFile);
		exit(192);
	}
	
	for(int iPol=0; iPol<pData->nPol; iPol++){
		for(int iMono=0; iMono<pData->nMono[iPol]; iMono++){
			fscanf(pFile, "%lf %lf %lf %*lf %*lf %*lf", pData->pos[iPol][iMono], pData->pos[iPol][iMono]+1, pData->pos[iPol][iMono]+2);
		}
	}
	fclose(pFile);
	return pData;
}

void PrintPasData(PasData* pData){
	printf("nPol = %i\n", pData->nPol);
	for(int iPol=0; iPol<pData->nPol; iPol++){
		printf("len = %i\n", pData->nMono[iPol]);
		printf("polId = %i\n", pData->polId[iPol]);
		printf("centromer = %i\n", pData->centromers[iPol]);
		for(int iMono=0; iMono<5 && iMono < pData->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++)
				printf("%10le ", pData->pos[iPol][iMono][k]);
			printf("\n");
		}
		printf("\n");
		for(int iMono=pData->nMono[iPol]-5; iMono<pData->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++)
				printf("%10le ", pData->pos[iPol][iMono][k]);
			printf("\n");
		}
		printf("\n");
	}
}