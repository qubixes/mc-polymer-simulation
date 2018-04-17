#include "pascal_lib.h"

Data* NewData(int maxNMono, int nPol, int L, int boundaryCond){
	Data* data    = malloc(sizeof(Data));
	data->nPol    = nPol;
	data->maxNMono= maxNMono;
	data->L       = L;
	
	data->tuv = malloc(sizeof(int**)*nPol);
	for(int i=0; i<nPol; i++){
		data->tuv[i] = malloc(sizeof(int*)*maxNMono);
		for(int j=0; j<maxNMono; j++)
			data->tuv[i][j] = malloc(sizeof(int)*3);
	}
	
	data->nMono    = malloc(sizeof(int)*nPol);
	data->polTypes = malloc(sizeof(int)*nPol);
	
	if(boundaryCond == BOUNDARY_PERIOD)
		data->TUV2Distance = &TUV2DistancePeriodic;
	else
		data->TUV2Distance = &TUV2DistanceStatic;
	
	return data;
}

Data* ReadData(char* file, int boundaryCond){
	int nPol, maxNMono, L;
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int i=0; i<3; i++)
		fscanf(pFile, "%*s %i", &L);
	
	fscanf(pFile, "%*s %i", &nPol);
	
	fscanf(pFile, "%*s %i", &maxNMono);
	Data* data = NewData(maxNMono, nPol, L, boundaryCond);
	char* str  = malloc(sizeof(char)*(maxNMono+1));
	
	data->nTotMono=0;
	for(int iPol=0; iPol<nPol; iPol++){
		int tuv[3];
		fscanf(pFile, "%*s %i %i %i %i %s", &data->nMono[iPol], tuv, tuv+1, tuv+2, str);
		data->nTotMono += data->nMono[iPol];
		for(int iMono=0; iMono<data->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++)
				data->tuv[iPol][iMono][k] = tuv[k];
			AddCharToTUV(str[iMono], tuv);
		}
		if(str[data->nMono[iPol]] == 'f') data->polTypes[iPol] = POL_TYPE_LIN;
		else data->polTypes[iPol] = POL_TYPE_RING;
	}
	fclose(pFile);
	return data;
}

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

void AddCharToTUV(char c, int tuv[3]){
	int bond = CharToHex(c);
	
	for(int i=0; i<3; i++){
		tuv[i] += (bond>>i)&0x1;
		tuv[i] -= (bond>>3)&0x1;
	}
}

void DblTUV2XYZ(double tuv[3], double xyz[3]){
	xyz[0] = (tuv[0]+tuv[1]-tuv[2])/sqrt(2);
	xyz[1] = (tuv[0]-tuv[1]       )/sqrt(2);
	xyz[2] = (              tuv[2])/sqrt(2);
}


double Squabs(double xyz[3]){
	int tot=0;
	for(int k=0; k<3; k++) tot += xyz[k]*xyz[k];
	return tot;
}

/// Calculate the distance between two points on the lattice with periodic boundary conditions.
/// 
/// It's looks a bit weird, but dtuv starts at the (-L,-L,-L) - (0,0,0) periodic box 
/// for slightly easier indexing.

double TUV2DistancePeriodic(double tuv1[3], double tuv2[3], int L){
	double dxyz[3];
	
	double dtuv[3], newtuv[3];
	
	
	
	for(int k=0; k<3; k++){
		while(tuv1[k] <  0) tuv1[k] += L;
		while(tuv1[k] >= L) tuv1[k] -= L;
		while(tuv2[k] <  0) tuv2[k] += L;
		while(tuv2[k] >= L) tuv2[k] -= L;
		dtuv[k] = tuv1[k]-tuv2[k]-L;
	}
	
	double minDist=L*L*L*L;
	
	newtuv[0]=dtuv[0];
	for(int dt=-L; dt <= L; dt += L, newtuv[0] += L){
		newtuv[1]=dtuv[1];
		for(int du=-L; du <= L; du += L, newtuv[1] += L){
			newtuv[2]=dtuv[2];
			for(int dv=-L; dv <= L; dv += L, newtuv[2] += L){
				DblTUV2XYZ(newtuv, dxyz);
				double newDist = Squabs(dxyz);
				minDist = MIN(minDist, newDist);
			}
		}
	}
	
	return minDist;
}

double TUV2DistanceStatic(double tuv1[3], double tuv2[3], int L){
	double dtuv[3], dxyz[3];
	
	for(int k=0; k<3; k++)
		dtuv[k] = tuv1[k]-tuv2[k];
	DblTUV2XYZ(dtuv, dxyz);
	return Squabs(dxyz);
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