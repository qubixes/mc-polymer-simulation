#include "denspol_io.h"

void WriteCS(CurState* cs, long t){
	char filePol[10000];
	char fileLat[10000];
	
	sprintf(filePol, "%s/t=%li_dev=%i.res", cs->ss.dir, t, 0);
	sprintf(fileLat, "%s/sav_t%li_dev=%i.res", cs->ss.dir, t, 0);
	
	WriteLatticeFile(cs, filePol);
	WriteMetaData(cs, fileLat);
}

int WriteLatticeFile(CurState* cs, char* file){
	int iPol;
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		assert(pFile);
	}
	fprintf(pFile, "LT= %i\n", cs->L);
	fprintf(pFile, "LU= %i\n", cs->L);
	fprintf(pFile, "LV= %i\n", cs->L);
	fprintf(pFile, "np= %i\n", cs->nPol);
	fprintf(pFile, "maxPolLength= %i\n", cs->polSize);
	
	for(iPol=0; iPol<cs->nPol; iPol++){
		WritePolymer(cs, iPol, pFile);
	}
	fclose(pFile);
	return 0;
}

int CSFromFile(CurState* cs, char* dir, long lastT){
	
	char file[3000];
	sprintf(file, "%s/t=%li_dev=0.res", dir, lastT);
	
	FILE* pFile  = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s for reading\n", file);
		exit(192);
	}
	
	int LX[3];
	for(int i=0; i<3; i++){
		fscanf(pFile, "%*s %i", LX+i);
		if(i>0 && LX[i] != LX[i-1]){
			printf("Only square lattices are supported: %i vs %i\n", LX[i], LX[i-1]);
			exit(192);
		}
	}
	int nPol, maxLength;
	fscanf(pFile, "%*s %i", &nPol);
	fscanf(pFile, "%*s %i", &maxLength);
	CSInit(cs, 1293744, nPol, maxLength, LX[0], dir);
	
	int t,u,v;
	char* strIn = malloc(sizeof(char)*(maxLength+1));
	int iPol=0;
	while(fscanf(pFile, "%*s %i", &cs->polSize)>0){
		fscanf(pFile, "%i %i %i", &t, &u, &v);
		fscanf(pFile, "%s", strIn);
		if(strlen(strIn) != cs->polSize){
			printf("Meh: %li vs %i\n", strlen(strIn), cs->polSize);
			exit(0);
		}
		int coor = TUV2Coor(t,u,v, cs->L);
		for(int iMono=0; iMono<cs->polSize; iMono++){
			int unit = CharToHex(strIn[iMono]);
			cs->coorPol[iPol][iMono] = coor;
			cs->unitPol[iPol][iMono] = unit;
			coor = AddUnitToCoor(unit, coor, cs->L);
		}
		cs->coorPol[iPol][cs->polSize] = coor;
		iPol++;
	}
	fclose(pFile);
	
	char latFile[3000];
	sprintf(latFile, "%s/sav_t%li_dev=0.res", dir, lastT);
	
	pFile = fopen(latFile, "r");
	if(!pFile){
		printf("Error opening file %s\n", latFile);
		exit(192);
	}
	
	for(int i=0; i<4; i++){
		fscanf(pFile, "%u", cs->rngState+i);
	}
	
	for(int i=0; i<cs->LSize; i++){
		fscanf(pFile, "%i %i", cs->topoState+i, cs->bondOcc+i);
	}
	fclose(pFile);
	return 0;
}

void WritePolymer(CurState* cs, int iPol, FILE* pFile){
	int L = cs->L;
#if POL_TYPE == POL_TYPE_LIN
	fprintf(pFile, "len= %i\n", cs->polSize+1);
#else
	fprintf(pFile, "len= %i\n", cs->polSize);
#endif
	fprintf(pFile, "%u  %u  %u\t", TCoor(cs->coorPol[iPol][0], L), UCoor(cs->coorPol[iPol][0], L), VCoor(cs->coorPol[iPol][0], L));
	
	for(int iMono=0; iMono<cs->polSize; iMono++)
		fprintf(pFile, "%x", cs->unitPol[iPol][iMono]);
#if POL_TYPE == POL_TYPE_LIN
	fprintf(pFile, "f");
#endif
	fprintf(pFile, "\n");
}

void WriteMetaData(CurState* cs, char* file){
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		assert(pFile);
	}
	
	for(int i=0; i<4; i++) fprintf(pFile, "%u ", cs->rngState[i]);
	fprintf(pFile, "\n");
	for(int coor=0; coor<cs->LSize; coor++){
		fprintf(pFile, "%i %i\n", cs->topoState[coor], cs->bondOcc[coor]);
	}
	fclose(pFile);
}



void WriteSimulationSettings(CurState* cs){
	SimulationSettings* ss = &cs->ss;
	char file[2000];
	
	sprintf(file, "%s/simulation_settings.txt", ss->dir);
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "Start_seed = %u\n", ss->seed);
#if POL_TYPE == POL_TYPE_RING
	fprintf(pFile, "Length = %i\n", ss->polSize);
#elif POL_TYPE == POL_TYPE_LIN
	fprintf(pFile, "Length = %i\n", ss->polSize+1);
#endif
#if POL_TYPE == POL_TYPE_RING
	fprintf(pFile, "Polytype = ring\n");
#else
	fprintf(pFile, "Polytype = lin\n");
#endif
	fprintf(pFile, "Density = %lf\n", ss->density);
	fprintf(pFile, "Latsize = %i\n", ss->L);
	fprintf(pFile, "Start_polysize = %i\n", ss->polSize);
	fprintf(pFile, "Double_step = %i\n", ss->dblStep);
	fprintf(pFile, "Interval = %li\n", ss->interval);
	fprintf(pFile, "Equilibrated = 0\n");
	fprintf(pFile, "Npol = %i\n", cs->nPol);
// 	if(cfg->polModel == SL_EQUAL)
// 		fprintf(pFile, "Polymodel = sl_equal\n");
// 	else if(cfg->polModel == SL_DOUBLE)
// 		fprintf(pFile, "Polymodel = sl_double\n");
// 	else if(cfg->polModel == SL_QUAD)
// 		fprintf(pFile, "Polymodel = sl_quad\n");
// 	else 
// 		fprintf(pFile, "Polymodel = ??\n");
	fprintf(pFile, "Executable = denspol\n");
	fprintf(pFile, "Bend_energy = %lf\n", ss->bendEnergy);
	fclose(pFile);
}

void TopoMapFromFile(LookupTables* lt, char* file){
	lt->mutTopo = malloc(sizeof(int*)*MAX_TOPO_STATES);
	for(int i=0; i<MAX_TOPO_STATES; i++){
		lt->mutTopo[i] = malloc(sizeof(int)*NMUTATOR);
		for(int j=0; j<NMUTATOR; j++)
			lt->mutTopo[i][j] = -1;
	}
	
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	int src, mut, dst, max=-1;
	while(fscanf(pFile, "%i %i %i", &src, &mut, &dst) == 3){
		lt->mutTopo[src][mut] = dst;
		if(src > max) max = src;
	}
	
	if(max>MAX_TOPO_STATES) {
		printf("Allocated not enough memory!\n");
		exit(0);
	}
	
	lt->nTopo = max+1;
	printf("Found %i topo states\n", lt->nTopo);
	for(int i=lt->nTopo; i<MAX_TOPO_STATES; i++){
		free(lt->mutTopo[i]);
	}
	lt->mutTopo = realloc(lt->mutTopo, sizeof(int*)*lt->nTopo);
}