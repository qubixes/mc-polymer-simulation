#include "denspol_io.h"

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

void WritePolymer(CurState* cs, int iPol, FILE* pFile){
	int L = cs->L;
	fprintf(pFile, "len= %i\n", cs->polSize);
	fprintf(pFile, "%u  %u  %u\t", TCoor(cs->coorPol[iPol][0], L), UCoor(cs->coorPol[iPol][0], L), VCoor(cs->coorPol[iPol][0], L));
	
	for(int iMono=0; iMono<cs->polSize; iMono++)
		fprintf(pFile, "%x", cs->unitPol[iPol][iMono]);
	fprintf(pFile, "\n");
}


void WriteSimulationSettings(CurState* cs){
	SimulationSettings* ss = &cs->ss;
	char file[2000];
	
	sprintf(file, "%s/simulation_settings.txt", ss->dir);
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "Start_seed = %u\n", ss->seed);
	fprintf(pFile, "Length = %i\n", ss->polSize);
// #if POL_TYPE == POL_RING
	fprintf(pFile, "Polytype = ring\n");
// #else
// 	fprintf(pFile, "Polytype = lin\n");
// #endif
	fprintf(pFile, "Density = %lf\n", ss->density);
	fprintf(pFile, "Latsize = %i\n", ss->L);
	fprintf(pFile, "Start_polysize = %i\n", ss->polSize);
	fprintf(pFile, "Double_step = %i\n", 0);
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