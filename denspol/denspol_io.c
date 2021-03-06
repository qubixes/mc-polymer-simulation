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
	fprintf(pFile, "maxPolLength= %i\n", cs->maxNMono);
	
	for(iPol=0; iPol<cs->nPol; iPol++){
		WritePolymer(cs, cs->pol+iPol, pFile);
	}
	fclose(pFile);
	return 0;
}

int ReadLengthFile(CurState* cs, LookupTables* lt){
	char* file = cs->ss.contactFile;
	FILE* pFile= fopen(file, "r");
	if(!pFile) printf("Error opening file %s for reading\n", file);
	
	char polType[200];
	cs->nPol=0;
	cs->maxNMono=0;
	cs->nTotMono=0;
	fscanf(pFile, "%*s %i", &cs->nPol);
	fclose(pFile);
	
	FillDefaultSS(&cs->ss);
	Seed(cs->rngState, cs->ss.seed);
	if(cs->ss.boundaryCond == BOUNDARY_PERIODIC)
		cs->AddUnitToCoor = &AddUnitToCoorPeriod;
	else if(cs->ss.boundaryCond == BOUNDARY_STATIC)
		cs->AddUnitToCoor = &AddUnitToCoorWBounds;
	LatticeInit(cs, lt);
	
	int* newNMono = ComputePolLengths(file, lt->nBondUsed, cs->ss.density);
	for(int iPol=0; iPol<cs->nPol; iPol++) cs->maxNMono = MAX(cs->maxNMono, newNMono[iPol]);
	AllocPolymers(cs);
	
	pFile = fopen(file, "r");
	for(int i=0; i<6; i++) fscanf(pFile, "%*s");
	for(int iPol=0; iPol<cs->nPol; iPol++){
		Polymer* pol = cs->pol+iPol;
		fscanf(pFile, "%s %*i", polType);
		pol->nMono = newNMono[iPol];
		cs->nTotMono += pol->nMono;
		if(!strcmp(polType, "lin")){
			pol->polType = POL_TYPE_LIN;
			pol->polSize = pol->nMono-1;
		}
		else{
			pol->polType = POL_TYPE_RING;
			pol->polSize = pol->nMono;
		}
	}
	fclose(pFile);
	return 0;
}

int CSFromFile(CurState* cs, LookupTables* lt, long lastT){
	char* dir = cs->ss.dir;
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
	cs->ss.L            = LX[0];
	cs->ss.seed         = 1294812;
	cs->ss.nPol         = nPol;
	cs->maxNMono        = maxLength;
	cs->ss.latticeShape = LATTICE_SHAPE_CUSTOM;
	cs->nPol            = cs->ss.nPol;
	LatticeInit(cs, lt);

	FillDefaultSS(&cs->ss);
	AllocPolymers(cs);
	if(cs->ss.boundaryCond == BOUNDARY_PERIODIC)
		cs->AddUnitToCoor = &AddUnitToCoorPeriod;
	else if(cs->ss.boundaryCond == BOUNDARY_STATIC)
		cs->AddUnitToCoor = &AddUnitToCoorWBounds;

	
	int t,u,v;
	char* strIn = malloc(sizeof(char)*(maxLength+1));
	Polymer* pol = cs->pol;
	while(fscanf(pFile, "%*s %i", &pol->nMono)>0){
		fscanf(pFile, "%i %i %i", &t, &u, &v);
		fscanf(pFile, "%s", strIn);
		if(strlen(strIn) != pol->nMono){
			printf("Error in reading polymer: length %li while expecting %i\n", strlen(strIn), pol->nMono);
			printf("%s\n", strIn);
			exit(0);
		}
		int coor = TUV2Coor(t,u,v, cs->L);
		for(int iMono=0; iMono<pol->nMono; iMono++){
			int unit = CharToHex(strIn[iMono]);
			pol->coorPol[iMono] = coor;
			pol->unitPol[iMono] = unit;
			int OOB;
			coor = cs->AddUnitToCoor(unit, coor, cs->L, &OOB);
			if(OOB){
				printf("Error: detecting out of bounds while reading file %s. Periodic boundary conditions?\n", file);
				exit(192);
			}
		}
		if(pol->unitPol[pol->nMono-1] == 0xf){
			pol->polType = POL_TYPE_LIN;
			pol->polSize = pol->nMono-1;
		}
		else{
			pol->polType = POL_TYPE_RING;
			pol->polSize = pol->nMono;
		}
		
		pol++;
	}
	fclose(pFile);
	
	cs->nTotMono=0;
	for(Polymer* pol=cs->pol; pol<cs->pol+cs->nPol; pol++)
		cs->nTotMono += pol->nMono;
	
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

void WritePolymer(CurState* cs, Polymer* pol, FILE* pFile){
	int L = cs->L;
	fprintf(pFile, "len= %i\n", pol->nMono);
	fprintf(pFile, "%u  %u  %u\t", TCoor(pol->coorPol[0], L), UCoor(pol->coorPol[0], L), VCoor(pol->coorPol[0], L));
	
	for(int iMono=0; iMono<pol->polSize; iMono++)
		fprintf(pFile, "%x", pol->unitPol[iMono]);
	if(pol->polType == POL_TYPE_LIN)
		fprintf(pFile, "f");
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
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	fprintf(pFile, "Start_seed = %u\n", ss->seed);
	fprintf(pFile, "Length = %i\n", ss->polSize);
	if(ss->polType == POL_TYPE_RING)
		fprintf(pFile, "Polytype = ring\n");
	else if(ss->polType == POL_TYPE_LIN)
		fprintf(pFile, "Polytype = lin\n");
	else
		fprintf(pFile, "Polytype = mixed\n");
	fprintf(pFile, "Density = %lf\n", ss->density);
	fprintf(pFile, "Latsize = %i\n", ss->L);
	fprintf(pFile, "Start_polysize = %i\n", ss->polSize);
	fprintf(pFile, "Double_step = %i\n", ss->dblStep);
	fprintf(pFile, "Interval = %li\n", ss->interval);
	fprintf(pFile, "Equilibrated = 0\n");
	fprintf(pFile, "Npol = %i\n", cs->nPol);
	fprintf(pFile, "Executable = denspol\n");
	fprintf(pFile, "Bend_energy = %lf\n", ss->bendEnergy);
	fprintf(pFile, "HP_strength = %lf\n", ss->hpStrength);
	fprintf(pFile, "Chain_crossing = %i\n", ss->useTopo);
#ifdef RELEASE
#define RELEASE_STR TOSTR(RELEASE)
	fprintf(pFile, "Release = %s\n", RELEASE_STR);
#else
	fprintf(pFile, "Release = unknown\n");
#endif
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
// 	printf("Found %i topo states\n", lt->nTopo);
	for(int i=lt->nTopo; i<MAX_TOPO_STATES; i++){
		free(lt->mutTopo[i]);
	}
	lt->mutTopo = realloc(lt->mutTopo, sizeof(int*)*lt->nTopo);
}

void WriteTopComp(LookupTables* lt, char* file){
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "nTopo= %i\n", lt->nTopoComp);
	for(int i=0; i<lt->nTopoComp; i++){
		fprintf(pFile, "%i %i %i", i, lt->topComp[i].sameTopo, lt->topComp[i].permBond);
		for(int iMut=0; iMut<NMUTATOR; iMut++)
			fprintf(pFile, " %i", lt->topComp[i].mutators[iMut]);
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

void ReadTopComp(LookupTables* lt, char* file){
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	fscanf(pFile, "%*s %i", &lt->nTopoComp);
	lt->topComp = malloc(sizeof(TopoCompact)*lt->nTopoComp);
	for(int iTopo=0; iTopo<lt->nTopoComp; iTopo++){
		fscanf(pFile, "%*i %i %i", &lt->topComp[iTopo].sameTopo, &lt->topComp[iTopo].permBond);
		for(int iMut=0; iMut<NMUTATOR; iMut++)
			fscanf(pFile, "%i", &lt->topComp[iTopo].mutators[iMut]);
	}
	fclose(pFile);
}

