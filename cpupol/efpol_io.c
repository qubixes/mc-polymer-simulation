#include "efpol.h"

int ReadLatticeFile(CurState* cs, char* file){
	FILE* pFile = fopen(file, "r");
	char* strIn = malloc(sizeof(char)*(cs->polSize+1));
	int len, iPol;
	int boxT, boxU, boxV, nPol, maxLength;
	uint coor;
	uint t, u, v;
	int LT=cs->con.L, LU=cs->con.L, LV=cs->con.L;
	
	fscanf(pFile, "%*s %i", &boxT);
	fscanf(pFile, "%*s %i", &boxU);
	fscanf(pFile, "%*s %i", &boxV);
	fscanf(pFile, "%*s %i", &nPol);
	fscanf(pFile, "%*s %i", &maxLength);
	
	if(boxT != LT || boxU != LU || boxV != LV || nPol != cs->nPol || maxLength != cs->polSize){
		printf("Error loading polymer file. Parameters not the same! Loaded [%i %i %i %i %i]\n", boxT, boxU, boxV, nPol, maxLength);
		return 10;
	}
	printf("maxLength = %i\n", maxLength);
	
	iPol=0;
	while(fscanf(pFile, "%*s %i", &len)>0){
		if(len != cs->polSize){
			printf("Error length = %i vs %i\n", len, cs->polSize);
			exit(0);
		}
		fscanf(pFile, "%u %u %u", &t, &u, &v);
		fscanf(pFile, "%s", strIn);
		if(strlen(strIn) != cs->polSize){
			printf("Meh: %li vs %i\n", strlen(strIn), cs->polSize);
			exit(0);
		}
		coor = TUVtoCoor(t,u,v, &cs->con);
		for(int iMono=0; iMono<cs->polSize; iMono++){
			cs->coorPol[iPol+iMono*cs->nPol] = coor;
			SetLattice(coor, cs);
			uint unit = CharToHex(strIn[iMono]);
			cs->unitPol[iPol/8+iMono*cs->intsPerMono] |= unit<<(4*(iPol%8));
			coor = AddUnitToCoor(unit, coor, &cs->con);
		}
		iPol++;
	}
	return 0;
}


int WriteLatticeFile(CurState* cs, char* file){
	int iPol;
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		assert(pFile);
	}
	fprintf(pFile, "LT= %i\n", cs->con.L);
	fprintf(pFile, "LU= %i\n", cs->con.L);
	fprintf(pFile, "LV= %i\n", cs->con.L);
	fprintf(pFile, "np= %i\n", cs->nPol);
	fprintf(pFile, "maxPolLength= %i\n", cs->polSize);
	
	for(iPol=0; iPol<cs->nPol; iPol++){
		WritePolymer(cs, iPol, pFile);
	}
	fclose(pFile);
	return 0;
}

void WritePolymer(CurState* cs, int iPol, FILE* pFile){
	fprintf(pFile, "len= %i\n", cs->polSize);
	fprintf(pFile, "%u  %u  %u\t", TCoor(cs->coorPol[iPol], &cs->con), UCoor(cs->coorPol[iPol], &cs->con), VCoor(cs->coorPol[iPol], &cs->con));
	
	for(int iMono=0; iMono<cs->polSize; iMono++)
		fprintf(pFile, "%x", (cs->unitPol[iPol/8+iMono*cs->intsPerMono] >> (4*(iPol%8)))&0xf);
	fprintf(pFile, "\n");
}

void PrintUnits(){
	for(int i=1; i<0xf; i++){
		if(!IsValid(i)) continue;
		printf("0x%x -> %i %i %i\n", i, XUnit(i), YUnit(i), ZUnit(i));
	}
}

void PrintTrans(){
	for(int src1=0; src1<=0xf; src1++){
		if(src1 == 0x3 || src1 == 0xc) continue;
		for(int src2=0; src2<=0xf; src2++){
			if(src2 == 0x3 || src2 == 0xc) continue;
			if( (src1^src2) == 0xf && (src1 != 0xf && src2 != 0xf)) continue;
			printf("[%x,%x] | ", src1, src2);
			int combo = src1*16+src2;
			if(tab.curNTrans[combo]==0){
				printf (" X\n");
				continue;
			}
			
			for(int i=0; i<tab.curNTrans[combo]; i++){
				uint dst1 = tab.transTable[combo*NOPT+i]&0xf;
				uint dst2 = (tab.transTable[combo*NOPT+i]>>4)&0xf;
				if(dst1 || dst2)
					printf("[%x,%x] ", dst1, dst2);
			}
			printf("\n");
		}
	}
	PrintUnits();
}

void WriteSimulationSettings(Config* cfg, CurState* cs, int iStep){
	char file[2000];
	
	sprintf(file, "%s/simulation_settings.txt", cs->dir);
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "Start_seed = %u\n", cfg->seed);
	fprintf(pFile, "Length = %i\n", cs->polSize);
#if POL_TYPE == POL_RING
	fprintf(pFile, "Polytype = ring\n");
#else
	fprintf(pFile, "Polytype = lin\n");
#endif
	fprintf(pFile, "Density = %lf\n", cfg->density);
	fprintf(pFile, "Latsize = %i\n", cs->con.L);
	fprintf(pFile, "Start_polysize = %i\n", cfg->initPolSize);
	fprintf(pFile, "Double_step = %i\n", iStep);
	fprintf(pFile, "Interval = %li\n", cfg->interval);
	fprintf(pFile, "Equilibrated = 0\n");
	fprintf(pFile, "Npol = %i\n", cs->nPol);
	if(cfg->polModel == SL_EQUAL)
		fprintf(pFile, "Polymodel = sl_equal\n");
	else if(cfg->polModel == SL_DOUBLE)
		fprintf(pFile, "Polymodel = sl_double\n");
	else if(cfg->polModel == SL_QUAD)
		fprintf(pFile, "Polymodel = sl_quad\n");
	else 
		fprintf(pFile, "Polymodel = ??\n");
	fprintf(pFile, "Executable = efpol\n");
	fclose(pFile);
}

void PrintCoor(int coor, Constants* con){
	printf("[%i %i %i]", TCoor(coor, con), UCoor(coor, con), VCoor(coor, con));
	if(coor&(~con->BOUND_MASK)) printf("???\n");
}

void PrintPol(int iPol, CurState* cs){
	for(int i=0; i<cs->polSize; i++){
		if(i%8==0) printf("-----------------------------\n");
		printf("[%x, %x]: ", cs->coorPol[iPol+i*cs->nPol], (cs->unitPol[iPol/8+i*cs->intsPerMono]>>((iPol%8)*4))&0xf);
		PrintCoor(cs->coorPol[iPol+i*cs->nPol], &cs->con);
		if((~cs->con.BOUND_MASK)&cs->coorPol[iPol+i*cs->nPol])printf("??\n");
		printf("\n");
	}
}


int CheckIntegrity(CurState* cs){
	uint coor;
	uint* checkLattice = malloc(sizeof(uint)*cs->con.LAT_ALLOC);
	uint unit=100;
	for(int i=0; i<cs->con.LAT_ALLOC; i++) checkLattice[i]=0;
	
	for(int iPol=0; iPol<cs->nPol; iPol++){
		coor = cs->coorPol[iPol];
// 		printf("startCoor = %i\n", coor);
		for(int mono=0; mono<cs->polSize; mono++){
			if (coor != cs->coorPol[iPol+mono*cs->nPol]){
				PrintPol(iPol, cs);
				printf("Error: unit/coor integrity fault, polId=%i, monoId=%i, unit = %u, coor = %x vs %x\n", iPol, mono, unit, coor, cs->coorPol[iPol+mono*cs->nPol]);
// 				for(int i=0; i<16; i++){
// 					printf("0x%x\n", LargeToSmallUnit(SmallToLargeUnit(i, &cs->con), &cs->con));
// 				}
				return 2;
			}
			unit = (cs->unitPol[iPol/8+mono*cs->intsPerMono]>>((iPol%8)*4))&0xf;
			if(!IsValid(unit) && unit!=0x0 && !(POL_TYPE==POL_LIN && unit==0xf && mono==cs->polSize-1)){
				printf("Invalid unit vector: %x\n", unit);
				PrintPol(iPol, cs);
				return 4;
			}
			if(unit && checkLattice[coor/32]&((uint)1<<(coor%32)) && !(POL_TYPE==POL_RING && mono==cs->polSize-1)){
				PrintPol(iPol, cs);
				printf("Lattice error, double occupancy!\n");
				printf("mono=%i, polId=%i, coor=%x\n", mono, iPol, coor);
				return 3;
			}
			
			if(unit)
				checkLattice[coor/32] |= ((uint)1<<(coor%32));
// 			printf("%x+%x => %x\n", unit, coor, AddUnitToCoor(unit, coor, &cs->con));
			coor = AddUnitToCoor(unit, coor, &cs->con);
// 			exit(0);
		}
#if POL_TYPE == POL_RING
		if(coor != cs->coorPol[iPol]){
			printf("Polymer broken!\n");
			return 4;
		}
#endif
	}
	
	for(int i=0; i<cs->con.LAT_ALLOC; i++){
		if(checkLattice[i] != cs->lattice[i]){
			printf("Lattice error, inconsistency!\n");
			return 4;
		}
	}
	free(checkLattice);
	return 0; 
}

void WriteGenom(Config* cfg, CurState* cs, Result* res){
	char file[10000];
	FILE* pFile;
	sprintf(file, "%s/genom_N=%i.dat", cfg->dir, cs->polSize);
	
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file: %s\n", file);
		exit(0);
	}
	
	for(int i=0; i<cs->polSize; i++){
		if(res->genCounts[i])
			fprintf(pFile, "%i %lf\n", i, res->gen[i]/res->genCounts[i]);
	}
	fclose(pFile);
}

void WriteResults(Config* cfg, CurState* cs, Result* res){
	WriteGenom(cfg, cs, res);
}
