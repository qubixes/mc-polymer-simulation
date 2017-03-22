#include "gpupol_io.h"

void PrintPol(Polymer* pol){
	int i;
	printf("Polymer information:\nlength = %i\n", pol->length);
	printf("Start: (%i %i %i)\n", pol->startTUV.t, pol->startTUV.u, pol->startTUV.v);
	for(i=pol->labelStart; i<pol->length; i++)
		printf("%x", pol->bonds[i]);
	for(i=0; i<pol->labelStart; i++)
		printf("%x", pol->bonds[i]);
	printf("\n\n");
}

void PrintBinary(uint num){
	int i;
	for(i=0; i<32; i++){
		if(!((32-i)%6)) printf(" ");
		printf("%u", (num>>(31-i))&0x1);
	}
	printf("    ");
}

void PrintSite(char* lattice, uint t, uint u, uint v){
	
	uint site;
	
	site= GetGpuSite(t,u,v, &sp);
	printf("%x\n", lattice[site]);
}

void WriteSimulationSettings(SimProperties* sp, SimState* ss){
	char file[2000];
	
	sprintf(file, "%s/simulation_settings.txt", sp->dirOut);
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "Start_seed = %u\n", sp->seed);
	fprintf(pFile, "Length = %i\n", sp->polLength);
#if POL_TYPE == POL_RING
	fprintf(pFile, "Polytype = ring\n");
#else
	fprintf(pFile, "Polytype = lin\n");
#endif
	fprintf(pFile, "Density = %lf\n", sp->density);
	fprintf(pFile, "Latsize = %i\n", sp->LT);
	fprintf(pFile, "Start_polysize = %i\n", sp->polLength);
	fprintf(pFile, "Double_step = %i\n", sp->double_step);
	fprintf(pFile, "Interval = %li\n", sp->writeInterval);
	fprintf(pFile, "Equilibrated = %i\n", sp->equilibrated);
	fprintf(pFile, "Npol = %i\n", ss->nPol);
// 	if(cfg->polModel == SL_EQUAL)
// 		fprintf(pFile, "Polymodel = sl_equal\n");
// 	else if(cfg->polModel == SL_DOUBLE)
	fprintf(pFile, "Polymodel = sl_equal\n");
// 	else if(cfg->polModel == SL_QUAD)
// 		fprintf(pFile, "Polymodel = sl_quad\n");
// 	else 
// 		fprintf(pFile, "Polymodel = ??\n");
#if LCELL == 3
	fprintf(pFile, "Executable = gpupol2\n");
#elif LCELL == 4
	fprintf(pFile, "Executable = gpupol3\n");
#else
	THIS_ERROR_MEANS_LCELL_IS_NOT_DEFINED_AND_IT_SHOULD_BE
#endif
	fprintf(pFile, "Fast_equilibration = %i\n", ((sp->fastEq)?1:0));
	fprintf(pFile, "FEQ_timestep = %li\n", sp->fastEq);
#ifdef RELEASE
#define TOSTR2(x) #x
#define TOSTR(x) TOSTR2(x)
// #define STRINGIZE2(s) #s
// #define STRINGIZE(s) STRINGIZE2(s)
#define RELEASE_STR TOSTR(RELEASE)
// 	char rel[]= RELEASE_STR;
	fprintf(pFile, "Release = %s\n", RELEASE_STR);
#else
	fprintf(pFile, "Release = unknown\n");
#endif
	fclose(pFile);
}

// void PrintSummary(SimProperties* sp, SimState* ss){
// 	int iPol;//, polIndex;
// 	int* revTranslate;
// // 	Polymer* pol;
// 	
// 	SimState* curState;
// 	for(int i=0; i<sp->nDevices; i++){
// 		curState = ss+i;
// 		printf("device %i: %i polymers\n", i, curState->nPol);
// 		revTranslate = (int*) malloc(sizeof(int)*curState->nPol);
// 		for(iPol=0; iPol<curState->nPol; iPol++){
// // 			origPol = pol;
// 			revTranslate[curState->polNumTranslate[iPol]] = iPol;
// // 			printf("%i <=> %i\n", dev->polNumTranslate[iPol], iPol);
// 		}
// 		for(iPol=0; iPol<curState->nPol; iPol++){
// 			polIndex = revTranslate[iPol];
// // 			pol = curState->pol+polIndex;
// // 			polIndex=pol;
// // 			printf("pol %i: CMS=(%.1lf, %.1lf, %.1lf), length=%i\n", iPol, pol->oldModes[0].x, pol->oldModes[0].y, pol->oldModes[0].z, pol->length);
// 		}
// 		free(revTranslate);
// 	}
// 	printf("------------------------------------\n\n");
// }

int ReadLatticeFile(SimState* ss, SimProperties* sp, char* file){
	printf("Reading lattice file\n");
	FILE* pFile;
	uint i;
	uint iPol;
	uint length;
// 	int latSide;
	char* spol;
	int* slArray;
	uint* bonds, *tBonds;
	uint t,u,v;
	int np;
	char buf[2]={'\0','\0'};
	
	pFile = fopen(file, "r");
	if(!pFile){
		printf("Error reading file %s\n", file);
		return 199;
	}
	
	int LT, LU, LV;
	
	assert(fscanf(pFile,"%*s %i", &LT)==1);
	assert(fscanf(pFile,"%*s %i", &LU)==1);
	assert(fscanf(pFile,"%*s %i", &LV)==1);
	
	if(LT%(WST*LCELL) || LU%(WSU*LCELL) || LV%(WSV*LCELL)){
		printf("Error: Box dimensions not correct. Work group dimensions: (%i, %i, %i) vs file (%i, %i, %i)\n", WST*LCELL, WSU*LCELL, WSV*LCELL, LT, LU, LV);
		exit(202);
	}
	
	SetBoxDimension(sp, LT, LU, LV);
	SimStateInit(ss, sp);
	
	fscanf(pFile, "%*s %i", &np);
	sp->maxPolLength=0;
	fscanf(pFile, "%*s %i", &sp->maxPolLength);
// 	printf("nPol = %i\nmaxLength=%i\n", np, sp->maxPolLength);
	bonds = (uint*) malloc(sizeof(uint)*sp->maxPolLength);
	ss->pol = (Polymer*) malloc(sizeof(Polymer)*np);
	spol = (char*) malloc(sizeof(char)*sp->maxPolLength);
	slArray = (int*) malloc(sizeof(int)*sp->maxPolLength);
	tBonds = (uint*) malloc(sizeof(uint)*sp->maxPolLength);
	for(iPol=0; iPol<np; iPol++){
		fscanf(pFile, "%*s %i\n", &length);
		fscanf(pFile, "%u %u %u", &t, &u, &v);
// 		printf("len=%i\n%u %u %u\n", length,t,u,v);
		fscanf(pFile, "%s", spol);
		for(i=0; i<length; i++){
			buf[0] = spol[i];
			sscanf(buf, "%x", bonds+i);
		}
		AddPolymerToSim(ss, sp,t,u,v, length, bonds);
	}
	fclose(pFile);
	free(spol); free(slArray);
	free(bonds); free(tBonds); 
	printf("Succesfully read file %s\n", file);
	printf("Finished Reading lattice file\n");
	return 0;
}

int WriteLatticeFile(SimProperties* sp, SimState* ss, char* file){
	int iPol;
	int* revTranslate = (int*) malloc(sizeof(int)*ss->nPol);
	
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		assert(pFile);
	}
	fprintf(pFile, "LT= %i\n", sp->LT);
	fprintf(pFile, "LU= %i\n", sp->LU);
	fprintf(pFile, "LV= %i\n", sp->LV);
	fprintf(pFile, "np= %i\n", ss->nPol);
	fprintf(pFile, "maxPolLength= %i\n", sp->maxPolLength);
	for(iPol=0; iPol<ss->nPol; iPol++){
		revTranslate[ss->polNumTranslate[iPol]] = iPol;
	}
	for(iPol=0; iPol<ss->nPol; iPol++){
		WriteRingPolymer(pFile, ss->pol+revTranslate[iPol]);
	}
	free(revTranslate); fclose(pFile);
	return 0;
}

void WriteRingPolymer(FILE* pFile, Polymer* pol){
	int i;
// 	uint unit;
	fprintf(pFile, "len= %i\n", pol->length);
	fprintf(pFile, "%u  %u  %u\t", pol->startTUV.t, pol->startTUV.u, pol->startTUV.v);
	
	for(i=pol->labelStart; i<pol->length; i++)
		fprintf(pFile, "%x", pol->bonds[i]);
	for(i=0; i<pol->labelStart; i++)
		fprintf(pFile, "%x", pol->bonds[i]);
		
	fprintf(pFile, "\n");
}
