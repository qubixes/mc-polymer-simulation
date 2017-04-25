#define IS_MAIN
#include "denspol.h"

int main(int argc, char** argv){
	LookupTables lt;
	CurState cs;
	Timer timer;
	char outFile[10000];
	
	cs.ss.L=100;
	cs.ss.polSize=100;

	cs.ss.density=cs.ss.polSize/(double)(cs.ss.L*cs.ss.L*cs.ss.L);
	cs.ss.seed=12846102;
	cs.ss.tMax=30000000;
	cs.ss.interval=100000;
	
	if(argc==9){
		cs.ss.seed=      (unsigned int)(atol(argv[1]));
		cs.ss.dir=       argv[2];
		cs.ss.density=   atof(argv[3]);
		cs.ss.tMax=      atof(argv[4]);
		cs.ss.interval=  atof(argv[5]);
		cs.ss.polSize=   atoi(argv[6]);
		cs.ss.L=         atoi(argv[7]);
		cs.ss.eeFile=    argv[8];
	}
	else{
		cs.ss.dir = malloc(sizeof(char)*100);
		cs.ss.eeFile = malloc(sizeof(char)*100);
		sprintf(cs.ss.dir, "./data/test");
		sprintf(cs.ss.eeFile, "./ee_topo.dat");
		printf("Using test values: (supply %i instead of %i arguments)\n", 8, argc);
	}
	
	CSInit(&cs, cs.ss.seed, cs.ss.density, cs.ss.polSize, cs.ss.L, cs.ss.dir);
	GenerateMutators(&lt, cs.ss.eeFile);
// 	CheckTopo(&lt);
// 	PrintMutators(&lt); //exit(192);
// 	StatTopo(&lt); //exit(0);
	GeneratePolymers(&cs, &lt);
	
// 	DoMCStep(748300, &cs, &lt);

	long curStep=0, totStep=cs.ss.tMax*cs.polSize*cs.nPol;
	for(long t=0; t<cs.ss.tMax; t+= cs.ss.interval){
		sprintf(outFile, "%s/t=%li_dev=%i.res", cs.ss.dir, t, 0);
		WriteLatticeFile(&cs, outFile);
		TimerStart(&timer);
		DoMCStep(cs.ss.interval, &cs, &lt);
		curStep += cs.ss.interval*cs.polSize*cs.nPol;
// 		PrintSystem(&cs);
		printf("\b \b\b \b\rCurrently done %.1lf%% at %.2lf M/s", 100*curStep/(double)totStep, 1e-6*cs.ss.interval*cs.polSize*cs.nPol/TimerElapsed(&timer));
		fflush(NULL);
	}
	printf("\n");
	sprintf(outFile, "%s/t=%li_dev=%i.res", cs.ss.dir, ((cs.ss.tMax-1)/cs.ss.interval+1)*cs.ss.interval, 0);
	WriteLatticeFile(&cs, outFile);
	WriteSimulationSettings(&cs);
// 	PrintMoveCounts(&lt);
	return 0;
}
