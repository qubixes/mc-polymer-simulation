#define IS_MAIN
#include "denspol.h"

int main(int argc, char** argv){
	LookupTables lt;
	CurState cs;
	Timer timer;
// 	char outFile[10000];
	
	cs.ss.L            = 10;
	cs.ss.polSize      = 100;
	cs.ss.density      = 7.0;
	cs.ss.seed         = 12846102;
	cs.ss.tMax         = 20000;
	cs.ss.interval     = 1000;
	cs.ss.bendEnergy   = 0.3;
	cs.ss.hpFile       = NULL;
	cs.ss.polIdShuffle = 0;
	cs.ss.hpStrength   = 0;
	
	if(argc>=11){
		cs.ss.seed       = (unsigned int)(atol(argv[1]));
		cs.ss.dir        = argv[2];
		cs.ss.density    = atof(argv[3]);
		cs.ss.tMax       = atof(argv[4]);
		cs.ss.interval   = atof(argv[5]);
		cs.ss.polSize    = atoi(argv[6]);
		cs.ss.L          = atoi(argv[7]);
		cs.ss.eeFile     = argv[8];
		cs.ss.bendEnergy = atof(argv[9]);
		cs.ss.dblStep    = atoi(argv[10]);
	}
	else{
		cs.ss.dir = malloc(sizeof(char)*100);
		cs.ss.eeFile = malloc(sizeof(char)*100);
		sprintf(cs.ss.dir   , "./data/test");
		sprintf(cs.ss.eeFile, "./ee_topo.dat");
		printf("Using test values: (supply %i+ instead of %i arguments)\n", 10, argc);
	}
	if(argc>=12)
		cs.ss.polIdShuffle= atoi(argv[11]);
	if(argc>=13) 
		cs.ss.hpFile = argv[12];
	if(argc>=14)
		cs.ss.hpStrength = atof(argv[13]);
	
	SimulationInit(&cs, &lt, cs.ss.seed, cs.ss.density, cs.ss.polSize, cs.ss.L, cs.ss.dir);
// 	PrintStraightTopo(&lt);
	long curStep=0, totStep=cs.ss.tMax*cs.polSize*cs.nPol;
	long endT = cs.curT + cs.ss.tMax;
	for(; cs.curT<endT; cs.curT += cs.ss.interval){
		WriteCS(&cs, cs.curT);
		TimerStart(&timer);
		DoMCStep(cs.ss.interval, &cs, &lt);
		curStep += cs.ss.interval*cs.polSize*cs.nPol;
// 		PrintSystem(&cs);
		printf("\b \b\b \b\rCurrently done %.1lf%% at %.2lf M/s", 100*curStep/(double)totStep, 1e-6*cs.ss.interval*cs.polSize*cs.nPol/TimerElapsed(&timer));
		fflush(NULL);
	}
	printf("\n");
	WriteCS(&cs, cs.curT);
	WriteSimulationSettings(&cs);
	CheckIntegrity(&cs, "After simulation");
// 	PrintMoveCounts(&lt);
	return 0;
}
