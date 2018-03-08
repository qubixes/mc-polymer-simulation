#define IS_MAIN
#include "denspol.h"

int main(int argc, char** argv){
	LookupTables lt;
	CurState cs;
	Timer timer;
	
	SetDefaultSS(&cs.ss);
	if(argc<2){
		cs.ss.dir = malloc(sizeof(char)*100);
		cs.ss.eeFile = malloc(sizeof(char)*100);
		sprintf(cs.ss.dir   , "./data/test");
		sprintf(cs.ss.eeFile, "./ee_topo_comp.dat");
		printf("Using test values: (supply %i+ instead of %i arguments)\n", 10, argc);
	}
	else{
		ReadArgumentsFromFile(&cs.ss, argv[1]);
	}
	
	SimulationInit(&cs, &lt);
	WriteSimulationSettings(&cs);
	
	long curStep=0, totStep=cs.ss.tMax*cs.nTotMono;
	long endT = cs.curT + cs.ss.tMax;
	for(; cs.curT<endT; cs.curT += cs.ss.interval){
		WriteCS(&cs, cs.curT);
		TimerStart(&timer);
		DoMCStep(cs.ss.interval, &cs, &lt);
		curStep += cs.ss.interval*cs.nTotMono;
// 		PrintSystem(&cs);
		printf("\b \b\b \b\rCurrently done %.1lf%% at %.2lf M/s", 100*curStep/(double)totStep, 1e-6*cs.ss.interval*cs.nTotMono/TimerElapsed(&timer));
		fflush(NULL);
	}
	printf("\n");
	WriteCS(&cs, cs.curT);
	WriteSimulationSettings(&cs);
	CheckIntegrity(&cs, "After simulation");
// 	WriteTopComp(&lt, "./ee_topo_comp.dat");
// 	PrintMoveCounts(&lt);
	return 0;
}
