#ifndef __DENSPOL_H_INCLUDED__
#define __DENSPOL_H_INCLUDED__
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
// #define N_TOPO_STATES 1950237
// #define N_TOPO_STATES 824062
#define MAX_TOPO_STATES 1950238
#define NMUTATOR 96
#define BEND_LVL 10

typedef struct SimulationSettings{
	double density;
	double bendEnergy;
	char* dir;
	char* eeFile;
	unsigned int seed;
	long tMax;
	long interval;
	int L;
	int polSize;
}SimulationSettings;

typedef struct CurState{
	int* topoState;
	
	int** unitPol;
	int** coorPol;
	
	int nPol;
	int polSize;
	
	int L;
	int LSize;
	unsigned int rngState[4];
	SimulationSettings ss;
}CurState;

typedef struct Topo{
	int bonds[6][2];
	int bondFree[6];
	int nBonds;
	
	struct Topo* mutTopo[NMUTATOR];
	struct Topo* next;
	struct SuperTopo* supTopo;
}Topo;


typedef struct SuperTopo{
	struct Topo* topo;
}SuperTopo;

typedef struct LookupTables{
	int** mutTopo;
	SuperTopo* sTopo;
	int nTopo;
	int unitDot[16][16];
	double bendProb[2*BEND_LVL];
	int mutator[16][16][4][3];
	int newUnits[16][16][4][2];
	int counts[16][16][4][4];
	int validUnits[12];
}LookupTables;

int nSupTopo;

#include "denspol_init.h"
#include "denspol_lib.h"
#include "denspol_mc.h"
#include "denspol_io.h"
#include "denspol_dbg.h"
#include "timer.h"


#endif
