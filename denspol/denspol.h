#ifndef __DENSPOL_H_INCLUDED__
#define __DENSPOL_H_INCLUDED__
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define POL_TYPE_RING 0
#define POL_TYPE_LIN 1
#ifndef POL_TYPE
	#define POL_TYPE POL_TYPE_RING
#endif

#define TRUE 1
#define FALSE 0
// #define TOPO_DENSE TRUE

#define MAX(X,Y) ((X>Y)?X:Y)
#define MAX_TOPO_STATES 1950238
#define NMUTATOR 96
#define BEND_LVL 10

#define NON_EXISTING -1
#define SAME_TOPO -2

/** 
 Notice: the polymer length is from now on defined as the number of bonds. This is different 
 from the GPUpol package, and also previous work. Especially in the context of denspol
 this makes much more sense, since the excluded volume is tested on the bonds, instead of
 the monomers.  
 Notice2: I don't think the above is true!
 Notice3: Yes it is?
 **/

typedef struct SimulationSettings{
	double density;
	double bendEnergy;
	char* dir;
	char* eeFile;
	char* hpFile;
	unsigned int seed;
	long tMax;
	long interval;
	int polIdShuffle;
	int L;
	int polSize;
	int dblStep;
	double hpStrength;
}SimulationSettings;

typedef struct CurState{
	int* topoState;
	int* bondOcc;
	
	int** unitPol;
	int** coorPol;
	
	int nPol;
	int polSize;
	
	int L;
	int LSize;
	long curT;
	unsigned int rngState[4];
	SimulationSettings ss;
}CurState;

typedef struct Topo{
	int bonds[6][2];
	int bondFree[6];
	int nBonds;
	int set;
	
	struct Topo* mutTopo[NMUTATOR];
	struct Topo* next;
	struct SuperTopo* supTopo;
}Topo;

typedef struct SuperTopo{
	struct Topo* topo;
	struct Topo* mutTopo[NMUTATOR];
	struct SuperTopo* topoSameBonds;
	int permBond;
	int id;
	int keyInserted;
}SuperTopo;

typedef struct KeyNode{
	int val;
	struct KeyNode* next;
	struct KeyNode* child;
	SuperTopo* sTopo;
}KeyNode;

typedef struct TopoCompact{
	int mutators[NMUTATOR];
	int sameTopo;
	int permBond;
}TopoCompact;

typedef struct HPTable{
	double** HPStrength;
	double*** distance;
	int** monoId;
	int* nInterHP;
}HPTable;

typedef struct LookupTables{
	int** mutTopo;
	TopoCompact* topComp;
	KeyNode* sameTopoTree;
	int nTopoComp;
 	int nTopo;
	int unitDot[16][16];
	double bendProb[2*BEND_LVL];
	int mutator[16][16][4][3];
	int newUnits[16][16][4][2];
	int counts[16][16][4][4];
	int validUnits[12];
	
	int mutIdTableDouble[16][16];
	int revMutTable[48][2];
	int revMutTableTriple[96][3];
	int mutIdTableTriple[16][16][16];
	HPTable hp;
}LookupTables;

int nSupTopo;

#include "file.h"
#include "denspol_init.h"
#include "denspol_lib.h"
#include "denspol_mc.h"
#include "denspol_io.h"
#include "denspol_dbg.h"
#include "timer.h"


#endif
