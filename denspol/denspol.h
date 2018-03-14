#ifndef __DENSPOL_H_INCLUDED__
#define __DENSPOL_H_INCLUDED__
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define POL_TYPE_RING 0
#define POL_TYPE_LIN 1
#define POL_TYPE_MIXED 2

#define MOVE_DIFFUSE 0
#define MOVE_TRANS_RING 1
#define MOVE_TRANS_LIN 2
#define MOVE_TOPO 3
#define MOVE_START 4
#define MOVE_END 6

#define BOUNDARY_PERIODIC 0 
#define BOUNDARY_STATIC 1
#define SS_DEFAULT (-1)

#define LATTICE_SHAPE_EMPTY 0
#define LATTICE_SHAPE_SPHERE 1
#define LATTICE_SHAPE_CUSTOM 2

#define TRUE 1
#define FALSE 0

#define MAX(X,Y) ((X>Y)?X:Y)
#define MAX_TOPO_STATES 1950238
#define NMUTATOR 96
#define BEND_LVL 10
// #ifndef __FROM_TOPO_COMP__
// 	#define __FROM_TOPO_COMP__ TRUE
// #endif

#define NON_EXISTING -1
#define SAME_TOPO -2

/** 
  * Notice: the polymer length is from now on defined as the number of bonds. This is different 
  * from the GPUpol package, and also previous work. Especially in the context of denspol
  * this makes much more sense, since the excluded volume is tested on the bonds, instead of
  * the monomers.  
  **/

typedef struct SimulationSettings{
	double density;
	double bendEnergy;
	char* dir;
	char* eeFile;
	char* contactFile;
	unsigned int seed;
	long tMax;
	long interval;
// 	int polIdShuffle;
	int L;
	int nPol;
	int polSize;
	int polType;
	int dblStep;
	int boundaryCond;
	int latticeShape;
	int useTopo;
	double hpStrength;
}SimulationSettings;

typedef struct Polymer{
	int* unitPol;
	int* coorPol;
	int polSize;
	int nMono;
	int polType;
	int origNMono;
}Polymer;

typedef struct CurState{
	int* topoState;
	int* bondOcc;
	
	Polymer* pol;
	int nPol;
	int maxNMono;
	int nTotMono;
	int (*AddUnitToCoor)(int,int,int,int*); ///Sometimes C++ is a little nicer...
	
// 	int nLatticeUsed;
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

typedef struct TripleHP{
	int iPol, iMono;
	double strength;
}TripleHP;

typedef struct HPTable{
	double*** distance;
	TripleHP*** inter;
	int** nInter;
}HPTable;

typedef struct LookupTables{
	int** moveChoiceTable; // MOVE_TYPE | POLYMER_ID | MONO_ID
	int nMoveChoice;
	int* latticeUsed;
	int nLatticeUsed;
	int nBondUsed;
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
