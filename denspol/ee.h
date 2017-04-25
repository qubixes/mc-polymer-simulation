#define _XOPEN_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "lattice.h"
#include "rng.h"
#include "timer.h"

typedef unsigned int uint;

#define TRUE 1
#define FALSE 0
#define END 99999
#define MAX_TOPO 10000000
#define NMUTATOR 96

#define MIN(X,Y) ((X<Y)?X:Y)
#define MAX(X,Y) ((X>Y)?X:Y)

typedef struct Config{
	int key[50];
	int pos[6][12];
	int bonds[6][12];
	int nBonds[6];
	int nPol;
	int nKey;
	struct Topo* topo;
}Config;

typedef struct Node{
	struct Node* next;
	struct Node* child;
	struct Topo* topo;
// 	Config* cfg;
	int val;
}Node;

typedef struct SuperTopo{
	struct Topo* topo;
	struct Topo* mutTopo[NMUTATOR];
	int id;
	int level;
}SuperTopo;

typedef struct Topo{
	Node* allConfigs;
	SuperTopo* supTopo;
	struct Topo* topoNext;
	int mutated;
}Topo;

typedef struct Mutator{
	int move[2];
	int bond;
	int initial;
}Mutator;

typedef struct ExactEnum{
	Lattice* lattice;
	MoveTable* mt;
	Topo* allTopo;
	Mutator allMutator[NMUTATOR];
	Node* allConfigs;
	int nMutator;
	int nConfig;
	int nTopo;
	int nNodes;
	int cornerMoves[60][5];
	int nCornerMoves;
	int maxDepth;
	int maxTopo;
}ExactEnum;

ExactEnum ee;

void LatticeSetHole(Lattice* lattice, MoveTable* mt);
void SetCornerMoves(ExactEnum* ee);
void EEInit(int maxDepth);
void EETopoInit(ExactEnum* ee);
void PrintCornerMoves(ExactEnum* ee);
void GenerateMutators(ExactEnum* ee);
void PrintMutators(ExactEnum* ee);
Node* NewNode(int val);
Node* InsertConfig(Node* node, int iKey, Config* cfg, int* inserted, Topo** destTopo);
int FindInternalCfg(ExactEnum* ee, Config* cfg, Topo** topo);
int MutateFindConfig(ExactEnum* ee, Topo* topo, Node* node, Config curConfig);
Topo* NewTopo();
Config* NewConfig();
void PrintConfig(Config* cfg);
void PrintTree(Node* node, int depth);
int ConfigOccupied(int pos, Config* cfg);
void InsertValue(int val, int iPos, int len, int* arr);
int CheckIntegrity(Config* cfg);
// void MutateTopo(ExactEnum* ee, Topo* topo);
void DirectMutateTopo(ExactEnum* ee);
void DeleteNodeTree(Node* node);
void DeleteConfig(Config* cfg);
// void WriteTopo(ExactEnum* ee, char* file);
void WriteSuperTopo(ExactEnum* ee, char* file);
void MergeSuperTopo(SuperTopo* sTopoA, SuperTopo* sTopoB);


