#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "denspol.h"

void CSInit(CurState* cs, unsigned int seed, int nPol, int polSize, int L, char* dir);
void GeneratePolymers(CurState* cs, LookupTables* lt);
void PrintMutators(LookupTables* lt);
void GenerateMutators(LookupTables* lt, char* file);
void UnitDotInit(LookupTables* lt, double bendEnergy);
void SuperTopoInit(LookupTables* lt);
void SimulationInit(CurState* cs, LookupTables* lt, unsigned int seed, double density, int polSize, int L, char* dir);
void PrintStraightTopo(LookupTables* lt);
double*** GenerateDistanceMatrix(int L);
void LoadHPFile(char* file, HPTable* hp, CurState* cs);
void ShufflePolymerIds(CurState* cs);