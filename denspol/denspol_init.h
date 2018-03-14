#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "denspol.h"

void ReadArgumentsFromFile(SimulationSettings* ss, char* file);
void GeneratePolymers(CurState* cs);
void PrintMutators(LookupTables* lt);
void GenerateMutators(LookupTables* lt);
void UnitDotInit(LookupTables* lt, double bendEnergy);
void SuperTopoInit(LookupTables* lt);
void SimulationInit(CurState* cs, LookupTables* lt);
void PrintStraightTopo(LookupTables* lt);
double*** GenerateDistanceMatrix(int L);
void LoadHPFile(char* file, HPTable* hp, CurState* cs);
void ShufflePolymerIds(CurState* cs, LookupTables* lt);
void SetDefaultSS(SimulationSettings* ss);
void FillDefaultSS(SimulationSettings* ss);
void AllocPolymers(CurState* cs);
void PopulateMoveList(CurState* cs, LookupTables* lt);
void LatticeInit(CurState* cs, LookupTables* lt);
int CSFromParameters(CurState* cs, LookupTables* lt);