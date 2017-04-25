
#include "rng.h"
#include "denspol.h"

int CheckMutation(int mutation, int coor, int* topoState, LookupTables* lt);
void PerformMutation(int mutation, int coor, int* topoState, LookupTables* lt);
int TransStep(CurState* cs, LookupTables* lt);
int DiffuseStep(CurState* cs);
double DoMCStep(long nStep, CurState* cs, LookupTables* lt);
void PrintSystem(CurState* cs);
void CheckIntegrity(CurState* cs, char* msg);