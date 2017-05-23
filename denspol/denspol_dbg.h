#include <stdlib.h>
#include <stdio.h>
#include "denspol.h"

int CheckIntegrity(CurState* cs, char* msg);
void PrintSystem(CurState* cs);
void PrintMutators(LookupTables* lt);
void PrintMoveCounts(LookupTables* lt);
void StatTopo(LookupTables* lt);
void CheckTopo(LookupTables* lt);
void PrintUnitDot(LookupTables* lt);
void ComputeSymmetry(int table[12], int depth, LookupTables* lt);
void PrintUnit(int unit);
