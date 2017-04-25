#include <stdlib.h>
#include <stdio.h>
#include "denspol.h"


void CheckIntegrity(CurState* cs, char* msg);
void PrintSystem(CurState* cs);
void PrintMutators(LookupTables* lt);
void PrintMoveCounts(LookupTables* lt);
void StatTopo(LookupTables* lt);
void CheckTopo(LookupTables* lt);