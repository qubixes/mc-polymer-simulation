#include <stdlib.h>
#include <stdio.h>
#include "denspol.h"

void CSInit(CurState* cs, unsigned int seed, double density, int polSize, int L, char* dir);
void GeneratePolymers(CurState* cs, LookupTables* lt);
void PrintMutators(LookupTables* lt);
void GenerateMutators(LookupTables* lt, char* file);
