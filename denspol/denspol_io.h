#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "denspol.h"


int WriteLatticeFile(CurState* cs, char* file);
void WritePolymer(CurState* cs, int iPol, FILE* pFile);
void WriteSimulationSettings(CurState* cs);
void TopoMapFromFile(LookupTables* lt, char* file);
