#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "denspol.h"


int WriteLatticeFile(CurState* cs, char* file);
void WritePolymer(CurState* cs, int iPol, FILE* pFile);
void WriteSimulationSettings(CurState* cs);
void TopoMapFromFile(LookupTables* lt, char* file);
int CSFromFile(CurState* cs, char* dir, long lastT);
void WriteCS(CurState* cs, long t);
void WriteMetaData(CurState* cs, char* file);
void WriteTopComp(LookupTables* lt, char* file);
void ReadTopComp(LookupTables* lt, char* file);
