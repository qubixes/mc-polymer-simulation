
#include "rng.h"
#include "denspol.h"
#include "denspol_lib.h"

int StartMove(CurState* cs, LookupTables* lt, Polymer* pol);
int EndMove(CurState* cs, LookupTables* lt, Polymer* pol);
int TransMoveLinear(CurState* cs, LookupTables* lt, Polymer* pol, int iMono);

int TransMoveRing(CurState* cs, LookupTables* lt, Polymer* pol, int iMono);

int CheckMutation(int mutation, int coor, int* topoState, LookupTables* lt);
void PerformMutation(int mutation, int coor, int* topoState, LookupTables* lt);
int DiffuseMove(CurState* cs, LookupTables* lt, Polymer* pol, int iMono);
int TopoMove(CurState* cs, LookupTables* lt);
double DoMCStep(long nStep, CurState* cs, LookupTables* lt);
