#include <stdio.h>
#include "denspol.h"
#define MIN(X,Y) ((X<Y)?X:Y)

int TUV2Coor(int t, int u, int v, int L);
int ValidateAddUnitVectors(int a, int b, int* c);
int AddUnitVectors(int a, int b);
int IsValid(int a);
int AddUnitToCoor(int unit, int coor, int L);
int TCoor(int coor, int L);
int UCoor(int coor, int L);
int VCoor(int coor, int L);
int UnitInProd(int unit1, int unit2);
int CharToHex(char c);
void PrintCoor(int coor, int L);
int MonoPol2Id(int iMono, int iPol, int polSize);
void Id2MonoPol(int monoId, int* iMono, int* iPol, int polSize);
void TUV2XYZ(int tuv[3], int xyz[3]);
void Coor2TUV(int coor, int tuv[3], int L);
int TestMoveHP(CurState* cs, LookupTables* lt, int iMono, int iPol, int newUnit);
int TopoMove(CurState* cs, LookupTables* lt);
