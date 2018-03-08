#include <stdio.h>
#include "denspol.h"
#define MIN(X,Y) ((X<Y)?X:Y)
#define TOSTR2(x) #x
#define TOSTR(x) TOSTR2(x)


int TUV2Coor(int t, int u, int v, int L);
int ValidateAddUnitVectors(int a, int b, int* c);
int AddUnitVectors(int a, int b);
int IsValid(int a);
int AddUnitToCoorPeriod(int unit, int coor, int L, int* OOB);
int AddUnitToCoorWBounds(int unit, int coor, int L, int* OOB);
int TCoor(int coor, int L);
int UCoor(int coor, int L);
int VCoor(int coor, int L);
int UnitInProd(int unit1, int unit2);
int CharToHex(char c);
void PrintCoor(int coor, int L);
int MonoPol2Id(int iMono, Polymer* pol);
void TUV2XYZ(int tuv[3], int xyz[3]);
void DTUV2XYZ(double tuv[3], double xyz[3]);
double Distance(double xyz1[3], double xyz2[3]);
void Coor2TUV(int coor, int tuv[3], int L);
int TestMoveHP(CurState* cs, LookupTables* lt, int iMono, Polymer* pol, int newUnit);
int TopoMove(CurState* cs, LookupTables* lt);
