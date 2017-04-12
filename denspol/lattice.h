#include <stdlib.h>
#include <stdio.h>

typedef struct Lattice{
	int* lat;
	int* mLat;
	int* labels;
	int* mLabels;
	int L;
	int LS;
}Lattice;

typedef struct MoveTable{
	int forwMoves[16][4][2];
	int backMoves[16][16];
	int validUnits[12];
}MoveTable;

Lattice* NewLattice(int L);
void PrintLattice(Lattice* lattice);
void SetMidLattice(int pos, int val, Lattice* lattice);
void SetLattice(int pos, int val, Lattice* lattice);
int GetMidLattice(int pos, Lattice* lattice);
int GetLattice(int pos, Lattice* lattice);
int AddUnitToCoor(int unit, int pos, Lattice* lattice);

int TUV2Coor(int t, int u, int v, Lattice* lattice);
int IsValid(int a);
int ValidateAddUnitVectors(int a, int b, int* c);
int AddUnitVectors(int a, int b);

MoveTable* NewMoveTable();
void PrintMoveTable(MoveTable* mt);
