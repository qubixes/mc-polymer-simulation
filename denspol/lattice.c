#include "lattice.h"

Lattice* NewLattice(int L){
	Lattice* lattice = malloc(sizeof(Lattice));

	lattice->L = L;
	lattice->LS = L*L*L;
	
	lattice->lat = malloc(sizeof(int)*L*L*L);
	lattice->mLat = lattice->lat + L*L*(L/2) + L*(L/2) + (L/2);
	lattice->labels = malloc(sizeof(int)*lattice->LS);
	lattice->mLabels = lattice->labels + L*L*(L/2) + L*(L/2) + (L/2);
	
	
	for(int i=0; i<lattice->LS; i++) lattice->lat[i]=0;
	return lattice;
}

void PrintLattice(Lattice* lattice){
	int L = lattice->L;
	for(int v=0; v<L; v++){
		for(int u=0; u<L; u++){
			for(int t=0; t<L; t++){
				int pos = TUV2Coor(t,u,v, lattice);
				printf("%i ", lattice->lat[pos]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n");
}

void SetMidLattice(int pos, int val, Lattice* lattice){
	lattice->mLat[pos] = val;
}

void SetLattice(int pos, int val, Lattice* lattice){
	lattice->lat[pos] = val;
}

int GetMidLattice(int pos, Lattice* lattice){
	return lattice->mLat[pos];
}

int GetLattice(int pos, Lattice* lattice){
	return lattice->lat[pos];
}

int AddUnitToCoor(int unit, int pos, Lattice* lattice){
	int L = lattice->L;
	int dw =  (unit>>3)&0x1;
	
	int dt =  (unit&0x1)     - dw;
	int du = ((unit>>1)&0x1) - dw;
	int dv = ((unit>>2)&0x1) - dw;
	
	return pos + dt + du*L + dv*L*L;
}

int TUV2Coor(int t, int u, int v, Lattice* lattice){
	int L = lattice->L;
	return t+L*u+L*L*v;
}

int IsValid(int a){
	return (a!=0 && a!=0x3 && a!=0xc && a!=0xf);
}

int ValidateAddUnitVectors(int a, int b, int* c){
	int r, valid;
	if((a|b) != 0xf && (a&b))
		return 0;
	r = (((a|b)==0xf)?(a&b):(a|b));
	valid = IsValid(r);
	*c = r;
	return valid;
}

int AddUnitVectors(int a, int b){
	return (((a|b)==0xf)?(a&b):(a|b));
}

MoveTable* NewMoveTable(){
	MoveTable* mt = malloc(sizeof(MoveTable));
	int nForwMoves[16];
	for(int i=0; i<16; i++){
		nForwMoves[i]=0;
		for(int j=0; j<16; j++)
			mt->backMoves[i][j]=0;
	}
	
	int unitC;
	int nValid=0;
	for(int unitA=1; unitA<0xf; unitA++){
		if(!IsValid(unitA)) continue;
		mt->validUnits[nValid++] = unitA;
		for(int unitB=unitA+1; unitB<0xf; unitB++){
			if(!IsValid(unitB)) continue;
			
			if(ValidateAddUnitVectors(unitA, unitB, &unitC)){
				mt->forwMoves[unitC][nForwMoves[unitC]  ][0] = unitA;
				mt->forwMoves[unitC][nForwMoves[unitC]++][1] = unitB;
				mt->forwMoves[unitC][nForwMoves[unitC]  ][0] = unitB;
				mt->forwMoves[unitC][nForwMoves[unitC]++][1] = unitA;
				mt->backMoves[unitA][unitB] = unitC;
				mt->backMoves[unitB][unitA] = unitC;
			}
		}
	}
	return mt;
}

void PrintMoveTable(MoveTable* mt){
	for(int unitA=1; unitA<0xf; unitA++){
		if(!IsValid(unitA)) continue;
		printf("%x => ", unitA);
		for(int i=0; i<4; i++) printf(" %x + %x ", mt->forwMoves[unitA][i][0], mt->forwMoves[unitA][i][1]);
		printf("\n");
	}
	
	for(int unitA=1; unitA<0xf; unitA++){
		for(int unitB=1; unitB<0xf; unitB++){
			if(mt->backMoves[unitA][unitB]){
				printf("%x + %x => %x\n", unitA, unitB, mt->backMoves[unitA][unitB]);
			}
		}
	}
}

