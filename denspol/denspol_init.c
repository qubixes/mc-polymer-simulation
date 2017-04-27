#include "denspol.h"

void CSInit(CurState* cs, unsigned int seed, double density, int polSize, int L, char* dir){
	char exec[1000];
	int nPol = (int)(L*L*L*density/(double)polSize+0.5);
	
	cs->ss.dir = dir;
	cs->L = L;
	cs->LSize = L*L*L;
	cs->nPol = nPol;
	cs->polSize = polSize;
	cs->ss.density = density;
	
	cs->topoState = malloc(sizeof(int)*cs->LSize);
	cs->unitPol = malloc(sizeof(int*)*cs->nPol);
	cs->coorPol = malloc(sizeof(int*)*cs->nPol);
	
	for(int iPol=0; iPol<cs->nPol; iPol++){
		cs->unitPol[iPol] = malloc(sizeof(int)*cs->polSize);
		cs->coorPol[iPol] = malloc(sizeof(int)*cs->polSize);
	}
	
	for(int i=0; i<cs->LSize; i++) cs->topoState[i]=0;
	
	Seed(cs->rngState, seed);
	
	sprintf(exec, "mkdir -p %s", dir);
	system(exec);
}

void GeneratePolymers(CurState* cs, LookupTables* lt){
	int delta;
	int L = cs->L;
	
	for(delta=L; delta>=1; delta--){
		int nTot = 1;
		nTot *= L/(delta+1);
		nTot *= L/delta;
		nTot *= L/(delta+1);
		if(nTot >= cs->nPol){
			break;
		}
	}
	if(!delta){
		printf("Failed to find partition for polymers\n");
		exit(102);
	}
	
// 	printf("Delta=%i\n", delta);
	
// 	int nPol=0;
	
	int units[3] = {0x4,0x1,0xa};
	int topoState[3];
	
	for(int i=0; i<3; i++){
		int j;
		for(j=0; j<4; j++){
			if(lt->newUnits[(~units[i])&0xf][0][j][0] != ((units[(i+1)%3])&0xf)) continue;
			if(lt->newUnits[(~units[i])&0xf][0][j][1] != ((units[(i+2)%3])&0xf)) continue;	
			break;
		}
// 		if(j==4) printf("????\n");
// 		printf("NewUnits = %i/%i\n", lt->newUnits[(~units[i])&0xf][0][j][0], lt->newUnits[(~units[i])&0xf][0][j][1]);
		topoState[(i+2)%3] = lt->mutTopo[0][lt->mutator[(~units[i])&0xf][0][j][1]];
	}
// 	for(int i=0; i<3; i++){
// 		printf("%i/%x ", topoState[i], units[i]);
// 	}
// 	printf("\n");
// 	exit(0);
	int iPol=0;
	for(int t=0; t<L && iPol<cs->nPol; t += delta+1){
		for(int u=0; u<L && iPol<cs->nPol; u += delta){
			for(int v=0; v<L && iPol<cs->nPol; v += delta+1){
				cs->unitPol[iPol][0] = 0x4;
				cs->unitPol[iPol][1] = 0x1;
				cs->unitPol[iPol][2] = 0xa;
				
				for(int iMono=3; iMono<cs->polSize; iMono++)
					cs->unitPol[iPol][iMono] = 0;
				int coor = TUV2Coor(t,u,v,L);
				for(int iMono=0; iMono<cs->polSize; iMono++){
					cs->coorPol[iPol][iMono] = coor;
					coor = AddUnitToCoor(cs->unitPol[iPol][iMono], coor, L);
				}
				
				for(int iMono=0; iMono<3; iMono++)
					cs->topoState[cs->coorPol[iPol][iMono]] = topoState[iMono];
				iPol++;
			}
		}
	}
}



void GenerateMutators(LookupTables* lt, char* file){
	int nMoves[16][16];
	int mutIdTable[16][16];
	
	for(int i=0; i<16; i++){
		for(int j=0; j<16; j++){
			nMoves[i][j]=0;
			mutIdTable[i][j]=-1;
			for(int k=0; k<4; k++){
				for(int l=0; l<3; l++){
					lt->mutator[i][j][k][l]=-1;
				}
				for(int l=0; l<2; l++)
					lt->newUnits[i][j][k][l]=-1;
				for(int l=0; l<4; l++)
					lt->counts[i][j][k][l]=0;

			}
		}
	}
	
	int mutId=0;
	for(int unitA=1; unitA<0xf; unitA++){
		if(!IsValid(unitA)) continue;
		for(int unitB=1; unitB<0xf; unitB++){
			if(!IsValid(unitB)) continue;
			int unitC;
			if(!ValidateAddUnitVectors(unitA, unitB, &unitC)) continue;
			
			mutIdTable[unitA][unitB] = mutId++;
// 			printf("%x -> %x = %i\n", unitA, unitB, mutId-1);
		}
	}
	
	for(int unitA=1; unitA<0xf; unitA++){
		if(!IsValid(unitA)) continue;
		for(int unitB=1; unitB<0xf; unitB++){
			if(!IsValid(unitB)) continue;
			int unitC;
			if(!ValidateAddUnitVectors(unitA, unitB, &unitC)) continue;
			lt->mutator [unitA][unitB][nMoves[unitA][unitB]][0] = 2*mutIdTable[(~unitA)&0xf][unitC];
			lt->mutator [unitA][unitB][nMoves[unitA][unitB]][1] = 2*mutIdTable[unitA][unitB]+1;
			lt->mutator [unitA][unitB][nMoves[unitA][unitB]][2] = 2*mutIdTable[unitB][(~unitC)&0xf];
			lt->newUnits[unitA][unitB][nMoves[unitA][unitB]][0] = 0;
			lt->newUnits[unitA][unitB][nMoves[unitA][unitB]][1] = unitC;
			nMoves[unitA][unitB]++;
			
			lt->mutator [unitA][unitB][nMoves[unitA][unitB]][0] = 2*mutIdTable[(~unitA)&0xf][unitC];
			lt->mutator [unitA][unitB][nMoves[unitA][unitB]][1] = 2*mutIdTable[unitA][unitB]+1;
			lt->mutator [unitA][unitB][nMoves[unitA][unitB]][2] = 2*mutIdTable[unitB][(~unitC)&0xf];
			lt->newUnits[unitA][unitB][nMoves[unitA][unitB]][0] = unitC;
			lt->newUnits[unitA][unitB][nMoves[unitA][unitB]][1] = 0;
			nMoves[unitA][unitB]++;
			
			lt->mutator [0][unitC][nMoves[0][unitC]][0] = 2*mutIdTable[(~unitC)&0xf][unitA];
			lt->mutator [0][unitC][nMoves[0][unitC]][1] = 2*mutIdTable[unitA][unitB]+1;
			lt->mutator [0][unitC][nMoves[0][unitC]][2] = 2*mutIdTable[unitC][(~unitB)&0xf];
			lt->newUnits[0][unitC][nMoves[0][unitC]][0] = unitA;
			lt->newUnits[0][unitC][nMoves[0][unitC]][1] = unitB;
			nMoves[0][unitC]++;
			
			lt->mutator [unitC][0][nMoves[unitC][0]][0] = 2*mutIdTable[(~unitC)&0xf][unitA];
			lt->mutator [unitC][0][nMoves[unitC][0]][1] = 2*mutIdTable[unitA][unitB]+1;
			lt->mutator [unitC][0][nMoves[unitC][0]][2] = 2*mutIdTable[unitC][(~unitB)&0xf];
			lt->newUnits[unitC][0][nMoves[unitC][0]][0] = unitA;
			lt->newUnits[unitC][0][nMoves[unitC][0]][1] = unitB;
			nMoves[unitC][0]++;
		}
	}
	
// 	for(int unitA=1; unitA<0xf; unitA++){
// 		if(!IsValid(unitA)) continue;
// 		for(int unitB=1; unitB<0xf; unitB++){
// 			if(!IsValid(unitB)) continue;
// 			if(nMoves[unitA][unitB] == 2){
// 				for(int k=0; k<2; k++){
// 					for(int i=0; i<3; i++)
// 						lt->mutator[unitA][unitB][2+k][i] = lt->mutator[unitA][unitB][k][i];
// 					for(int i=0; i<2; i++)
// 						lt->newUnits[unitA][unitB][2+k][i] = lt->newUnits[unitA][unitB][k][i];
// 				}
// 			}
// 		}
// 	}
// 	
	TopoMapFromFile(lt,file);
	
}
