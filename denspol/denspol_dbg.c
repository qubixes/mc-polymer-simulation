#include "denspol_dbg.h"
void CheckIntegrity(CurState* cs, char* msg){
	int L = cs->L;
	for(int iPol=0; iPol<cs->nPol; iPol++){
		int coor = cs->coorPol[iPol][0];
		for(int iMono=0; iMono<cs->polSize; iMono++){
			coor= AddUnitToCoor(cs->unitPol[iPol][iMono], coor, L);
			if(coor != cs->coorPol[iPol][(iMono+1)%cs->polSize]){
				printf("Error in coordinates/unit conflageration\n");
				exit(192);
			}
			if(cs->unitPol[iPol][iMono]+cs->unitPol[iPol][(iMono+1)%cs->polSize] == 0xf){
				printf("Error in units: hairpin configuration\n");
				printf("At: %s\n", msg);
				exit(192);
			}
		}
	}
}

void PrintSystem(CurState* cs){
	for(int iPol=0; iPol<cs->nPol; iPol++){
		for(int iMono=0; iMono<cs->polSize; iMono++){
			printf("%x", cs->unitPol[iPol][iMono]);
		}
		
		for(int iMono=0; iMono<cs->polSize; iMono++){
			printf(" %i ", cs->topoState[cs->coorPol[iPol][iMono]]);
		}
		printf("\n");
	}
}

void PrintMutators(LookupTables* lt){
	for(int unitA=0; unitA<0xf; unitA++){
		for(int unitB=0; unitB<0xf; unitB++){
			for(int i=0; i<4; i++){
				int* mut = lt->mutator[unitA][unitB][i];
				int* newUnits = lt->newUnits[unitA][unitB][i];
				if(mut[0] < 0){if(i>0) printf("\n"); break; }
				if(i==0) printf("[%x %x]: ", unitA, unitB);
				printf(" [%x %x] (%i %i %i)", newUnits[0], newUnits[1], mut[0], mut[1], mut[2]);
				if(i==3) printf("\n");
			}
		}
	}
	printf("\n");
}

void PrintMoveCounts(LookupTables* lt){
	for(int unitA=0; unitA<0xf; unitA++){
		if(!IsValid(unitA) && unitA) continue;
		for(int unitB=0; unitB<0xf; unitB++){
			if(!IsValid(unitB) && unitB) continue;
			for(int rand=0; rand<4; rand++){
				if(lt->counts[unitA][unitB][rand][0]){
					printf("%x + %x => %x + %x: ", unitA, unitB, lt->newUnits[unitA][unitB][rand][0], lt->newUnits[unitA][unitB][rand][1]);
					for(int k=0; k<4; k++) printf("%i ", lt->counts[unitA][unitB][rand][k]);
					printf("\n");
				}
			}
		}
	}
}

void StatTopo(LookupTables* lt){
	int depth=0, nextDepthStart=0, maxNext=-1;
// 	int maxDepth=6;
	int maxDepth=11;
	
	int counts[maxDepth+1][NMUTATOR];
	for(int depth=0; depth<=maxDepth; depth++){
		for(int i=0; i<NMUTATOR; i++) counts[depth][i]=0;
	}
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		for(int i=0; i<NMUTATOR; i++){
			if(lt->mutTopo[iTopo][i]>=0){
				counts[depth][i]++;
				if(lt->mutTopo[iTopo][i] > maxNext) maxNext = lt->mutTopo[iTopo][i];
			}
		}
		if(iTopo == nextDepthStart){
			depth++;
			nextDepthStart=maxNext;
		}
	}
	
	
	for(int i=0; i<NMUTATOR; i++){
		printf("#[%i] = ", i);
		for(int depth=0; depth<=maxDepth; depth++){
			printf("%i ", counts[depth][i]);
		}
		printf("\n");
	}
}

void CheckTopo(LookupTables* lt){
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			int newTopo = lt->mutTopo[iTopo][iMut];
			if(newTopo>=0){
				int balance=0;
				for(int jMut=0; jMut<NMUTATOR; jMut++){
					int backTopo = lt->mutTopo[newTopo][jMut];
					if(backTopo == iTopo){ balance=1; break; }
				}
				if(!balance){
					printf("Error: no detailed balance!\n");
					printf("%i + %i -> %i??\n", iTopo, iMut, newTopo);
					for(int jMut=0; jMut<NMUTATOR; jMut++){
						int backTopo = lt->mutTopo[newTopo][jMut];
						if(backTopo >= 0){
							printf("%i + %i = %i\n", newTopo, jMut, backTopo);
						}
					}
					exit(0);
				}
			}
		}
	}
}
