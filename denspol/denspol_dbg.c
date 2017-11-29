#include "denspol_dbg.h"
int CheckIntegrity(CurState* cs, char* msg){
	int L = cs->L;
	for(int iPol=0; iPol<cs->nPol; iPol++){
		int coor = cs->coorPol[iPol][0];
		for(int iMono=0; iMono<cs->polSize; iMono++){
			coor= AddUnitToCoor(cs->unitPol[iPol][iMono], coor, L);
#if POL_TYPE == POL_TYPE_LIN
			int jMono = iMono+1;
#elif POL_TYPE == POL_TYPE_RING
			int jMono = (iMono+1)%cs->polSize;
#endif
			if(coor != cs->coorPol[iPol][jMono]){
				printf("Error: coor[%i]+unit[%i] != coor[%i]\n", iMono, iMono, iMono+1);
// 				printf("%i +%i != %i\n", cs->coorPol
				printf("At: %s\n", msg);
				exit(192);
			}
#if POL_TYPE == POL_TYPE_RING
			if(cs->unitPol[iPol][iMono]+cs->unitPol[iPol][(iMono+1)%cs->polSize] == 0xf){
				printf("Error in units: hairpin configuration\n");
				printf("At: %s\n", msg);
				exit(192);
			}
#elif POL_TYPE == POL_TYPE_LIN
			if(iMono<cs->polSize-1 && cs->unitPol[iPol][iMono]+cs->unitPol[iPol][(iMono+1)%cs->polSize] == 0xf){
				printf("Error in units: hairpin configuration\n");
				printf("At: %s\n", msg);
				exit(192);
			}
#endif
// #if TOPO_DENSE == TRUE
// #if POL_TYPE == POL_TYPE_LIN
// 			if(
			if(cs->unitPol[iPol][iMono]){
				int jMono=(iMono+1)%cs->polSize;
				while(!cs->unitPol[iPol][jMono]){
					jMono = (jMono+1)%cs->polSize;
				}
				if(!(cs->bondOcc[cs->coorPol[iPol][iMono]]&(1<<(cs->unitPol[iPol][iMono]^0xf)))){
					printf("Bond not stored!\n");
					printf("At: %s\n", msg);
					printf("iMono = %i, polSize=%i\n", iMono, cs->polSize);
					printf("bondOcc: %x, units[%x,%x]\n", cs->bondOcc[coor], cs->unitPol[iPol][iMono], cs->unitPol[iPol][jMono]);
					printf("iMono = %i, iPol = %i, coor = (", (iMono+1)%cs->polSize, iPol); PrintCoor(cs->coorPol[iPol][(iMono-1+cs->polSize)%cs->polSize],L); printf(", "); PrintCoor(cs->coorPol[iPol][iMono],L); printf(", "); PrintCoor(coor,L); printf(", "); PrintCoor(cs->coorPol[iPol][(iMono+2)%cs->polSize], L); printf(")\n");
					exit(192);
				}
#if POL_TYPE == POL_TYPE_LIN
				if(jMono <= iMono) continue;
#endif
				if(!(cs->bondOcc[cs->coorPol[iPol][jMono]]&(1<<cs->unitPol[iPol][iMono]))){
					printf("Bond not stored 2\n");
					printf("iPol = %i, iMono = %i, jMono = %i, coor = %i\n", iPol, iMono, jMono, cs->coorPol[iPol][jMono]);
					PrintCoor(coor, L); printf("\n");
					printf("BondOcc=%x, unit=%x, topo=%x\n", cs->bondOcc[cs->coorPol[iPol][jMono]], cs->unitPol[iPol][iMono], cs->topoState[cs->coorPol[iPol][jMono]]);
					printf("At: %s\n", msg);
					exit(192);
				}
			}
// #endif
		}
	}
// #if TOPO_DENSE == TRUE
#if POL_TYPE == POL_TYPE_LIN
	int cycle = 2*cs->polSize;
#elif POL_TYPE == POL_TYPE_RING
	int cycle = cs->polSize;
#endif

	int nTrueBonds=0;
	for(int iPol=0; iPol<cs->nPol; iPol++){
		for(int iMono=0; iMono<cs->polSize; iMono++){
			int iCoor1= cs->coorPol[iPol][iMono];
			int iCoor2= cs->coorPol[iPol][(iMono+1)%cycle];
			if(iCoor1 == iCoor2) continue;
			nTrueBonds++;
			for(int jPol=0; jPol<cs->nPol; jPol++){
				for(int jMono=0; jMono<cs->polSize; jMono++){
					if(iPol == jPol && iMono == jMono) continue;
					int jCoor1 = cs->coorPol[jPol][jMono];
					int jCoor2 = cs->coorPol[jPol][(jMono+1)%cycle];
					if(jCoor1 == jCoor2) continue;
					if((iCoor1 == jCoor1 && iCoor2 == jCoor2) || (iCoor1 == jCoor2 && iCoor2 == jCoor1)){
						printf("Error: bond occupied!\n");
						printf("coor: (%i %i) vs (%i %i)\n", iCoor1, iCoor2, jCoor1, jCoor2);
						printf("(%i, %i) vs (%i, %i)\n", iPol, iMono, jPol, jMono);
						printf("At: %s\n", msg);
						exit(192);
					}
				}
			}
		}
	}
	
	int nBondOcc=0;
	for(int i=0; i<cs->LSize; i++){
		for(int j=0; j<16; j++){
			if(cs->bondOcc[i]&(1<<j)) nBondOcc++;
		}
	}
	
	if(nBondOcc != 2*nTrueBonds){
		printf("Error: Number of bonds is not equal to the true number of bonds: %i vs %i\n", nBondOcc, nTrueBonds);
		printf("At: %s\n", msg);
		exit(192);
	}
// #endif
	return 0;
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

void PrintUnit(int unit){
	char ret[5];
	int len=0;
	
	if(unit&0x1) ret[len++]='t';
	if(unit&0x2) ret[len++]='u';
	if(unit&0x4) ret[len++]='v';
	if(unit&0x8) ret[len++]='w';
	if(!unit) sprintf(ret, "0");
	else ret[len++] = '\0';
	printf("%s", ret);
}

void PrintUnitDot(LookupTables* lt){
	for(int i=0; i<16; i++){
		if(!IsValid(i)) continue;
		for(int j=0; j<16; j++){
			if(!IsValid(j)) continue;
			PrintUnit(i); printf(" * "); PrintUnit(j); 
			printf(" = %i\n", lt->unitDot[i][j]);
		}
	}
	printf("\n");
}



void ComputeSymmetry(int table[12], int depth, LookupTables* lt){
	
	if(depth == 12){
		for(int i=0; i<12; i++){
			PrintUnit(lt->validUnits[i]); printf(" -> "); PrintUnit(lt->validUnits[table[i]]);
			printf("\n");
		}
		printf("\n");
		return;
	}
	
	int oldUnit = lt->validUnits[depth];
	for(int i=0; i<12; i++){
		int newUnit = lt->validUnits[i];
		int valid=1;
		for(int j=0; j<depth && valid; j++){
			int oldCor = lt->unitDot[oldUnit][lt->validUnits[j]];
			int newCor = lt->unitDot[newUnit][lt->validUnits[table[j]]];
			if(newCor != oldCor) valid=0;
		}
		if(valid){
			table[depth] = i;
			ComputeSymmetry(table, depth+1, lt);
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


