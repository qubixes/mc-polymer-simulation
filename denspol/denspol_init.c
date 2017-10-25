#include "denspol.h"


void SimulationInit(CurState* cs, LookupTables* lt, unsigned int seed, double density, int polSize, int L, char* dir){
	long lastT = GetLastT(dir);
	
	UnitDotInit(lt, cs->ss.bendEnergy);
	GenerateMutators(lt, cs->ss.eeFile);
	SuperTopoInit(lt);
	lt->hp->distance = GenerateDistanceMatrix(L);
	
	if(lastT<0){
		cs->curT=0;
		int nPol = (int)(L*L*L*density/(double)polSize+0.5);
		CSInit(cs, seed, nPol, polSize, L, dir);
		GeneratePolymers(cs, lt);
	}
	else{
		printf("Starting from t=%li\n", lastT);
		cs->curT=lastT;
		CSFromFile(cs, dir, lastT);
	}
	CheckIntegrity(cs, "After construction");
}

void LoadHPFile(char* file, HPTable* hp, CurState* cs){
	int polSize = cs->polSize;
	FILE* pFile = fopen(file, "r");
	int nMono;
	fscanf(pFile, "%*s %i", &nMono);
	double ratio = (polSize+1)/(double)nMono;
	
	
	hp->HPStrength = malloc(sizeof(double*)*(polSize+1));
	hp->monoId     = malloc(sizeof(int*)*(polSize+1));
	hp->nInterHP   = malloc(sizeof(int)*(polSize+1));
	for(int iPol=0; iPol<=polSize; iPol++){
		hp->HPStrength[iPol] = malloc(sizeof(double)*(polSize+1));
		hp->monoId[iPol]     = malloc(sizeof(int)*(polSize+1));
		hp->nInterHP[iPol]   = 0;
	}
	int iMono, jMono, iPol, jPol;
	double strength;
	while( !(feof(pFile)) && fscanf(pFile, "%i %i %i %i %lf", &iMono, &iPol, &jMono, &jPol, &strength) == 3){
		int iMonoId = MonoPol2Id((int)(iMono*ratio), iPol, cs);
		int jMonoId = MonoPol2Id((int)(jMono*ratio), jPol, cs);
		
		hp->monoId[iMonoId][hp->nInterHP[iMonoId]] = jMonoId;
		hp->HPStrength[iMonoId][hp->nInterHP[iMonoId]++] = strength;
		hp->monoId[jMonoId][hp->nInterHP[jMonoId]] = iMonoId;
		hp->HPStrength[jMonoId][hp->nInterHP[jMonoId]++] = strength;
	}
}



void CSInit(CurState* cs, unsigned int seed, int nPol, int polSize, int L, char* dir){
	char exec[1000];
	
	cs->ss.dir = dir;
	cs->L = L;
	cs->LSize = L*L*L;
	cs->nPol = nPol;
	cs->polSize = polSize;
	
	cs->topoState = malloc(sizeof(int)*cs->LSize);
	cs->bondOcc = malloc(sizeof(int)*cs->LSize);
	cs->unitPol = malloc(sizeof(int*)*cs->nPol);
	cs->coorPol = malloc(sizeof(int*)*cs->nPol);
	
	for(int iPol=0; iPol<cs->nPol; iPol++){
		cs->unitPol[iPol] = malloc(sizeof(int)*(cs->polSize+1));
		cs->coorPol[iPol] = malloc(sizeof(int)*(cs->polSize+1));
	}
	
	for(int i=0; i<cs->LSize; i++){
		cs->topoState[i]=0;
		cs->bondOcc[i]=0;
	}
	
	Seed(cs->rngState, seed);
	
	sprintf(exec, "mkdir -p %s", dir);
	system(exec);
}




void GenerateRingPolymers(CurState* cs, LookupTables* lt){
	int delta;
	int L = cs->L;
	
	for(delta=L; delta>=1; delta--){
		int nTot = 1;
		nTot *= L/delta;
		nTot *= L/delta;
		nTot *= L/delta;
		if(nTot >= cs->nPol){
			break;
		}
	}
	if(!delta){
		printf("Failed to find partition for polymers\n");
		exit(102);
	}
	
	
	int units[3] = {0x4,0x1,0xa};
	int topoState[3];
	
	for(int i=0; i<3; i++){
		int j;
		for(j=0; j<4; j++){
			if(lt->newUnits[(~units[i])&0xf][0][j][0] != ((units[(i+1)%3])&0xf)) continue;
			if(lt->newUnits[(~units[i])&0xf][0][j][1] != ((units[(i+2)%3])&0xf)) continue;	
			break;
		}
		topoState[(i+2)%3] = lt->mutTopo[0][lt->mutator[(~units[i])&0xf][0][j][1]];
	}
	int iPol=0;
	for(int t=0; t+delta/2<L && iPol<cs->nPol; t += delta){
		for(int u=0; u+delta/2<L && iPol<cs->nPol; u += delta){
			for(int v=0; v+delta/2<L && iPol<cs->nPol; v += delta){
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
				
				for(int iMono=0; iMono<3; iMono++){
// #if TOPO_DENSE == FALSE
// 					cs->topoState[cs->coorPol[iPol][iMono]] = topoState[iMono];
// #endif
					int bondOcc = 1<<cs->unitPol[iPol][(iMono-1+3)%3];
					bondOcc |= 1<<(cs->unitPol[iPol][iMono]^0xf);
// 					if(cs->coorPol[iPol][iMono] == 0){
// 						printf("iPol = %i\n", iPol);
// 						printf("(%i,%i,%i)\n", t,u,v);
// 						printf("%i: (1<<(%i, %i))\n", iMono, cs->unitPol[iPol][(iMono-1+3)%3], cs->unitPol[iPol][iMono]^0xf);
// 						
// 					}
					cs->bondOcc[cs->coorPol[iPol][iMono]] |= bondOcc;
				}
				iPol++;
			}
		}
	}
}

void GenerateLinearPolymers(CurState* cs, LookupTables* lt){
	int delta;
	int L = cs->L;
	
	for(delta=L; delta>=1; delta--){
		int nTot = 1;
		nTot *= L/delta;
		nTot *= L/delta;
		nTot *= L/delta;
		if(nTot >= cs->nPol){
			break;
		}
	}
	if(!delta){
		printf("Failed to find partition for polymers\n");
		exit(102);
	}
	
	int iPol=0;
	for(int t=0; t+delta/2<L && iPol<cs->nPol; t += delta){
		for(int u=0; u+delta/2<L && iPol<cs->nPol; u += delta){
			for(int v=0; v+delta/2<L && iPol<cs->nPol; v += delta){
				cs->unitPol[iPol][0] = 0x1;
				
				for(int iMono=1; iMono<cs->polSize; iMono++)
					cs->unitPol[iPol][iMono] = 0;
				cs->unitPol[iPol][cs->polSize] = 0xf;
				
				int coor = TUV2Coor(t,u,v,L);
				for(int iMono=0; iMono<=cs->polSize; iMono++){
					cs->coorPol[iPol][iMono] = coor;
					coor = AddUnitToCoor(cs->unitPol[iPol][iMono], coor, L);
				}
				
				int bondOcc0 = 1<<(cs->unitPol[iPol][0]^0xf);
				int bondOcc1 = 1<<(cs->unitPol[iPol][0]);
				
				cs->bondOcc[cs->coorPol[iPol][0]] |= bondOcc0;
				cs->bondOcc[cs->coorPol[iPol][1]] |= bondOcc1;
				
				iPol++;
			}
		}
	}
}

void GeneratePolymers(CurState* cs, LookupTables* lt){
#if POL_TYPE == POL_TYPE_LIN
	GenerateLinearPolymers(cs,lt);
#elif POL_TYPE == POL_TYPE_RING
	GenerateRingPolymers(cs,lt);
#endif
}

void GenerateMutators(LookupTables* lt, char* file){
	int nMoves[16][16];
	int mutIdTable[16][16];
	
	for(int i=0; i<16; i++){
		for(int j=0; j<16; j++){
			nMoves[i][j]=0;
			mutIdTable[i][j]=-1;
			lt->mutIdTableDouble[i][j]=NON_EXISTING;
			for(int k=0; k<16; k++)
				lt->mutIdTableTriple[i][j][k]=NON_EXISTING;
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
// 	int mutIdTriple=0;
	for(int unitA=1; unitA<0xf; unitA++){
		if(!IsValid(unitA)) continue;
		for(int unitB=1; unitB<0xf; unitB++){
			if(!IsValid(unitB)) continue;
			int unitAB;
			if(!ValidateAddUnitVectors(unitA, unitB, &unitAB)) continue;
			
			lt->mutIdTableDouble[unitA][unitB] = mutId;
			lt->revMutTable[mutId][0] = unitA;
			lt->revMutTable[mutId][1] = unitB;
			mutIdTable[unitA][unitB] = mutId++;
			
			int nTriple=0;
			for(int unitC=1; unitC<0xf; unitC++){
				if(!IsValid(unitC)) continue;
				int unitBC;
				if(!ValidateAddUnitVectors(unitB, unitC, &unitBC)) continue;
				int unitABC;
				if(unitC == ((~unitAB)&0xf)){
// 					printf("[NE]: ");
// 					PrintUnit(unitA); printf(" => "); PrintUnit(unitB); printf(" => "); PrintUnit(unitBC); printf("\n");

// 					printf("[Non-existing]: %x => %x => %x\n", unitA, unitB, unitC);
					continue;
				}
// 				printf("%x -> %x = %x\n", unitA, unitB, unitC);
				if(ValidateAddUnitVectors(unitA, unitBC, &unitABC)){
					lt->mutIdTableTriple[unitA][unitB][unitBC] = SAME_TOPO;
// 					printf("[ST]: ");
// 					PrintUnit(unitA); printf(" => "); PrintUnit(unitB); printf(" => "); PrintUnit(unitBC); printf("\n");
					continue;
				}
				lt->mutIdTableTriple[unitA][unitB][unitBC] = (mutId-1)+48*nTriple;
				
				lt->revMutTableTriple[(mutId-1)+48*nTriple][0] = unitA;
				lt->revMutTableTriple[(mutId-1)+48*nTriple][1] = unitB;
				lt->revMutTableTriple[(mutId-1)+48*nTriple][2] = unitBC;
				
				nTriple++;
// 				PrintUnit(unitA); printf(" => "); PrintUnit(unitB); printf(" => "); PrintUnit(unitBC); printf("\n");
// 				printf("%x => %x => %x\n", unitA, unitB, unitBC);
			}
		}
	}
// 	printf("nDouble = %i, nTriple = %i\n", mutId, mutIdTriple);
	
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
	
	for(int unitA=1; unitA<0xf; unitA++){
		if(!IsValid(unitA)) continue;
		for(int unitB=1; unitB<0xf; unitB++){
			if(!IsValid(unitB)) continue;
			if(nMoves[unitA][unitB] == 2){
				for(int k=0; k<2; k++){
					for(int i=0; i<3; i++)
						lt->mutator [unitA][unitB][2+k][i] = lt->mutator [unitA][unitB][k][i];
					for(int i=0; i<2; i++)
						lt->newUnits[unitA][unitB][2+k][i] = lt->newUnits[unitA][unitB][k][i];
				}
			}
		}
	}
	
	TopoMapFromFile(lt,file);
}

void UnitDotInit(LookupTables* lt, double bendEnergy){
	
	int nValid=0;
	for(int i=0; i<16; i++){
		if(IsValid(i)){
			lt->validUnits[nValid++]=i;
		}
	}
	
	for(int i=0; i<16; i++){
		for(int j=0; j<16; j++){
			if(IsValid(i) && IsValid(j)){
				lt->unitDot[i][j] = UnitInProd(i,j);
			}
			else
				lt->unitDot[i][j] = -100;
		}
	}
	for(int dBend=-BEND_LVL; dBend<BEND_LVL; dBend++){
		if(dBend < 0)
			lt->bendProb[dBend+BEND_LVL] = exp(bendEnergy*dBend);
		else
			lt->bendProb[dBend+BEND_LVL] = 1;
	}
}

SuperTopo* NewSuperTopo(Topo* topo){
	SuperTopo* sTopo = malloc(sizeof(SuperTopo));
	sTopo->topo = topo;
	sTopo->id = NON_EXISTING;
	sTopo->permBond=NON_EXISTING;
	for(int i=0; i<NMUTATOR; i++) sTopo->mutTopo[i] = NULL;
	return sTopo;
}

void ClearTopo(Topo* topo){
	topo->nBonds=0;
	for(int i=0; i<NMUTATOR; i++) topo->mutTopo[i]=NULL;
	topo->next = NULL;
	topo->supTopo = NewSuperTopo(topo);
	topo->set=0;
}

Topo* NewTopo(){
	Topo* topo = malloc(sizeof(Topo));
	topo->nBonds=0;
	for(int i=0; i<NMUTATOR; i++) topo->mutTopo[i]=NULL;
	topo->next = NULL;
	topo->set = 0;
	topo->supTopo = NewSuperTopo(topo);
	return topo;
}

void CopyTopo(Topo* src, Topo* dst){
	
	dst->nBonds=src->nBonds;
	for(int i=0; i<src->nBonds; i++){
		for(int k=0; k<2; k++){
			dst->bonds[i][k] = src->bonds[i][k];
		}
		dst->bondFree[i] = src->bondFree[i];
	}
	
}

void MergeTopo(Topo* topoA, Topo* topoB){
	if(!(topoA && topoB)) return;
	
	SuperTopo *sTopoA = topoA->supTopo;
	SuperTopo *sTopoB = topoB->supTopo;
	
	if(sTopoA == sTopoB) return;
	
	Topo* curTopo=sTopoA->topo;
	
	while(curTopo->next){
		curTopo->supTopo = sTopoB;
		curTopo = curTopo->next;
	}
	
	curTopo->supTopo = sTopoB;
	curTopo->next = sTopoB->topo;
	
	sTopoB->topo = sTopoA->topo;
	
	free(sTopoA);
	nSupTopo--;
}

void PrintTopo(Topo* topo){
	if(!topo->nBonds) printf("[empty]\n");
	else{
		for(int iBond=0; iBond<topo->nBonds; iBond++){
			printf("["); PrintUnit(topo->bonds[iBond][0]); printf(" => "); PrintUnit(topo->bonds[iBond][1]); printf("]\n");
		}
	}
	printf("\n");
}

void PrintSTopo(SuperTopo* sTopo){
	
	int minBonds=999;
	Topo* topoMin = NULL;
	
	Topo* topo = sTopo->topo;
	
	while(topo){
		if(topo->nBonds <minBonds){
			topoMin = topo;
			minBonds = topo->nBonds;
		}
		topo = topo->next;
	}
	PrintTopo(topoMin);
}

void MutateTopo(Topo* topo, int iMut, Topo* destTopo, LookupTables* lt){
	
	int *newBond = lt->revMutTable[iMut/2];
	CopyTopo(topo, destTopo); 
	
	if(iMut%2 == 1){
		int exists=0;
		for(int iBond=0; iBond<topo->nBonds; iBond++){
			if( destTopo->bonds[iBond][0] == newBond[0] && destTopo->bonds[iBond][1] == newBond[1] ){
				destTopo->nBonds--;
				for(int k=0; k<2; k++)
					destTopo->bonds[iBond][k] = destTopo->bonds[destTopo->nBonds][k];
				exists = 1;
			}
		}
		
		if(!exists){
			for(int k=0; k<2; k++)
				destTopo->bonds[destTopo->nBonds][k] = newBond[k];
			destTopo->nBonds++;
		}
	}
	else{
		int inserted=0;
		for(int iBond=0; iBond<topo->nBonds && ! inserted; iBond++){
			for(int k=0; k<2 && !inserted; k++){
				if(destTopo->bonds[iBond][k] == (newBond[0]^(k*0xf)) ){
					destTopo->bonds[iBond][k] = newBond[1]^((0x1^k)*0xf);
					inserted = 1;
				}
			}
		}
		if(!inserted){
			printf("Error: mutation failed to insert bond\n");
			PrintTopo(topo); 
			PrintTopo(destTopo);
			PrintUnit(newBond[0]); printf(" --> "); PrintUnit(newBond[1]); printf("\n");
			exit(0);
		}
	}
}

void UpdateTopoMutations(Topo* topo, LookupTables* lt){
	
	int permBond=0;
	for(int iBond=0; iBond<topo->nBonds; iBond++){
		int mut = lt->mutIdTableDouble[topo->bonds[iBond][0]][topo->bonds[iBond][1]];
		if(mut>=0 && topo->mutTopo[2*mut+1]) continue;
		permBond |= 1 << topo->bonds[iBond][0];
		permBond |= 1 << (topo->bonds[iBond][1]^0xf);
	}
	
	if(topo->supTopo->permBond >= 0 && topo->supTopo->permBond != permBond){
		printf("Making a mistake while figuring out the perma bonds\n");
		printf("%x vs %x\n", topo->supTopo->permBond, permBond);
		PrintTopo(topo);
		exit(0);
	}
	topo->supTopo->permBond = permBond;
	
	
	for(int iMut=0; iMut<NMUTATOR; iMut++){
		if(!topo->mutTopo[iMut]) continue;
		
		Topo* destTopo = topo->mutTopo[iMut];
		if(!destTopo->set){
			MutateTopo(topo, iMut, destTopo, lt);
		}

		if(iMut%2 == 1) /// if it only adds/removes a link, it does not alter the super topo.
			continue;
		
		int start = lt->revMutTable[iMut/2][0];
		int end = lt->revMutTable[iMut/2][1];
		
		int foundBond=0;
		for(int iBond=0; iBond<topo->nBonds && !foundBond; iBond++){
			int oldMut = lt->mutIdTableDouble[topo->bonds[iBond][0]][topo->bonds[iBond][1]];
			for(int k=0; k<2 && !foundBond; k++){
				if((topo->bonds[iBond][k]^(k*0xf)) == start){
					int mutation=-1;
					if(oldMut>=0 && topo->mutTopo[2*oldMut+1]){
						if(k==0)
							mutation = lt->mutIdTableTriple[end^0xf][start^0xf][topo->bonds[iBond][1]];
						if(k==1)
							mutation = lt->mutIdTableTriple[topo->bonds[iBond][0]][start^0xf][end];
					}
					else{
						mutation = iMut/2; 
					}
					if(mutation == NON_EXISTING){
						printf("Uh oh! Not a valid mutation!\n");
						PrintUnit(topo->bonds[iBond][0]); printf(" -> "); PrintUnit(topo->bonds[iBond][1]); printf(" + "); PrintUnit(start); printf(" -> "); PrintUnit(end);
						printf("\nk=%i\n", k);
						exit(0);
					}
					else if (mutation != SAME_TOPO){
						if(topo->supTopo->mutTopo[mutation] && topo->supTopo->mutTopo[mutation]->supTopo != destTopo->supTopo){
							printf("Overwriting mutation!\n");
							exit(192);
						}
						topo->supTopo->mutTopo[mutation] = destTopo;
					}
					foundBond=1;
				}
			}
		}
	}
	
}

void PrintPermBond(int permBond){
	printf("(");
	for(int iB=0; iB<16; iB++){
		if(permBond & (1<<iB)){
			PrintUnit(iB); printf(" ");
		}
	}
	printf(")");
}

void PrintStraightTopo(LookupTables* lt){
	FILE* pFile = fopen("straight_topo.dat", "w");
	for(int i=0; i<lt->nTopoComp; i++){
		int permBond = lt->topComp[i].permBond;
		int bonds[2], nBonds=0;
		for(int iB=0; iB<16; iB++){
			if(permBond & (1<<iB)){
				if(nBonds>=2){
					nBonds=3; break;
				}
				else
					bonds[nBonds++] = iB;
			}
		}
		
		if(nBonds==2 && bonds[0] == 15-bonds[1]){
// 			fprintf(pFile, "%x %x %i\n", bonds[0], bonds[1], i);
			
			int found=0;
			for(int iMut=0; iMut<NMUTATOR && !found; iMut++){
				int topo1 = lt->topComp[i].mutators[iMut];
				if(topo1 < 0 ) continue;
				
				for(int jMut=0; jMut<NMUTATOR && !found; jMut++){
					int topo2 = lt->topComp[topo1].mutators[jMut];
					if(topo2 != 0) continue;
					
// 					PrintPermBond(lt->topComp[topo1].permBond);
// 					PrintUnit(bonds[0]); printf(" "); PrintUnit(bonds[1]); PrintPermBond(lt->topComp[i].permBond); printf(" => ");
// 					for(int k=0; k<2; k++){
// 						PrintUnit(lt->revMutTable[iMut][k]); printf(" ");
// 					}
// 					PrintPermBond(lt->topComp[topo1].permBond);
// 					printf("\n");
// 					for(int k=0; k<2; k++){
// 						PrintUnit(lt->revMutTable[jMut][k]); printf(" ");
// 					}
// 					printf("\n");
// 					printf("%i (%i)=> %i (%i)=> %i\n", i, iMut, topo1, jMut, topo2);
					
					int tripleMut=-1;
					for(int kMut=0; kMut<NMUTATOR; kMut++){
						if(lt->topComp[0].mutators[kMut] == topo1){
// 							for(int k=0; k<3; k++){
// 								PrintUnit(lt->revMutTableTriple[kMut][k]); printf(" ");
// 							}
							tripleMut = kMut;
						}
					}
// 					printf("\n\n");
					
					for(int k=0; k<2; k++){
						for(int l=0; l<3; l+=2){
							if(lt->revMutTableTriple[tripleMut][l] == bonds[k]){
								fprintf(pFile, "%x %i\n", bonds[k], i);
// 								printf("Correct bond = %i\n", bonds[k]);
								found=1;
							}
						}
					}
				}
			}
// 			printf("---------------------------\n");
// 			PrintUnit(bonds[0]); printf(" -> "); PrintUnit(bonds[1]); printf(" = %i\n", i); 
		}
	}
	exit(0);
}

void SuperTopoInit(LookupTables* lt){
	nSupTopo = lt->nTopo;
// 	printf("nTopo = %i\nnSupTopo = %i\n", lt->nTopo, nSupTopo);

	Topo* allTopo = malloc(sizeof(Topo)*lt->nTopo);
	
	for(int i=0; i<lt->nTopo; i++) ClearTopo(allTopo+i);
	
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(lt->mutTopo[iTopo][iMut] == -1) continue;
			allTopo[iTopo].mutTopo[iMut] = allTopo + lt->mutTopo[iTopo][iMut];
		}
	}
	
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(iMut%2 == 1){
				MergeTopo(allTopo+iTopo, allTopo[iTopo].mutTopo[iMut]);
			}
		}
	}
// 	printf("nTopo = %i\nnSupTopo = %i\n", lt->nTopo, nSupTopo);
	allTopo[0].set=1;
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		UpdateTopoMutations(allTopo+iTopo, lt);
	}
	
// 	for(int i=0; i<NMUTATOR; i++){
// 		if(allTopo[0].supTopo->mutSupTopo[i]){
// 			printf("%i\n", i);
// 		}
// 		else{
// 			printf("N/A\n");
// 		}
// 	}
	
	SuperTopo** sTopoArray = malloc(sizeof(SuperTopo*)*nSupTopo);
	lt->topComp = malloc(sizeof(SuperTopo)*nSupTopo);
	sTopoArray[0] = allTopo[0].supTopo;
	sTopoArray[0]->id = 0;
	for(int i=0; i<nSupTopo; i++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			lt->topComp[i].mutators[iMut] = NON_EXISTING;
		}
	}
	
	int nIds=1;
	for(int i=0; i<nSupTopo; i++){
		lt->topComp[i].permBond = sTopoArray[i]->permBond;
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(!sTopoArray[i]->mutTopo[iMut]) continue;
			SuperTopo* newSTopo= sTopoArray[i]->mutTopo[iMut]->supTopo;
			if(newSTopo->id <0 ){
				newSTopo->id = nIds;
				sTopoArray[nIds++] = newSTopo;
			}
			lt->topComp[i].mutators[iMut] = newSTopo->id;
		}
	}
	
	lt->nTopoComp = nIds;
	
	for(int i=0; i<nSupTopo; i++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			int newTopo = lt->topComp[i].mutators[iMut];
			if(newTopo < 0 ) continue;
			int foundBack=0;
			for(int jMut=0; jMut<NMUTATOR && !foundBack; jMut++){
				if(lt->topComp[newTopo].mutators[jMut] == i)
					foundBack=1;
			}
			if(!foundBack){
				PrintTopo(sTopoArray[i]->topo);
				for(int iMut=0; iMut<NMUTATOR; iMut++){
					if(sTopoArray[i]->mutTopo[iMut]){
						PrintSTopo(sTopoArray[i]->mutTopo[iMut]->supTopo);
					}
				}
				printf("==============================\n");
				PrintTopo(sTopoArray[newTopo]->topo);
				for(int iMut=0; iMut<NMUTATOR; iMut++){
					if(sTopoArray[newTopo]->mutTopo[iMut]){
						PrintSTopo(sTopoArray[newTopo]->mutTopo[iMut]->supTopo);
					}
				}
				printf("Error contructing topComp, no detailed balance\n");
				exit(192);
			}
		}
	}
}

double*** GenerateDistanceMatrix(int L){
	
	double*** distances = (double***)malloc(sizeof(double**)*(2*L)) + L;
	for(int i=-(L-1); i<=L-1; i++){
		distances[i] = (double**) malloc(sizeof(double*)*(2*L)) + L;
		for(int j=-(L-1); j<=L-1; j++){
			distances[i][j] = (double*) malloc(sizeof(double)*(2*L)) + L;
			for(int k=-(L-1); k<=L-1; k++){
				distances[i][j][k]=0;
			}
		}
	}
	
	int tuv[3], xyz[3], dist;
	
	for(int t=0; t<L; t++){
		for(int u=0; u<L; u++){
			for(int v=0; v<L; v++){
				int minDist=L*L*L*L;
				
				for(int dt=0; dt<=L; dt += L){
					for(int du=0; du<=L; du += L){
						for(int dv=0; dv<=L; dv += L){
							tuv[0] = t-dt;
							tuv[1] = u-du;
							tuv[2] = v-dv;
							
							TUV2XYZ(tuv, xyz);
							dist=0;
							for(int k=0; k<3; k++) dist += xyz[k];
							minDist = MIN(dist, minDist);
						}
					}
				}
				
				for(int dt=0; dt<=L; dt += L){
					for(int du=0; du<=L; du += L){
						for(int dv=0; dv<=L; dv += L){
							distances[t-dt][u-du][v-dv] = minDist/2.0;
						}
					}
				}
			}
		}
	}
	
	return distances;
}



