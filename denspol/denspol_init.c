#include "denspol.h"


void SimulationInit(CurState* cs, LookupTables* lt, unsigned int seed, double density, int polSize, int L, char* dir){
	long lastT = GetLastT(dir);
// #if POL_TYPE == POL_TYPE_LIN
// 	cs->polSize = polSize-1;
// 	cs->ss.poLSize = polSize-1;
// #else
	cs->polSize = polSize;
// #endif
	int nPol = (int)(L*L*L*density/(double)polSize+0.5);
	cs->nPol = nPol;
	if(cs->ss.hpFile){
		LoadHPFile(cs->ss.hpFile, &lt->hp, cs);
	}
	else{
		lt->hp.nInterHP = malloc(sizeof(int)*(polSize+1)*nPol);
		for(int i=0; i<(polSize+1)*nPol; i++) lt->hp.nInterHP[i]=0;
	}
	

	UnitDotInit(lt, cs->ss.bendEnergy);
	GenerateMutators(lt);
#if __FROM_TOPO_COMP__ == TRUE
	ReadTopComp(lt, cs->ss.eeFile);
#else
	TopoMapFromFile(lt, cs->ss.eeFile);
	SuperTopoInit(lt);
#endif
	
	lt->hp.distance = GenerateDistanceMatrix(L);
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
	if(cs->ss.hpFile) LoadHPFile(cs->ss.hpFile, &lt->hp, cs);
	if(cs->ss.polIdShuffle){
		printf("Shuffling polymer\n");
		ShufflePolymerIds(cs);
	}
	CheckIntegrity(cs, "After construction");
}

void AddInteraction(HPTable* hp, int iMono, int iPol, int jMono, int jPol, int nMono, double strength, CurState* cs){
// 	if(iPol != jPol) return;
	
#if POL_TYPE == POL_TYPE_RING
	double ratio = (cs->polSize)/(double)nMono;
#else
	double ratio = (cs->polSize+1)/(double)nMono;
#endif
	strength *= ratio*cs->ss.hpStrength;
	
	int iMonoId = MonoPol2Id((int)(iMono*ratio+0.5)%cs->polSize, iPol, cs->polSize);
	int jMonoId = MonoPol2Id((int)(jMono*ratio+0.5)%cs->polSize, jPol, cs->polSize);
// 	printf("iMonoId = %i, iMono = %i, ratio = %lf\n", iMonoId, iMono, ratio);
	
	if(iMonoId == jMonoId) return;
	int exist=0;
	for(int indexI=0; indexI<hp->nInterHP[iMonoId] && !exist; indexI++){
		if(hp->monoId[iMonoId][indexI] == jMonoId){
			hp->HPStrength[iMonoId][indexI] += strength;
			exist=1;
		}
	}
	if(!exist){
		hp->monoId[iMonoId][hp->nInterHP[iMonoId]] = jMonoId;
		hp->HPStrength[iMonoId][hp->nInterHP[iMonoId]++] = strength;
	}
	
	exist=0;
	for(int indexJ=0; indexJ<hp->nInterHP[jMonoId] && !exist; indexJ++){
		if(hp->monoId[jMonoId][indexJ] == iMonoId){
			hp->HPStrength[jMonoId][indexJ] += strength;
			exist=1;
		}
	}
	if(!exist){
		hp->monoId[jMonoId][hp->nInterHP[jMonoId]] = iMonoId;
		hp->HPStrength[jMonoId][hp->nInterHP[jMonoId]++] = strength;
	}
}

void LoadHPFile(char* file, HPTable* hp, CurState* cs){
	int polSize = cs->polSize;
	FILE* pFile = fopen(file, "r");
	int nMono;
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	fscanf(pFile, "%*s %i", &nMono);
	
	double ratio = (cs->polSize+1)/(double)nMono;
	int maxInter = MAX(12, 12/(ratio*ratio));

	
	hp->HPStrength = (double**)malloc(sizeof(double*)*(polSize+1)*cs->nPol);
	hp->monoId     = (int**)malloc(sizeof(int*)*(polSize+1)*cs->nPol);
	hp->nInterHP   = (int*)malloc(sizeof(int)*(polSize+1)*cs->nPol);
	for(int monoId=0; monoId<(polSize+1)*cs->nPol; monoId++){
		hp->HPStrength[monoId] = (double*)malloc(sizeof(double)*maxInter);
		hp->monoId[monoId]     = (int*)malloc(sizeof(int)*maxInter);
		hp->nInterHP[monoId]   = 0;
	}
	int iMono, jMono, iPol, jPol;
	double strength;
	while( !(feof(pFile)) && fscanf(pFile, "%i %i %i %i %lf", &iPol, &iMono, &jPol, &jMono, &strength) == 5){
		AddInteraction(hp, iMono, iPol, jMono, jPol, nMono, strength, cs);
	}
	int monoId;
	iPol=0;
	for(int iMono=0; iMono<cs->polSize; iMono++){
		monoId = MonoPol2Id(iMono, iPol, polSize);
		for(int iId=0; iId<hp->nInterHP[monoId]; iId++){
			int newMono, newPol;
			Id2MonoPol(hp->monoId[monoId][iId], &newMono, &newPol, polSize);
// 			double strength = hp->HPStrength[monoId][iId];
// 			printf("%i %i %i %lf\n", iMono, newPol, newMono, strength);
		}
	}
// 	exit(0);
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

/// Use a Fisher-Yates shuffle to randomize the polymer ID's.
/// 
/// Obviously, most of the time this is not what you want. 
/// One application is if you want to reconstruct a configuration, but not start
/// with polymers in the same position, which would be akin to cheating.

void ShufflePolymerIds(CurState* cs){
	int* tUnit, *tCoor;
	
	for(int i=0; i<cs->nPol-1; i++){
		int j= i+DRng(cs->rngState)*(cs->nPol-i);
		if(i != j){
			tUnit = cs->unitPol[i];
			cs->unitPol[i] = cs->unitPol[j];
			cs->unitPol[j] = tUnit;
			
			tCoor = cs->coorPol[i];
			cs->coorPol[i] = cs->coorPol[j];
			cs->coorPol[j] = tCoor;
		}
	}
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
	
	int iPol=0;
	for(int t=0; t+delta/2<L && iPol<cs->nPol; t += delta){
		for(int u=0; u+delta/2<L && iPol<cs->nPol; u += delta){
			for(int v=0; v+delta/2<L && iPol<cs->nPol; v += delta){
				for(int j=0; j<3; j++)
					cs->unitPol[iPol][j] = units[j];
				for(int iMono=3; iMono<cs->polSize; iMono++)
					cs->unitPol[iPol][iMono] = 0;
				int coor = TUV2Coor(t,u,v,L);
				for(int iMono=0; iMono<cs->polSize; iMono++){
					cs->coorPol[iPol][iMono] = coor;
					coor = AddUnitToCoor(cs->unitPol[iPol][iMono], coor, L);
				}
				
				for(int iMono=0; iMono<3; iMono++){
					int bondOcc = 1<<cs->unitPol[iPol][(iMono-1+3)%3];
					bondOcc |= 1<<(cs->unitPol[iPol][iMono]^0xf);
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

void GenerateMutators(LookupTables* lt){
	int nMoves[16][16];
	int mutIdTable[16][16];
	
	lt->sameTopoTree = NULL;
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


double*** GenerateDistanceMatrix(int L){
	double*** distances = (double***)malloc(sizeof(double**)*(2*L)) + L;
	for(int i=-L; i<=L-1; i++){
		distances[i] = (double**) malloc(sizeof(double*)*(2*L)) + L;
		for(int j=-L; j<=L-1; j++){
			distances[i][j] = (double*) malloc(sizeof(double)*(2*L)) + L;
			for(int k=-L; k<=L-1; k++){
				distances[i][j][k]=-1;
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
							for(int k=0; k<3; k++) dist += xyz[k]*xyz[k];
							minDist = MIN(dist, minDist);
						}
					}
				}
				
				for(int dt=0; dt<=L; dt += L){
					for(int du=0; du<=L; du += L){
						for(int dv=0; dv<=L; dv += L){
							distances[t-dt][u-du][v-dv] = minDist/2.0;
// 							printf("%i %i %i: %lf\n", t-dt, u-du, v-dv, minDist/2.0);
						}
					}
				}
			}
		}
	}
	
	for(int t=-L; t<L; t++){
		for(int u=-L; u<L; u++){
			for(int v=-L; v<L; v++){
				if(distances[t][u][v] == -1){
					printf("Error: distance not filled in\n");
					exit(192);
				}
			}
		}
	}
	
// 	exit(0);
	return distances;
}



