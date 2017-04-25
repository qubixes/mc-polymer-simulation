#include "denspol_mc.h"

int CheckMutation(int mutation, int coor, int* topoState, LookupTables* lt){
	int topo = topoState[coor];
	if(lt->mutTopo[topo][mutation]<0) return 1;
	else return 0;
}

void PerformMutation(int mutation, int coor, int* topoState, LookupTables* lt){
	int topo = topoState[coor];
	topoState[coor] = lt->mutTopo[topo][mutation];
// 	if(lt->mutTopo[topo][mutation] < 0){ printf("????\n"); exit(192); }
}

int TransStep(CurState* cs, LookupTables* lt){
	int mono = DRng(cs->rngState)*cs->nPol*cs->polSize;
	int iMono = mono%cs->polSize;
	int iPol = mono/cs->polSize;
	int L = cs->L;
	int* mutations;
	
// 	printf("iPol=%i, iMono=%i\n", iPol, iMono);
	int unit1 = cs->unitPol[iPol][iMono];
	int coor[3] = {cs->coorPol[iPol][iMono],-1,-1};
	int unit2 = cs->unitPol[iPol][(iMono+1)%cs->polSize];
// 	printf("Faulted?\n");
	int rand = 4*DRng(cs->rngState);

	mutations = lt->mutator[unit1][unit2][rand];
	if(mutations[0]<0) return 0;
	int* newUnits = lt->newUnits[unit1][unit2][rand];
	
	///Forward move
	if(!(unit1 && unit2)){
		coor[1] = AddUnitToCoor(newUnits[0], coor[0], L);
		coor[2] = cs->coorPol[iPol][(iMono+2)%cs->polSize];
		
// 		if(coor[1] == coor[2]){
// 			printf("Wtf?\n");
// 			exit(192);
// 		}
		
// 		printf("new: %x [%x %x] %x, coor=(%i %i %i)\n", cs->unitPol[iPol][(iMono-1+cs->polSize)%cs->polSize], newUnits[0], newUnits[1], cs->unitPol[iPol][(iMono+2)%cs->polSize], coor[0], coor[1], coor[2]);
// 		printf("old: [%x %x] \n", unit1, unit2);
		
		for(int k=0; k<3; k++){
			lt->counts[unit1][unit2][rand][k]++;
// 			printf("CheckMutation: %i , coor=%i, topo=%i\n", mutations[k], coor[k], cs->topoState[coor[k]]);
			if(CheckMutation(mutations[k], coor[k], cs->topoState, lt)){
// 				exit(0);
				return 0;
			}	
// 			printf("Succes\n");
		}
// 		printf("topo[%i] = %i\n", coor[0], cs->topoState[coor[0]]);
// 		printf("\n");
		for(int k=0; k<3; k++){
			PerformMutation(mutations[k], coor[k], cs->topoState, lt);
		}
// 		printf("Newtopo[%i] = %i\n", coor[0], cs->topoState[coor[0]]);

// 		///Check back mutations:
// 		
// 		int slot = (unit1)?1:0;
// 		int* bMutations = lt->mutator[newUnits[0]][newUnits[1]][slot];
// 		for(int k=0; k<3; k++){
// 			if(CheckMutation(bMutations[k], coor[k], cs->topoState, lt)){
// 				printf("Error: no back mutation, failed at %i\n", k);
// 				printf("%x + %x => %x + %x\n", unit1, unit2, newUnits[0], newUnits[1]);
// 				for(int i=0; i<3; i++){
// 					printf("%i vs %i\n", mutations[i], bMutations[i]);
// 				}
// 				exit(0);
// 			}
// 		}
// 			
		
// 		lt->counts[unit1][unit2][rand][3]++;
		
		cs->unitPol[iPol][iMono] = newUnits[0];
		cs->unitPol[iPol][(iMono+1)%cs->polSize] = newUnits[1];
		cs->coorPol[iPol][(iMono+1)%cs->polSize] = coor[1];
		return 1;
// 		CheckIntegrity(cs, "Trans1");
	}
	///Backward move
	else {
		coor[1] = cs->coorPol[iPol][(iMono+1)%cs->polSize];
		coor[2] = cs->coorPol[iPol][(iMono+2)%cs->polSize];
		
		for(int k=0; k<3; k++){
			lt->counts[unit1][unit2][rand][k]++;
			if(CheckMutation(mutations[k], coor[k], cs->topoState, lt))
				return 0;
		}
		
		for(int k=0; k<3; k++){
			PerformMutation(mutations[k], coor[k], cs->topoState, lt);
		}
		
// 		///Check forward mutations:
// 		
// 		int slot;
// 		for(slot=0; slot<4; slot++){
// 			if(lt->newUnits[newUnits[0]][newUnits[1]][slot][0] == unit1) break;
// 		}
// 		if(slot==4){ printf("Wtf?\n"); exit(0); }
// 		int* bMutations = lt->mutator[newUnits[0]][newUnits[1]][slot];
// 		for(int k=0; k<3; k++){
// 			if(CheckMutation(bMutations[k], coor[k], cs->topoState, lt)){
// 				printf("Error: no forward mutation, failed at %i\n", k);
// 				printf("%x + %x => %x + %x\n", unit1, unit2, newUnits[0], newUnits[1]);
// 				for(int i=0; i<3; i++){
// 					printf("%i vs %i\n", mutations[i], bMutations[i]);
// 				}
// 				exit(0);
// 			}
// 		}
		
		
		
// 		lt->counts[unit1][unit2][rand][3]++;
		
		cs->unitPol[iPol][iMono] = newUnits[0];
		cs->unitPol[iPol][(iMono+1)%cs->polSize] = newUnits[1];
		cs->coorPol[iPol][(iMono+1)%cs->polSize] = (newUnits[0])?coor[2]:coor[0];
// 		printf("new: [%x %x], coor=(%i %i %i)\n", newUnits[0], newUnits[1], 
// 		CheckIntegrity(cs, "Trans2");
		return 1;
	}
	return 1;
}

int DiffuseStep(CurState* cs){
	int mono = DRng(cs->rngState)*cs->nPol*cs->polSize;
	int iMono = mono%cs->polSize;
	int iPol = mono/cs->polSize;
	
	int iPrev = (iMono-1+cs->polSize)%cs->polSize;
	
// 	printf("iPrev=%i, iMono=%i, iPol=%i\n", iPrev, iMono, iPol);
	
	if(cs->unitPol[iPol][iMono] && !cs->unitPol[iPol][iPrev]){
		cs->unitPol[iPol][iPrev] = cs->unitPol[iPol][iMono];
		cs->unitPol[iPol][iMono] = 0;
		cs->coorPol[iPol][iMono] = cs->coorPol[iPol][(iMono+1)%cs->polSize];
// 		CheckIntegrity(cs, "Dif1");
// 		exit(0);
		return 1;
	}
	else if(!cs->unitPol[iPol][iMono] && cs->unitPol[iPol][iPrev]){
		cs->unitPol[iPol][iMono] = cs->unitPol[iPol][iPrev];
		cs->unitPol[iPol][iPrev] = 0;
		cs->coorPol[iPol][iMono] = cs->coorPol[iPol][iPrev];
// 		CheckIntegrity(cs, "Dif2");
		return 1;
	}
	return 0;
}

double MeasSl(CurState* cs){
	long nSl=0;
	for(int iPol=0; iPol<cs->nPol; iPol++){
		for(int iBond=0; iBond<cs->polSize; iBond++){
			if(!cs->unitPol[iPol][iBond]) nSl++;
		}
	}
	return nSl/(double)(cs->nPol*cs->polSize);
}

double DoMCStep(long nStep, CurState* cs, LookupTables* lt){
// 	long nAcc=0;
	long nAccTrans=0;
	long nAccDiff=0;
// 	long interval=10000;
// 	long nMeas=0;
// 	double slRat=0;
	
	for(int iStep=0; iStep<cs->polSize*cs->nPol*nStep; iStep++){
		nAccTrans += TransStep(cs, lt);
		nAccDiff += DiffuseStep(cs);
// 		if(iStep%interval == interval-1){
// 			slRat += MeasSl(cs);
// 			nMeas++;
// 		}
	}
	
// 	slRat /= nMeas;
	double ratTrans = nAccTrans/(double)(cs->polSize*cs->nPol*nStep);
// 	double ratDiff = nAccDiff/(double)(cs->polSize*cs->nPol*nStep);
// 	printf("rat = {%lf %lf}\n", ratTrans, ratDiff);
// 	printf("sl= %lf\n", slRat);
// 	printf("True density = %lf\n", cs->density*(1-slRat));
// 	CheckIntegrity(cs, "Check after MCSteps\n");
	return ratTrans;
}



