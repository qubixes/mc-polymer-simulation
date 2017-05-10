#include "denspol_mc.h"

int CheckMutation(int mutation, int coor, int* topoState, LookupTables* lt){
	int topo = topoState[coor];
	if(lt->mutTopo[topo][mutation]<0) return 1;
	else return 0;
}

void PerformMutation(int mutation, int coor, int* topoState, LookupTables* lt){
	int topo = topoState[coor];
	topoState[coor] = lt->mutTopo[topo][mutation];
}

int TransStep(CurState* cs, LookupTables* lt){
	int mono = DRng(cs->rngState)*cs->nPol*cs->polSize;
	int iMono = mono%cs->polSize;
	int iPol = mono/cs->polSize;
	int L = cs->L;
	int* mutations;
	
	int unit1 = cs->unitPol[iPol][iMono];
	int coor[3] = {cs->coorPol[iPol][iMono],-1,-1};
	int unit2 = cs->unitPol[iPol][(iMono+1)%cs->polSize];
	int rand = 4*DRng(cs->rngState);
	
	mutations = lt->mutator[unit1][unit2][rand];
	if(mutations[0]<0) return 0;
	int* newUnits = lt->newUnits[unit1][unit2][rand];
	
	///Figure out the bonds outside the two we have.
	int jMono=(iMono-1+cs->polSize)%cs->polSize;
	while(!cs->unitPol[iPol][jMono]) jMono = (jMono-1+cs->polSize)%cs->polSize;
	int unit0=cs->unitPol[iPol][jMono];
	jMono=(iMono+2)%cs->polSize;
	while(!cs->unitPol[iPol][jMono]) jMono = (jMono+1)%cs->polSize;
	int unit3=cs->unitPol[iPol][jMono];
	
	
	///Forward move
	if(!(unit1 && unit2)){
		
		int dBend=0;
		dBend -= lt->unitDot[unit0][unit1|unit2];
		dBend -= lt->unitDot[unit1|unit2][unit3];
		dBend += lt->unitDot[unit0][newUnits[0]];
		dBend += lt->unitDot[newUnits[0]][newUnits[1]];
		dBend += lt->unitDot[newUnits[1]][unit3];
		
		if(lt->bendProb[dBend+BEND_LVL] < 1 && DRng(cs->rngState) > lt->bendProb[dBend+BEND_LVL]) 
			return 0;
		
		
		coor[1] = AddUnitToCoor(newUnits[0], coor[0], L);
		coor[2] = cs->coorPol[iPol][(iMono+2)%cs->polSize];
		
		for(int k=0; k<3; k++){
			lt->counts[unit1][unit2][rand][k]++;
			if(CheckMutation(mutations[k], coor[k], cs->topoState, lt)){
				return 0;
			}
		}
		for(int k=0; k<3; k++){
			PerformMutation(mutations[k], coor[k], cs->topoState, lt);
		}
		
		cs->unitPol[iPol][iMono] = newUnits[0];
		cs->unitPol[iPol][(iMono+1)%cs->polSize] = newUnits[1];
		cs->coorPol[iPol][(iMono+1)%cs->polSize] = coor[1];
		return 1;
	}
	///Backward move
	else {
		
		int dBend=0;
		dBend -= lt->unitDot[unit0][unit1];
		dBend -= lt->unitDot[unit1][unit2];
		dBend -= lt->unitDot[unit2][unit3];
		dBend += lt->unitDot[unit0][newUnits[0]|newUnits[1]];
		dBend += lt->unitDot[newUnits[0]|newUnits[1]][unit3];
		
		if(lt->bendProb[dBend+BEND_LVL] < 1 && DRng(cs->rngState) > lt->bendProb[dBend+BEND_LVL]) 
			return 0;
		
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
		
		cs->unitPol[iPol][iMono] = newUnits[0];
		cs->unitPol[iPol][(iMono+1)%cs->polSize] = newUnits[1];
		cs->coorPol[iPol][(iMono+1)%cs->polSize] = (newUnits[0])?coor[2]:coor[0];
		return 1;
	}
	return 1;
}

int DiffuseStep(CurState* cs){
	int mono = DRng(cs->rngState)*cs->nPol*cs->polSize;
	int iMono = mono%cs->polSize;
	int iPol = mono/cs->polSize;
	
	int iPrev = (iMono-1+cs->polSize)%cs->polSize;
	
	if(cs->unitPol[iPol][iMono] && !cs->unitPol[iPol][iPrev]){
		cs->unitPol[iPol][iPrev] = cs->unitPol[iPol][iMono];
		cs->unitPol[iPol][iMono] = 0;
		cs->coorPol[iPol][iMono] = cs->coorPol[iPol][(iMono+1)%cs->polSize];
		return 1;
	}
	else if(!cs->unitPol[iPol][iMono] && cs->unitPol[iPol][iPrev]){
		cs->unitPol[iPol][iMono] = cs->unitPol[iPol][iPrev];
		cs->unitPol[iPol][iPrev] = 0;
		cs->coorPol[iPol][iMono] = cs->coorPol[iPol][iPrev];
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

void MeasBends(CurState* cs, LookupTables* lt, long counts[4]){
	
	for(int iPol=0; iPol<cs->nPol; iPol++){
		
		for(int iBond=0; iBond<cs->polSize; iBond++){
			if(cs->unitPol[iPol][iBond] == 0) continue;
			int jBond=(iBond+1)%cs->polSize; 
			while(cs->unitPol[iPol][jBond] == 0) jBond= (jBond+1)%cs->polSize;
			
			int unitI = cs->unitPol[iPol][iBond];
			int unitJ = cs->unitPol[iPol][jBond];
			int bend = lt->unitDot[unitI][unitJ]+1;
			if(bend<0 || bend>=4){ printf("Error in bend finding!\n"); exit(0);}
			counts[bend]++;
		}
	}
}

double DoMCStep(long nStep, CurState* cs, LookupTables* lt){
	long nAccTrans=0;
	long nAccDiff=0;
// 	long bendCounts[4]={0,0,0,0};
// 	double sl=0;
// 	long nMeas=0;
	
	for(int iStep=0; iStep<cs->polSize*cs->nPol*nStep; iStep++){
		nAccTrans += TransStep(cs, lt);
		nAccDiff += DiffuseStep(cs);
// 		if(iStep%(cs->polSize*cs->nPol) == cs->polSize*cs->nPol-1){
// 			MeasBends(cs, lt, bendCounts);
// 			sl += MeasSl(cs);
// 			nMeas++;
// 		}
	}
	
// 	long totBends=0;
// 	for(int i=0; i<4; i++)
// 		totBends += bendCounts[i];
// 	printf("\n");
// 	double avgBend=0;
// 	for(int i=0; i<4; i++){
// 		avgBend += (i-1)*bendCounts[i]/(double)totBends;
// 		printf("%i %.2lf\n", i-1, bendCounts[i]/(double)totBends);
// 	}
// 	printf("Avg= %.3lf\n", avgBend);
// 	printf("Sl = %.2lf\n", sl/nMeas);
	double ratTrans = nAccTrans/(double)(cs->polSize*cs->nPol*nStep);
	printf("Trans moves accepted: %.2lf %%\n", 100*ratTrans);
	return ratTrans;
}



