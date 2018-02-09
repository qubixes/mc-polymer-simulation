#include "denspol_mc.h"

/// The "noTopo" flag is raised if the backbone of the polymer is length one 
/// which in this case means the 1st unit vector is the only non-zero unit vector
/// Obviously this means that nothing topologically is happening. 

int StartMove(CurState* cs, LookupTables* lt, int iPol){
	if(cs->unitPol[iPol][0] == 0x0) return 0;
	
	int unit0=cs->unitPol[iPol][0];
	int coor0=cs->coorPol[iPol][0];
	int coor1=cs->coorPol[iPol][1];
	
	int newUnit = Rng(cs->rngState)&0xf;
	if(!IsValid(newUnit)) return 0;
	
	int newCoor = AddUnitToCoor(newUnit^0xf, coor1, cs->L);
	int jMono = 1;
	while(jMono<cs->polSize && !cs->unitPol[iPol][jMono]) jMono++;
	int noTopo= (jMono==cs->polSize)?1:0;
	int unit1 = (jMono==cs->polSize)?0xf:cs->unitPol[iPol][jMono];
	
	if(cs->bondOcc[newCoor]&(1<<(newUnit^0xf))) return 0;
	
	if(!noTopo){
		int dBend=0;
		dBend -= lt->unitDot[unit0][unit1];
		dBend += lt->unitDot[newUnit][unit1];
		if(lt->bendProb[dBend+BEND_LVL] < 1 && DRng(cs->rngState) > lt->bendProb[dBend+BEND_LVL])
			return 0;
		
		int mut=NON_EXISTING;
		if(lt->topComp[cs->topoState[coor1]].permBond&(1<<unit0))
			mut = lt->mutIdTableDouble[unit0][newUnit^0xf];
		else
			mut = lt->mutIdTableTriple[newUnit][unit0^0xf][unit1];
		if(mut == NON_EXISTING) return 0;
		if(mut >= 0 && lt->topComp[cs->topoState[coor1]].mutators[mut] == NON_EXISTING)
			return 0;
		
#ifdef __HP_ENABLED__
		int monoMove = 0;
		int unitMove = AddUnitVectors(unit0, newUnit^0xf);
		if(!TestMoveHP(cs,lt,monoMove, iPol, unitMove)) 
			return 0;
#endif
		if(mut != SAME_TOPO)
			cs->topoState[coor1] = lt->topComp[cs->topoState[coor1]].mutators[mut];
	}
	
	cs->bondOcc[coor0]   ^= 1<<(unit0^0xf);
	cs->bondOcc[coor1]   ^= 1<<unit0;
	cs->bondOcc[newCoor] |= 1<<(newUnit^0xf);
	cs->bondOcc[coor1]   |= 1<<newUnit;
	
	cs->unitPol[iPol][0] = newUnit;
	cs->coorPol[iPol][0] = newCoor;
// 	CheckIntegrity(cs, "After start move");
	return 1;
}

int EndMove(CurState* cs, LookupTables* lt, int iPol){
	int mono0=cs->polSize;
	int mono1=cs->polSize-1;
		
	if(cs->unitPol[iPol][mono1] == 0x0) return 0;
	
	int unit0=cs->unitPol[iPol][mono1];
	int coor0=cs->coorPol[iPol][mono0];
	int coor1=cs->coorPol[iPol][mono1];
	
	int newUnit = Rng(cs->rngState)&0xf;
	if(!IsValid(newUnit)) return 0;
	
	int newCoor = AddUnitToCoor(newUnit, coor1, cs->L);
	int jMono = mono1-1;
	while(jMono>=0 && !cs->unitPol[iPol][jMono]) jMono--;
	int noTopo=(jMono<0)?1:0;
	int unit1 = (jMono<0)?0xf:cs->unitPol[iPol][jMono];
	
	if(cs->bondOcc[newCoor]&(1<<newUnit)) return 0;

	
	if(!noTopo){
		int dBend=0;
		dBend -= lt->unitDot[unit1][unit0];
		dBend += lt->unitDot[unit1][newUnit];
		if(lt->bendProb[dBend+BEND_LVL] < 1 && DRng(cs->rngState) > lt->bendProb[dBend+BEND_LVL])
			return 0;
		
		int mut = NON_EXISTING;
		if(lt->topComp[cs->topoState[coor1]].permBond&(1<<(unit0^0xf)))
			mut = lt->mutIdTableDouble[unit0^0xf][newUnit];
		else
			mut = lt->mutIdTableTriple[unit1][unit0][newUnit];
		if(mut == NON_EXISTING) return 0;
		if(mut >= 0 && lt->topComp[cs->topoState[coor1]].mutators[mut] == NON_EXISTING)
			return 0;
		
#ifdef __HP_ENABLED__
		int monoMove = mono0;
		int unitMove = AddUnitVectors(unit0^0xf, newUnit);
		if(!TestMoveHP(cs,lt,monoMove, iPol, unitMove)) 
			return 0;
#endif
		
		if(mut != SAME_TOPO)
			cs->topoState[coor1] = lt->topComp[cs->topoState[coor1]].mutators[mut];
	}
	
// 	cs->bondOcc[coor0]   ^= 1<<(unit0);
	
	cs->bondOcc[coor1]   ^= 1<<(unit0^0xf);
	cs->bondOcc[coor1]   |= 1<<(newUnit^0xf);
	cs->bondOcc[coor0]   ^= 1<<unit0;
	cs->bondOcc[newCoor] |= 1<<newUnit;
	
	cs->unitPol[iPol][mono1] = newUnit;
	cs->coorPol[iPol][mono0] = newCoor;
// 	CheckIntegrity(cs, "After end move");
	return 1;
}

/** 0 <= iMono < polSize-1
  * 
  **/

int TransStepCompactLinear(CurState* cs, LookupTables* lt, int iMono, int iPol){
// 	int mono = DRng(cs->rngState)*cs->nPol*cs->polSize;
// 	int iMono = mono%cs->polSize;
// 	int iPol = mono/cs->polSize;
	int L = cs->L;
	int* mutations;
	
	int mut[3];
	int unit1 = cs->unitPol[iPol][iMono];
	int coor[3] = {cs->coorPol[iPol][iMono],-1,-1};
	int unit2 = cs->unitPol[iPol][iMono+1];
	int rand = 4*DRng(cs->rngState);
	
	mutations = lt->mutator[unit1][unit2][rand];
	if(mutations[0]<0) return 0;
	int* newUnits = lt->newUnits[unit1][unit2][rand];
	
	///Figure out the bonds outside the two we have.
	int jMono=iMono-1;
	while(jMono >= 0 && !cs->unitPol[iPol][jMono]) jMono--;
	int unit0ringed=(jMono<0)?1:0;
	int unit0=unit0ringed?0xf:cs->unitPol[iPol][jMono];
	jMono=iMono+2;
	while(jMono<cs->polSize && !cs->unitPol[iPol][jMono]) jMono++;
	int unit3ringed=(jMono>=cs->polSize)?1:0;
	int unit3=unit3ringed?0xf:cs->unitPol[iPol][jMono];
	
	
	///Forward move
	if(!(unit1 && unit2)){
		
		int dBend=0;
		if(!unit0ringed){
			dBend -= lt->unitDot[unit0][unit1|unit2];
			dBend += lt->unitDot[unit0][newUnits[0]];
		}
		if(!unit3ringed){
			dBend -= lt->unitDot[unit1|unit2][unit3];
			dBend += lt->unitDot[newUnits[1]][unit3];
		}
		dBend += lt->unitDot[newUnits[0]][newUnits[1]];
		
		if(lt->bendProb[dBend+BEND_LVL] < 1 && DRng(cs->rngState) > lt->bendProb[dBend+BEND_LVL])
			return 0;
		
		
		coor[1] = AddUnitToCoor(newUnits[0], coor[0], L);
		coor[2] = cs->coorPol[iPol][iMono+2];
		
		///Check if new position is empty
		if(cs->bondOcc[coor[1]]&((1<<newUnits[0])|(1<<(newUnits[1]^0xf)))) return 0;
		
		///Check if the left position is able to move
		if(!unit0ringed){
			if(lt->topComp[cs->topoState[coor[0]]].permBond&(1<<unit0))
				mut[0] = lt->mutIdTableDouble[(unit1|unit2)^0xf][newUnits[0]];
			else
				mut[0] = lt->mutIdTableTriple[unit0][unit1|unit2][newUnits[0]];
			if(mut[0] == NON_EXISTING) return 0;
			if(mut[0] >= 0 && lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]] == NON_EXISTING)
				return 0;
		}
		
		///Check if right position is able to move
		if(!unit3ringed){
			if(lt->topComp[cs->topoState[coor[2]]].permBond&(1<<(unit3^0xf)))
				mut[1] = lt->mutIdTableDouble[(unit1|unit2)][newUnits[1]^0xf];
			else
				mut[1] = lt->mutIdTableTriple[newUnits[1]][(unit1|unit2)^0xf][unit3];
			if(mut[1] == NON_EXISTING) return 0;
			if(mut[1] >= 0 && lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]] == NON_EXISTING)
				return 0;
		}
		
#ifdef __HP_ENABLED__
		int monoMove = (iMono+1)%cs->polSize;
		int unitMove = unit2?newUnits[0]:(newUnits[1]^0xf);
		if(!TestMoveHP(cs,lt,monoMove, iPol, unitMove)) 
			return 0;
#endif
		
		///Change the topological state of the left position
		if(!unit0ringed && mut[0] != SAME_TOPO)
			cs->topoState[coor[0]] = lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]];
		
		///Change the topological state of the right position
		if(!unit3ringed && mut[1] != SAME_TOPO)
			cs->topoState[coor[2]] = lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]];
		
		
		cs->bondOcc[coor[0]] ^= 1<<((unit1|unit2)^0xf);
		cs->bondOcc[coor[0]] |= 1<<(newUnits[0]^0xf);
		cs->bondOcc[coor[1]] |= 1<< newUnits[0];
		cs->bondOcc[coor[1]] |= 1<<(newUnits[1]^0xf);
		cs->bondOcc[coor[2]] ^= 1<<(unit1|unit2);
		cs->bondOcc[coor[2]] |= 1<< newUnits[1];
		
		cs->unitPol[iPol][iMono] = newUnits[0];
		cs->unitPol[iPol][iMono+1] = newUnits[1];
		cs->coorPol[iPol][iMono+1] = coor[1];
// 		CheckIntegrity(cs, "After Forward Move");
// 			printf("iMono = %i, iPol = %i, coor = (%i,%i,%i)\n", iMono+1, iPol, coor[0], coor[1], coor[2]);
// 			exit(192);
// 		}
		return 1;
	}
	///Backward move
	else {
		
		int dBend=0;
		if(!unit0ringed){
			dBend -= lt->unitDot[unit0][unit1];
			dBend += lt->unitDot[unit0][newUnits[0]|newUnits[1]];
		}
		if(!unit3ringed){
			dBend -= lt->unitDot[unit2][unit3];
			dBend += lt->unitDot[newUnits[0]|newUnits[1]][unit3];
		}
		dBend -= lt->unitDot[unit1][unit2];

		if(lt->bendProb[dBend+BEND_LVL] < 1 && DRng(cs->rngState) > lt->bendProb[dBend+BEND_LVL]) 
			return 0;
		
		
		
		coor[1] = cs->coorPol[iPol][iMono+1];
		coor[2] = cs->coorPol[iPol][iMono+2];
		
		///Check if we can remove the middle monomer
		if(lt->topComp[cs->topoState[coor[1]]].permBond&(1<<unit1)) return 0;
		
		///Check if the new bond is already occupied
		if(cs->bondOcc[coor[2]]&(1<<(newUnits[0]|newUnits[1]))) return 0;
		
		///Check if the left bond can be modified
		if(!unit0ringed){
			if(lt->topComp[cs->topoState[coor[0]]].permBond&(1<<(unit1^0xf)))
				mut[0] = lt->mutIdTableDouble[unit1^0xf][newUnits[0]|newUnits[1]];
			else
				mut[0] = lt->mutIdTableTriple[unit0][unit1][newUnits[0]|newUnits[1]];
			if(mut[0] == NON_EXISTING) return 0;
			if(mut[0] != SAME_TOPO && lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]] == NON_EXISTING) 
				return 0;
		}
		
		///Check if the right bond can be modified
		if(!unit3ringed){
			if(lt->topComp[cs->topoState[coor[2]]].permBond&(1<<unit2))
				mut[1] = lt->mutIdTableDouble[unit2][(newUnits[0]|newUnits[1])^0xf];
			else
				mut[1] = lt->mutIdTableTriple[newUnits[0]|newUnits[1]][unit2^0xf][unit3];
			if(mut[1] == NON_EXISTING) return 0;
			if(mut[1] != SAME_TOPO && lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]] == NON_EXISTING)
				return 0;
		}
		
#ifdef __HP_ENABLED__
		int monoMove = iMono+1;
		int unitMove = newUnits[0]?unit2:(unit1^0xf);
		if(!TestMoveHP(cs,lt,monoMove, iPol, unitMove))
			return 0;
#endif
		
		///Change topological state of left side
		if(!unit0ringed && mut[0] != SAME_TOPO)
			cs->topoState[coor[0]] = lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]];
		
		///Change topological state of right side
		if(!unit3ringed && mut[1] != SAME_TOPO)
			cs->topoState[coor[2]] = lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]];
		
		cs->bondOcc[coor[0]] ^= 1<<(unit1^0xf);
		cs->bondOcc[coor[0]] |= 1<<((newUnits[0]|newUnits[1])^0xf);
		cs->bondOcc[coor[1]] ^= 1<<unit1;
		cs->bondOcc[coor[1]] ^= 1<<(unit2^0xf);
		cs->bondOcc[coor[2]] |= 1<<(newUnits[0]|newUnits[1]);
		cs->bondOcc[coor[2]] ^= 1<<unit2;
		
		cs->unitPol[iPol][iMono  ] = newUnits[0];
		cs->unitPol[iPol][iMono+1] = newUnits[1];
		cs->coorPol[iPol][iMono+1] = (newUnits[0])?coor[2]:coor[0];
// 		CheckIntegrity(cs, "After backward move");
// 			printf("Coor = (%i %i %i)\n", coor[0], coor[1], coor[2]);
// 			exit(192);
// 		};
		return 1;
	}
	return 1;
}

int DiffuseStepLinear(CurState* cs, LookupTables* lt){
	int mono = DRng(cs->rngState)*cs->nPol*(cs->polSize-1);
	int iMono = (mono%(cs->polSize-1))+1;
	int iPol = mono/(cs->polSize-1);
	
	int iPrev = iMono-1;
	
	if(cs->unitPol[iPol][iMono] && !cs->unitPol[iPol][iPrev]){
#ifdef __HP_ENABLED__
		int unitMove = cs->unitPol[iPol][iMono];
		if(!TestMoveHP(cs,lt,iMono, iPol, unitMove))
			return 0;
#endif
		cs->unitPol[iPol][iPrev] = cs->unitPol[iPol][iMono];
		cs->unitPol[iPol][iMono] = 0;
		cs->coorPol[iPol][iMono] = cs->coorPol[iPol][iMono+1];
		return 1;
	}
	else if(!cs->unitPol[iPol][iMono] && cs->unitPol[iPol][iPrev]){
#ifdef __HP_ENABLED__
		int unitMove = cs->unitPol[iPol][iPrev]^0xf;
		if(!TestMoveHP(cs,lt,iMono, iPol, unitMove)) 
			return 0;
#endif
		cs->unitPol[iPol][iMono] = cs->unitPol[iPol][iPrev];
		cs->unitPol[iPol][iPrev] = 0;
		cs->coorPol[iPol][iMono] = cs->coorPol[iPol][iPrev];
		return 1;
	}
// 	CheckIntegrity(cs, "After diffusive move");
	return 0;
}

double DoMCStep(long nStep, CurState* cs, LookupTables* lt){
	long nAccTrans=0;
	long nAccDiff=0;
	
	for(long iStep=0; iStep<cs->polSize*cs->nPol*nStep; iStep++){
		int mono = DRng(cs->rngState)*cs->nPol*(cs->polSize+1);
		int iMono = mono%(cs->polSize+1);
		int iPol = mono/(cs->polSize+1);
		
		if(iMono==cs->polSize-1)
			StartMove(cs, lt, iPol);
		else if (iMono == cs->polSize)
			EndMove(cs, lt, iPol);
		else
			TransStepCompactLinear(cs, lt, iMono, iPol);
#ifdef __TOPO_MOVE_ENABLED__
		TopoMove(cs,lt);
#endif
		nAccDiff += DiffuseStepLinear(cs, lt);
// 		printf("Step %i/%i\n", iStep, cs->polSize*cs->nPol*nStep);
	}
// 	CheckIntegrity(cs, "After step");
	double ratTrans = nAccTrans/(double)(cs->polSize*cs->nPol*nStep);
	return ratTrans;
}


