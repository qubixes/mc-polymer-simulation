#include "denspol_mc.h"

/// The "noTopo" flag is raised if the backbone of the polymer is length one 
/// which in this case means the 1st unit vector is the only non-zero unit vector
/// Obviously this means that nothing topologically is happening. 

int StartMove(CurState* cs, LookupTables* lt, Polymer* pol){
	if(pol->unitPol[0] == 0x0) return 0;
	
	int unit0=pol->unitPol[0];
	int coor0=pol->coorPol[0];
	int coor1=pol->coorPol[1];
	
	int newUnit = Rng(cs->rngState)&0xf;
	if(!IsValid(newUnit)) return 0;
	
	int newCoor, OOB;
	newCoor = cs->AddUnitToCoor(newUnit^0xf, coor1, cs->L, &OOB);
	if(OOB) return 0;
	
	int jMono = 1;
	while(jMono<pol->polSize && !pol->unitPol[jMono]) jMono++;
	int noTopo= (jMono==pol->polSize)?1:0;
	int unit1 = (jMono==pol->polSize)?0xf:pol->unitPol[jMono];
	
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
		if(!TestMoveHP(cs,lt,monoMove, pol, unitMove)) 
			return 0;
#endif
		if(mut != SAME_TOPO)
			cs->topoState[coor1] = lt->topComp[cs->topoState[coor1]].mutators[mut];
	}
	
	cs->bondOcc[coor0]   ^= 1<<(unit0^0xf);
	cs->bondOcc[coor1]   ^= 1<<unit0;
	cs->bondOcc[newCoor] |= 1<<(newUnit^0xf);
	cs->bondOcc[coor1]   |= 1<<newUnit;
	
	pol->unitPol[0] = newUnit;
	pol->coorPol[0] = newCoor;
// 	CheckIntegrity(cs, "After start move");
	return 1;
}

int EndMove(CurState* cs, LookupTables* lt, Polymer* pol){
	int mono0=pol->polSize;
	int mono1=pol->polSize-1;
		
	if(pol->unitPol[mono1] == 0x0) return 0;
	
	int unit0=pol->unitPol[mono1];
	int coor0=pol->coorPol[mono0];
	int coor1=pol->coorPol[mono1];
	
	int newUnit = Rng(cs->rngState)&0xf;
	if(!IsValid(newUnit)) return 0;
	
	int newCoor, OOB;
	newCoor = cs->AddUnitToCoor(newUnit, coor1, cs->L, &OOB);
	if(OOB) return 0;
	
	int jMono = mono1-1;
	while(jMono>=0 && !pol->unitPol[jMono]) jMono--;
	int noTopo=(jMono<0)?1:0;
	int unit1 = (jMono<0)?0xf:pol->unitPol[jMono];
	
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
		if(!TestMoveHP(cs,lt,monoMove, pol, unitMove)) 
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
	
	pol->unitPol[mono1] = newUnit;
	pol->coorPol[mono0] = newCoor;
// 	CheckIntegrity(cs, "After end move");
	return 1;
}

/** 0 <= iMono < polSize-1
  * 
  **/

int TransMoveLinear(CurState* cs, LookupTables* lt, Polymer* pol, int iMono){
	int L = cs->L;
	int* mutations;
	
	int mut[3];
	int unit1 = pol->unitPol[iMono];
	int coor[3] = {pol->coorPol[iMono],-1,-1};
	int unit2 = pol->unitPol[iMono+1];
	int rand = 4*DRng(cs->rngState);
	
	mutations = lt->mutator[unit1][unit2][rand];
	if(mutations[0]<0) return 0;
	int* newUnits = lt->newUnits[unit1][unit2][rand];
	
	///Figure out the bonds outside the two we have.
	int jMono=iMono-1;
	while(jMono >= 0 && !pol->unitPol[jMono]) jMono--;
	int unit0ringed=(jMono<0)?1:0;
	int unit0=unit0ringed?0xf:pol->unitPol[jMono];
	jMono=iMono+2;
	while(jMono<pol->polSize && !pol->unitPol[jMono]) jMono++;
	int unit3ringed=(jMono>=pol->polSize)?1:0;
	int unit3=unit3ringed?0xf:pol->unitPol[jMono];
	
	
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
		
		int OOB;
		coor[1] = cs->AddUnitToCoor(newUnits[0], coor[0], L, &OOB);
		if(OOB) return 0;
		
		coor[2] = pol->coorPol[iMono+2];
		
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
		int monoMove = (iMono+1)%pol->polSize;
		int unitMove = unit2?newUnits[0]:(newUnits[1]^0xf);
		if(!TestMoveHP(cs,lt,monoMove, pol, unitMove)) 
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
		
		pol->unitPol[iMono]   = newUnits[0];
		pol->unitPol[iMono+1] = newUnits[1];
		pol->coorPol[iMono+1] = coor[1];
// 		CheckIntegrity(cs, "After Forward Linear Move");
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
		
		
		
		coor[1] = pol->coorPol[iMono+1];
		coor[2] = pol->coorPol[iMono+2];
		
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
		if(!TestMoveHP(cs,lt,monoMove, pol, unitMove))
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
		
		pol->unitPol[iMono  ] = newUnits[0];
		pol->unitPol[iMono+1] = newUnits[1];
		pol->coorPol[iMono+1] = (newUnits[0])?coor[2]:coor[0];
// 		CheckIntegrity(cs, "After Backward Linear Move");
 		return 1;
	}
	return 1;
}




