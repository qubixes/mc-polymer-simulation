#include "denspol_mc.h"

/// Forward transverse move: unit1+unit2 => newUnits[0]+newUnits[1]

int TransMoveRing(CurState* cs, LookupTables* lt, Polymer* pol, int iMono){
	int L = cs->L;
	int* mutations;
	
	int mut[3];
	int unit1 = pol->unitPol[iMono];
	int coor[3] = {pol->coorPol[iMono],-1,-1};
	int unit2 = pol->unitPol[(iMono+1)%pol->polSize];
	int rand = 4*DRng(cs->rngState);
	
	mutations = lt->mutator[unit1][unit2][rand];
	if(mutations[0]<0) return 0;
	int* newUnits = lt->newUnits[unit1][unit2][rand];
	
	///Figure out the bonds outside the two we have.
	int jMono=(iMono-1+pol->polSize)%pol->polSize;
	while(!pol->unitPol[jMono]) jMono = (jMono-1+pol->polSize)%pol->polSize;
	int unit0=pol->unitPol[jMono];
	jMono=(iMono+2)%pol->polSize;
	while(!pol->unitPol[jMono]) jMono = (jMono+1)%pol->polSize;
	int unit3=pol->unitPol[jMono];
	
	
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
	
		int OOB;
		coor[1] = cs->AddUnitToCoor(newUnits[0], coor[0], L, &OOB);
		if(OOB) return 0;
		
		coor[2] = pol->coorPol[(iMono+2)%pol->polSize];
		
		///Check if new position is empty
		if(cs->bondOcc[coor[1]]&((1<<newUnits[0])|(1<<(newUnits[1]^0xf)))) return 0;
		
		///Check if the left position is able to move
		if(lt->topComp[cs->topoState[coor[0]]].permBond&(1<<unit0))
			mut[0] = lt->mutIdTableDouble[(unit1|unit2)^0xf][newUnits[0]];
		else
			mut[0] = lt->mutIdTableTriple[unit0][unit1|unit2][newUnits[0]];
		if(mut[0] == NON_EXISTING) return 0;
		if(mut[0] >= 0 && lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]] == NON_EXISTING)
			return 0;
		
		///Check if right position is able to move
		if(lt->topComp[cs->topoState[coor[2]]].permBond&(1<<(unit3^0xf)))
			mut[1] = lt->mutIdTableDouble[(unit1|unit2)][newUnits[1]^0xf];
		else
			mut[1] = lt->mutIdTableTriple[newUnits[1]][(unit1|unit2)^0xf][unit3];
		if(mut[1] == NON_EXISTING) return 0;
		if(mut[1] >= 0 && lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]] == NON_EXISTING)
			return 0;
		
#ifdef __HP_ENABLED__
		int monoMove = (iMono+1)%pol->polSize;
		int unitMove = unit2?newUnits[0]:(newUnits[1]^0xf);
		if(!TestMoveHP(cs,lt,monoMove, pol, unitMove)) 
			return 0;
// 		printf("success forward!\n");
#endif
		///Change the topological state of the left position
		if(mut[0] != SAME_TOPO)
			cs->topoState[coor[0]] = lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]];
		
		///Change the topological state of the right position
		if(mut[1] != SAME_TOPO)
			cs->topoState[coor[2]] = lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]];
		

		
		cs->bondOcc[coor[0]] ^= 1<<((unit1|unit2)^0xf);
		cs->bondOcc[coor[0]] |= 1<<(newUnits[0]^0xf);
		cs->bondOcc[coor[1]] |= 1<< newUnits[0];
		cs->bondOcc[coor[1]] |= 1<<(newUnits[1]^0xf);
		cs->bondOcc[coor[2]] ^= 1<<(unit1|unit2);
		cs->bondOcc[coor[2]] |= 1<< newUnits[1];
		
		pol->unitPol[iMono] = newUnits[0];
		pol->unitPol[(iMono+1)%pol->polSize] = newUnits[1];
		pol->coorPol[(iMono+1)%pol->polSize] = coor[1];
// 		CheckIntegrity(cs, "After Forward Linear Move");
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
		
		
		
		coor[1] = pol->coorPol[(iMono+1)%pol->polSize];
		coor[2] = pol->coorPol[(iMono+2)%pol->polSize];
		
		///Check if we can remove the middle monomer
		if(lt->topComp[cs->topoState[coor[1]]].permBond&(1<<unit1)) return 0;
		
		///Check if the new bond is already occupied
		if(cs->bondOcc[coor[2]]&(1<<(newUnits[0]|newUnits[1]))) return 0;
		
		///Check if the left bond can be modified
		if(lt->topComp[cs->topoState[coor[0]]].permBond&(1<<(unit1^0xf)))
			mut[0] = lt->mutIdTableDouble[unit1^0xf][newUnits[0]|newUnits[1]];
		else
			mut[0] = lt->mutIdTableTriple[unit0][unit1][newUnits[0]|newUnits[1]];
		if(mut[0] == NON_EXISTING) return 0;
		if(mut[0] != SAME_TOPO && lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]] == NON_EXISTING) 
			return 0;
		
		
		///Check if the right bond can be modified
		if(lt->topComp[cs->topoState[coor[2]]].permBond&(1<<unit2))
			mut[1] = lt->mutIdTableDouble[unit2][(newUnits[0]|newUnits[1])^0xf];
		else
			mut[1] = lt->mutIdTableTriple[newUnits[0]|newUnits[1]][unit2^0xf][unit3];
		if(mut[1] == NON_EXISTING) return 0;
		if(mut[1] != SAME_TOPO && lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]] == NON_EXISTING)
			return 0;
		
#ifdef __HP_ENABLED__
		int monoMove = ((iMono+1)%pol->polSize);
		int unitMove = newUnits[0]?unit2:(unit1^0xf);
		if(!TestMoveHP(cs,lt,monoMove, pol, unitMove))
			return 0;
// 		printf("success backward!\n");
#endif
		
		///Change topological state of left side
		if(mut[0] != SAME_TOPO)
			cs->topoState[coor[0]] = lt->topComp[cs->topoState[coor[0]]].mutators[mut[0]];
		
		///Change topological state of right side
		if(mut[1] != SAME_TOPO)
			cs->topoState[coor[2]] = lt->topComp[cs->topoState[coor[2]]].mutators[mut[1]];
		
		cs->bondOcc[coor[0]] ^= 1<<(unit1^0xf);
		cs->bondOcc[coor[0]] |= 1<<((newUnits[0]|newUnits[1])^0xf);
		cs->bondOcc[coor[1]] ^= 1<<unit1;
		cs->bondOcc[coor[1]] ^= 1<<(unit2^0xf);
		cs->bondOcc[coor[2]] |= 1<<(newUnits[0]|newUnits[1]);
		cs->bondOcc[coor[2]] ^= 1<<unit2;
		
		pol->unitPol[iMono] = newUnits[0];
		pol->unitPol[(iMono+1)%pol->polSize] = newUnits[1];
		pol->coorPol[(iMono+1)%pol->polSize] = (newUnits[0])?coor[2]:coor[0];
// 		CheckIntegrity(cs, "After Backward Linear Move");
		return 1;
	}
	return 1;
}

