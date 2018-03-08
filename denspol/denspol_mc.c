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

int DiffuseStep(CurState* cs, LookupTables* lt, Polymer* pol, int iMono){
	int iPrev = (iMono-1+pol->nMono)%pol->nMono;
	if(pol->unitPol[iMono] && !pol->unitPol[iPrev]){
#ifdef __HP_ENABLED__
		int unitMove = pol->unitPol[iMono];
		if(!TestMoveHP(cs,lt,iMono, pol, unitMove))
			return 0;
#endif
		pol->unitPol[iPrev] = pol->unitPol[iMono];
		pol->unitPol[iMono] = 0;
		pol->coorPol[iMono] = pol->coorPol[(iMono+1)%pol->nMono];
// 		char msg[2000];
// 		sprintf(msg, "After Diffuse Move(a), iMono = %i, iPol = %li\n", iMono, pol-cs->pol);
// 		CheckIntegrity(cs, msg);
		return 1;
	}
	else if(!pol->unitPol[iMono] && pol->unitPol[iPrev]){
#ifdef __HP_ENABLED__
		int unitMove = pol->unitPol[iPrev]^0xf;
		if(!TestMoveHP(cs,lt,iMono, pol, unitMove)) 
			return 0;
#endif
		pol->unitPol[iMono] = pol->unitPol[iPrev];
		pol->unitPol[iPrev] = 0;
		pol->coorPol[iMono] = pol->coorPol[iPrev];
// 		char msg[2000];
// 		sprintf(msg, "After Diffuse Move(b), iMono = %i, iPol = %li\n", iMono, pol-cs->pol);
// 		CheckIntegrity(cs, msg);
		return 1;
	}
	return 0;
}

/// Time reversibility is broken here, so if you care, change it!

int TopoMove(CurState* cs, LookupTables* lt){
	int coor = lt->latticeUsed[(int)(lt->nLatticeUsed*DRng(cs->rngState))];
	
	int ret=-1;
	
	if(lt->topComp[cs->topoState[coor]].sameTopo != cs->topoState[coor])
		ret=1;
	else
		ret=0;
	
	cs->topoState[coor] = lt->topComp[cs->topoState[coor]].sameTopo;
// 	CheckIntegrity(cs, "After Topo Move");
	return ret;
}


double DoMCStep(long nStep, CurState* cs, LookupTables* lt){
	for(long iStep=0; iStep<nStep*lt->nMoveChoice; iStep++){
		int iMove= DRng(cs->rngState)*lt->nMoveChoice; 
		int move = lt->moveChoiceTable[iMove][0];
		Polymer* pol = cs->pol+lt->moveChoiceTable[iMove][1];
		int iMono= lt->moveChoiceTable[iMove][2];
		if(move == MOVE_DIFFUSE)
			DiffuseStep(cs, lt, pol, iMono);
		else if (move == MOVE_TRANS_RING)
			TransMoveRing(cs, lt, pol, iMono);
		else if (move == MOVE_TRANS_LIN)
			TransMoveLinear(cs, lt, pol, iMono);
		else if (move == MOVE_TOPO)
			TopoMove(cs, lt);
		else if (move == MOVE_START)
			StartMove(cs, lt, pol);
		else if (move == MOVE_END)
			EndMove(cs, lt, pol);
		else
			printf("Error; invalid move\n");
	}
	return 0;
}