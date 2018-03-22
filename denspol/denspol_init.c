#include "denspol.h"

void ReadArgumentsFromFile(SimulationSettings* ss, char* file){
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s (config)\n", file);
		exit(192);
	}
	char strProp[10000], strVal[10000];
	while(fscanf(pFile, "%s %s", strProp, strVal) == 2){
		if(!strcmp(strProp, "L"))
			ss->L = atoi(strVal);
		else if(!strcmp(strProp, "seed"))
			ss->seed = (unsigned int)(atol(strVal));
		else if(!strcmp(strProp, "tmax"))
			ss->tMax = atof(strVal);
		else if(!strcmp(strProp, "npol"))
			ss->nPol = atoi(strVal);
		else if(!strcmp(strProp, "density"))
			ss->density = atof(strVal);
		else if(!strcmp(strProp, "polsize"))
			ss->polSize = atoi(strVal);
		else if(!strcmp(strProp, "dblstep"))
			ss->dblStep = atoi(strVal);
		else if(!strcmp(strProp, "usetopo"))
			ss->useTopo = atoi(strVal);
		else if(!strcmp(strProp, "interval"))
			ss->interval = atof(strVal);
		else if(!strcmp(strProp, "bendenergy"))
			ss->bendEnergy = atof(strVal);
		else if(!strcmp(strProp, "hpstrength"))
			ss->hpStrength = atof(strVal);
		else if(!strcmp(strProp, "dir")){
			ss->dir = malloc(sizeof(char)*(strlen(strVal)+1));
			strcpy(ss->dir, strVal);
		}
		else if(!strcmp(strProp, "eefile")){
			ss->eeFile = malloc(sizeof(char)*(strlen(strVal)+1));
			strcpy(ss->eeFile, strVal);
		}
		else if(!strcmp(strProp, "contactfile")){
			ss->contactFile = malloc(sizeof(char)*(strlen(strVal)+1));
			strcpy(ss->contactFile, strVal);
		}
		else if(!strcmp(strProp, "poltype")){
			if(!strcmp(strVal, "lin") || !strcmp(strVal, "linear"))
				ss->polType = POL_TYPE_LIN;
			else if(!strcmp(strVal, "ring"))
				ss->polType = POL_TYPE_RING;
			else
				ss->polType = POL_TYPE_MIXED;
		}
		else if(!strcmp(strProp, "latticeshape")){
			if(!strcmp(strVal, "sphere"))
				ss->latticeShape = LATTICE_SHAPE_SPHERE;
			else
				ss->latticeShape = LATTICE_SHAPE_EMPTY;
		}
		else if(!strcmp(strProp, "boundarycondition")){
			if(!strcmp(strVal, "periodic"))
				ss->boundaryCond = BOUNDARY_PERIODIC;
			else if(!strcmp(strVal, "static"))
				ss->boundaryCond = BOUNDARY_STATIC;
		}
		else{
			printf("Warning: unknown argument (%s)\n", strProp);
		}
	}
}

void SetDefaultSS(SimulationSettings* ss){
	ss->L            = SS_DEFAULT;
	ss->polSize      = SS_DEFAULT;
	ss->density      = SS_DEFAULT;
	ss->seed         = SS_DEFAULT;
	ss->tMax         = SS_DEFAULT;
	ss->interval     = SS_DEFAULT;
	ss->bendEnergy   = SS_DEFAULT;
	ss->hpStrength   = SS_DEFAULT;
	ss->hpStrength   = SS_DEFAULT;
	ss->boundaryCond = SS_DEFAULT;
	ss->polType      = SS_DEFAULT;
	ss->dblStep      = SS_DEFAULT;
	ss->latticeShape = SS_DEFAULT;
	ss->useTopo      = SS_DEFAULT;
	ss->eeFile       = NULL;
	ss->contactFile  = NULL;
	ss->dir          = NULL;
}

void FillDefaultSS(SimulationSettings* ss){
	if(ss->L            == SS_DEFAULT) ss->L            = 10;
	if(ss->polSize      == SS_DEFAULT) ss->polSize      = 200;
	if(ss->density      == SS_DEFAULT) ss->density      = 7.2;
	if(ss->seed         == SS_DEFAULT) ss->seed         = 94619234;
	if(ss->tMax         == SS_DEFAULT) ss->tMax         = 1e4;
	if(ss->interval     == SS_DEFAULT) ss->interval     = 1e2;
	if(ss->bendEnergy   == SS_DEFAULT) ss->bendEnergy   = 0.3;
	if(ss->hpStrength   == SS_DEFAULT) ss->hpStrength   = 0.35;
	if(ss->boundaryCond == SS_DEFAULT) ss->boundaryCond = BOUNDARY_PERIODIC;
	if(ss->polType      == SS_DEFAULT) ss->polType      = POL_TYPE_RING;
	if(ss->dblStep      == SS_DEFAULT) ss->dblStep      = 0;
	if(ss->latticeShape == SS_DEFAULT) ss->latticeShape = LATTICE_SHAPE_EMPTY;
	if(ss->useTopo      == SS_DEFAULT) ss->useTopo      = 0;
}

void SimulationInit(CurState* cs, LookupTables* lt){
	UnitDotInit(lt, cs->ss.bendEnergy);
	GenerateMutators(lt);
	ReadTopComp(lt, cs->ss.eeFile);
	
	char exec[3000];
	sprintf(exec, "mkdir -p %s", cs->ss.dir);
	system(exec);
	
	long lastT = GetLastT(cs->ss.dir);
	if(lastT>=0){
		printf("Starting from t=%li\n", lastT);
		cs->curT=lastT;
		CSFromFile(cs, lt, lastT);
	}
	else{
		cs->curT=0;
		if(cs->ss.contactFile){
			ReadLengthFile(cs, lt);
		}
		else{
			CSFromParameters(cs,lt);
		}
		GeneratePolymers(cs);
	}
	
	cs->ss.polType = cs->pol[0].polType;
	for(Polymer* pol = cs->pol; pol<cs->pol+cs->nPol; pol++)
		if(pol->polType != cs->ss.polType)
			cs->ss.polType = POL_TYPE_MIXED;
	cs->ss.polSize = cs->maxNMono;
	
	if(cs->ss.contactFile){
		LoadHPFile(cs->ss.contactFile, &lt->hp, cs);
	}
	if(!cs->ss.contactFile){
		lt->hp.nInter = malloc(sizeof(int*)*cs->nPol);
		for(int iPol=0; iPol<cs->nPol; iPol++){
			lt->hp.nInter[iPol] = malloc(sizeof(int)*cs->pol[iPol].nMono);
			for(int iMono=0; iMono<cs->pol[iPol].nMono; iMono++){
				lt->hp.nInter[iPol][iMono]=0;
			}
		}
	}
	if(cs->ss.boundaryCond == BOUNDARY_PERIODIC)
		lt->hp.distance = GenerateDistanceMatrix(cs->L);
	else
		lt->hp.distance = GenerateBoundedDistanceMatrix(cs->L);
	
	lt->nMoveChoice=0;
	PopulateMoveList(cs, lt);
	
	CheckIntegrity(cs, "After construction");
}

void LatticeInit(CurState* cs, LookupTables* lt){
	cs->L = cs->ss.L;
	cs->LSize = cs->ss.L*cs->ss.L*cs->ss.L;
	cs->topoState = malloc(sizeof(int)*cs->LSize);
	cs->bondOcc   = malloc(sizeof(int)*cs->LSize);
	
	if(cs->ss.latticeShape == LATTICE_SHAPE_SPHERE)
		lt->nLatticeUsed = SetLatticeSphere(cs->topoState, cs->bondOcc, cs->L, &lt->nBondUsed);
	else if(cs->ss.latticeShape == LATTICE_SHAPE_EMPTY)
		lt->nLatticeUsed = SetLatticeEmpty(cs->topoState, cs->bondOcc, cs->L, &lt->nBondUsed);
}

void AddInteraction(HPTable* hp, int iMonoOrig, int iPol, int jMonoOrig, int jPol, double strength, CurState* cs){
	strength *= cs->ss.hpStrength;
	Polymer* polI = cs->pol+iPol;
	Polymer* polJ = cs->pol+jPol;
	double magRatioI = MagnificationRatio(polI->nMono, polI->origNMono, polI->polType);
	double magRatioJ = MagnificationRatio(polJ->nMono, polJ->origNMono, polJ->polType);
	
	int iMono = (int)(iMonoOrig*magRatioI+0.5)%polI->nMono;
	int jMono = (int)(jMonoOrig*magRatioJ+0.5)%polJ->nMono;
	
	if(iMono == jMono && iPol == jPol) return;
	int exist=0;
	for(int indexI=0; indexI<hp->nInter[iPol][iMono] && !exist; indexI++){
		TripleHP* interaction = &hp->inter[iPol][iMono][indexI];
		if(interaction->iMono == jMono && interaction->iPol == jPol){
			interaction->strength += strength;
			exist=1;
		}
	}
	if(!exist){
		hp->inter[iPol][iMono][hp->nInter[iPol][iMono]  ].iMono     = jMono;
		hp->inter[iPol][iMono][hp->nInter[iPol][iMono]  ].iPol      = jPol;
		hp->inter[iPol][iMono][hp->nInter[iPol][iMono]++].strength  = strength;
		
		hp->inter[jPol][jMono][hp->nInter[jPol][jMono]  ].iMono     = iMono;
		hp->inter[jPol][jMono][hp->nInter[jPol][jMono]  ].iPol      = iPol;
		hp->inter[jPol][jMono][hp->nInter[jPol][jMono]++].strength  = strength;
	}
	else{
		for(int indexJ=0; indexJ<hp->nInter[jPol][jMono]; indexJ++){
			TripleHP* interaction = &hp->inter[jPol][jMono][indexJ];
			if(interaction->iMono == iMono && interaction->iPol == iPol){
				interaction->strength += strength;
				break;
			}
		}
	}
}

void LoadHPFile(char* file, HPTable* hp, CurState* cs){
// 	int polSize = cs->polSize;
	FILE* pFile = fopen(file, "r");
	int maxNMono;
	int nPol;
	long nContacts=-1;
	
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	fscanf(pFile, "%*s %i", &nPol);
	fscanf(pFile, "%*s %i", &maxNMono);
	fscanf(pFile, "%*s %li", &nContacts); 
	
// 	printf("N = (%i %i)\n", cs->nMono, nMono);
	if(nPol != cs->nPol){
		printf("Error: wrong number of polymers in contact file (%i vs assumed %i)\n", nPol, cs->nPol);
		exit(192);
	}
	
	///Assume that linear/ring configuration is okay...
	for(int iPol=0; iPol<nPol; iPol++){
		fscanf(pFile, "%*s %i", &cs->pol[iPol].origNMono);
	}
	double ratio = cs->pol[0].nMono/(double)cs->pol[0].origNMono;
	int maxInter = MAX(12, 12/(ratio*ratio));
	
	hp->inter    = (TripleHP***)malloc(sizeof(TripleHP**)*cs->nPol);
	hp->nInter   = (int**)      malloc(sizeof(int*)      *cs->nPol);
	for(int iPol=0; iPol<cs->nPol; iPol++){
		hp->inter[iPol]  = malloc(sizeof(TripleHP*)*cs->pol[iPol].nMono);
		hp->nInter[iPol] = malloc(sizeof(int)      *cs->pol[iPol].nMono);
		for(int iMono=0; iMono<cs->pol[iPol].nMono; iMono++){
			hp->inter[iPol][iMono] = malloc(sizeof(TripleHP)*maxInter);
			hp->nInter[iPol][iMono]=0;
		}
	}
	int iMono, jMono, iPol, jPol;
	double strength;
	long nContactsFound=0;
	while( !(feof(pFile)) && fscanf(pFile, "%i %i %i %i %lf", &iPol, &iMono, &jPol, &jMono, &strength) == 5){
		AddInteraction(hp, iMono, iPol, jMono, jPol, strength, cs);
		nContactsFound++;
	}
	
	if(nContactsFound != nContacts){
		printf("Error in the number of contacts in file %s\n%li vs %li\n", file, nContacts, nContactsFound);
		exit(192);
	}
	
	double strPrefac = MIN(1, 0.874*cs->nTotMono/(double)nContacts);
	
	for(int iPol=0; iPol<cs->nPol; iPol++){
		for(int iMono=0; iMono<cs->pol[iPol].nMono; iMono++){
			for(int iContact=0; iContact<hp->nInter[iPol][iMono]; iContact++){
				hp->inter[iPol][iMono][iContact].strength *= strPrefac;
			}
		}
	}
}

void AllocPolymers(CurState* cs){
	cs->pol       = malloc(sizeof(Polymer)*cs->nPol);
	for(int iPol=0; iPol<cs->nPol; iPol++){
		cs->pol[iPol].unitPol = malloc(sizeof(int)*(cs->maxNMono));
		cs->pol[iPol].coorPol = malloc(sizeof(int)*(cs->maxNMono));
		for(int iMono=0; iMono<cs->maxNMono; iMono++){
			cs->pol[iPol].unitPol[iMono] = 0xf;
			cs->pol[iPol].coorPol[iMono] = -1;
		}
	}
}

int CSFromParameters(CurState* cs, LookupTables* lt){
	FillDefaultSS(&cs->ss);
	Seed(cs->rngState, cs->ss.seed);
	if(cs->ss.boundaryCond == BOUNDARY_PERIODIC)
		cs->AddUnitToCoor = &AddUnitToCoorPeriod;
	else if(cs->ss.boundaryCond == BOUNDARY_STATIC)
		cs->AddUnitToCoor = &AddUnitToCoorWBounds;
	LatticeInit(cs, lt);
	
	int nMono = cs->ss.polSize+((cs->ss.polType == POL_TYPE_LIN)?1:0);
	cs->nPol = (int)((cs->ss.density/6.0)*lt->nBondUsed/(double)nMono+0.5);
	cs->maxNMono = nMono;
	cs->nTotMono = cs->nPol*nMono;
	AllocPolymers(cs);
	for(int iPol=0; iPol<cs->nPol; iPol++){
		Polymer* pol = cs->pol+iPol;
		pol->polType = cs->ss.polType;
		pol->nMono = nMono;
		pol->polSize = cs->ss.polSize;
	}
	return 0;
}

void PopulateMoveList(CurState* cs, LookupTables* lt){
	
	if(!lt->nMoveChoice){
		lt->nMoveChoice=0;
		lt->nLatticeUsed=0;
		for(int coor=0; coor<cs->LSize; coor++){
			if(cs->topoState[coor] >= 0)
				lt->nLatticeUsed++;
		}
		
		if(cs->ss.useTopo)
			lt->nMoveChoice += lt->nLatticeUsed;
		
		for(Polymer* pol = cs->pol; pol<cs->pol+cs->nPol; pol++)
			lt->nMoveChoice += 2*pol->polSize;
		
		lt->moveChoiceTable = malloc(sizeof(int*)*lt->nMoveChoice);
		for(int i=0; i<lt->nMoveChoice; i++)
			lt->moveChoiceTable[i] = malloc(sizeof(int)*3);
		lt->latticeUsed = malloc(sizeof(int)*lt->nLatticeUsed);
	}
	
	int curMove=0;
	int curCoor=0;
	if(cs->ss.useTopo){
		for(int coor=0; coor<cs->LSize; coor++){
			if(cs->topoState[coor] >= 0){
				lt->moveChoiceTable[curMove  ][0] = MOVE_TOPO;
				lt->moveChoiceTable[curMove  ][1] = NON_EXISTING;
				lt->moveChoiceTable[curMove++][2] = NON_EXISTING;
				lt->latticeUsed[curCoor++] = coor;
			}
		}
	}
	for(Polymer* pol = cs->pol; pol<cs->pol+cs->nPol; pol++){
		int iPol = pol-cs->pol;
		if(pol->polType == POL_TYPE_LIN){
			for(int iMono=1; iMono<pol->nMono-1; iMono++){
				lt->moveChoiceTable[curMove  ][0] = MOVE_TRANS_LIN;
				lt->moveChoiceTable[curMove  ][1] = iPol;
				lt->moveChoiceTable[curMove++][2] = iMono;
				
				lt->moveChoiceTable[curMove  ][0] = MOVE_DIFFUSE;
				lt->moveChoiceTable[curMove  ][1] = iPol;
				lt->moveChoiceTable[curMove++][2] = iMono;
			}
			lt->moveChoiceTable[curMove  ][0] = MOVE_START;
			lt->moveChoiceTable[curMove  ][1] = iPol;
			lt->moveChoiceTable[curMove++][2] = NON_EXISTING;
			
			lt->moveChoiceTable[curMove  ][0] = MOVE_END;
			lt->moveChoiceTable[curMove  ][1] = iPol;
			lt->moveChoiceTable[curMove++][2] = NON_EXISTING;
		}
		else{
			for(int iMono=0; iMono<pol->nMono; iMono++){
				lt->moveChoiceTable[curMove  ][0] = MOVE_TRANS_RING;
				lt->moveChoiceTable[curMove  ][1] = iPol;
				lt->moveChoiceTable[curMove++][2] = iMono;
				
				lt->moveChoiceTable[curMove  ][0] = MOVE_DIFFUSE;
				lt->moveChoiceTable[curMove  ][1] = iPol;
				lt->moveChoiceTable[curMove++][2] = iMono;
			}
		}
	}
}

int CheckPolymerPartition(CurState* cs, int delta){
	int nFound=0;
	
	int units[] = {0x0, 0x4, 0x5};
	for(int t=0; t+delta/2<cs->L; t+=delta){
		for(int u=0; u+delta/2<cs->L; u+=delta){
			for(int v=0; v+delta/2<cs->L; v+=delta){
				int coor = TUV2Coor(t,u,v,cs->L);
				int OOB;
				int available=1;
				for(int iUnit=0; iUnit<3 && available; iUnit++){
					int newCoor = cs->AddUnitToCoor(units[iUnit], coor, cs->L, &OOB);
					if(OOB || cs->topoState[newCoor] <0) available = 0; 
				}
				nFound += available;
			}
		}
	}
	return nFound;
}

int CheckLinearPolymerPartition(CurState* cs, int delta){
	int nFound=0;
	
	for(int t=0; t+delta/2<cs->L; t+=delta){
		for(int u=0; u+delta/2<cs->L; u+=delta){
			for(int v=0; v+delta/2<cs->L; v+=delta){
				int coor = TUV2Coor(t,u,v,cs->L);
				int OOB;
				if(cs->topoState[coor] < 0) continue;
				
				for(int unit=0x1; unit<=0x7; unit++){
					if(!IsValid(unit)) continue;
					int newCoor = cs->AddUnitToCoor(unit, coor, cs->L, &OOB);
					if(!(OOB || cs->topoState[newCoor] <0))
						nFound++;
				}
			}
		}
	}
	return nFound;
}

int* GenerateRandomPolOrder(int nPol, unsigned int rng[4]){
	int* polOrder = malloc(sizeof(int)*nPol);
	for(int i=0; i<nPol; i++) polOrder[i]=i;
	int tPol;
	for(int i=0; i<nPol-1; i++){
		int j= i+DRng(rng)*(nPol-i);
		if(i != j){
			tPol = polOrder[i];
			polOrder[i] = polOrder[j];
			polOrder[j] = tPol;
		}
	}
	return polOrder;
}

void GenerateLinearPolymers(CurState* cs){
	int delta;
	int L = cs->L;
	
	for(delta=L; delta>=1; delta--){
		int nTot = 1;
		nTot *= L/delta;
		nTot *= L/delta;
		nTot *= L/delta;
		if(nTot >= cs->nPol && CheckLinearPolymerPartition(cs, delta) >= cs->nPol){
			break;
		}
	}
	if(!delta){
		printf("Error finding polymer partition, too many polymers/too dense (linear): got %i vs need %i.\n", CheckLinearPolymerPartition(cs, 1), cs->nPol);
		exit(192);
	}
	int iPol=0;
	int* polOrder = GenerateRandomPolOrder(cs->nPol, cs->rngState);
	
	for(int unit=0x1; unit<0x8; unit++){
		for(int t=0; t+delta/2<L && iPol<cs->nPol; t += delta){
			for(int u=0; u+delta/2<L && iPol<cs->nPol; u += delta){
				for(int v=0; v+delta/2<L && iPol<cs->nPol; v += delta){
					Polymer* pol = cs->pol+polOrder[iPol];
					
					int OOB;
					int available=1;
					int coor = TUV2Coor(t,u,v,cs->L);
					if(cs->topoState[coor] < 0) continue;
					int newCoor = cs->AddUnitToCoor(unit, coor, cs->L, &OOB);
					if(OOB || cs->topoState[newCoor] <0) available = 0; 
					
					if(!available) continue;
					
					pol->unitPol[0] = unit;
					for(int iMono=1; iMono<pol->polSize; iMono++)
						pol->unitPol[iMono] = 0x0;
					pol->unitPol[pol->polSize] = 0xf;
					
					for(int iMono=0; iMono<=pol->polSize; iMono++){
						pol->coorPol[iMono] = coor;
						coor = cs->AddUnitToCoor(pol->unitPol[iMono], coor, L, &OOB);
					}
					
					int bondOcc0 = 1<<(pol->unitPol[0]^0xf);
					int bondOcc1 = 1<<(pol->unitPol[0]);
					
					cs->bondOcc[pol->coorPol[0]] |= bondOcc0;
					cs->bondOcc[pol->coorPol[1]] |= bondOcc1;
					
					iPol++;
				}
			}
		}
	}
	free(polOrder);
}

void GeneratePolymers(CurState* cs){
	
	int onlyLinear=1;
	for(int iPol=0; iPol<cs->nPol && onlyLinear; iPol++){
		if(cs->pol[iPol].polType != POL_TYPE_LIN) onlyLinear=0;
	}
	
	if(onlyLinear){
		GenerateLinearPolymers(cs);
		return;
	}
	
	int delta;
	int L = cs->L;
	
	for(delta=L; delta>=1; delta--){
		int nTot = 1;
		nTot *= L/delta;
		nTot *= L/delta;
		nTot *= L/delta;
		if(nTot >= cs->nPol && CheckPolymerPartition(cs, delta) >= cs->nPol){
			break;
		}
	}
	if(!delta){
		printf("Error finding polymer partition, too many polymers/too dense.\n");
		exit(192);
	}
	
	int availUnits[]= {0x0,0x4,0x5};
	int ringUnits[] = {0x4,0x1,0xa};
	int linUnits[]  = {0x4};
	
	
	int* polOrder = GenerateRandomPolOrder(cs->nPol, cs->rngState);
	int iPol=0; 
	Polymer* pol;
	for(int t=0; t+delta/2<L && iPol<cs->nPol; t += delta){
		for(int u=0; u+delta/2<L && iPol<cs->nPol; u += delta){
			for(int v=0; v+delta/2<L && iPol<cs->nPol; v += delta){
				pol = cs->pol+polOrder[iPol];
				int OOB;
				int available=1;
				int coor = TUV2Coor(t,u,v,cs->L);
				for(int iUnit=0; iUnit<3 && available; iUnit++){
					int newCoor = cs->AddUnitToCoor(availUnits[iUnit], coor, cs->L, &OOB);
					if(OOB || cs->topoState[newCoor] <0) available = 0; 
				}
				if(!available) continue;
// 				printf("polType = %i\n", pol->polType);
				
				if(pol->polType == POL_TYPE_RING){
					for(int j=0; j<3; j++)
						pol->unitPol[j] = ringUnits[j];
					for(int iMono=3; iMono<pol->polSize; iMono++)
						pol->unitPol[iMono] = 0;
					int coor = TUV2Coor(t,u,v,L);
					for(int iMono=0; iMono<pol->polSize; iMono++){
						pol->coorPol[iMono] = coor;
						coor = cs->AddUnitToCoor(pol->unitPol[iMono], coor, L, &OOB);
						if(OOB) printf("Error: out of bounds...\n");
					}
					
					for(int iMono=0; iMono<3; iMono++){
						int bondOcc = 1<<pol->unitPol[(iMono-1+3)%3];
						bondOcc |= 1<<(pol->unitPol[iMono]^0xf);
						if(cs->bondOcc[pol->coorPol[iMono]] & bondOcc){
							printf("Uh oh, iMono = %i\n", iMono);
							printf("at (t,u,v): "); PrintCoor(pol->coorPol[iMono], cs->L);
							printf("\ntopo = %i, bondOcc = %i\n", cs->topoState[pol->coorPol[iMono]], cs->bondOcc[pol->coorPol[iMono]]);
							exit(192);
						}
						cs->bondOcc[pol->coorPol[iMono]] |= bondOcc;
					}
					iPol++;
				}
				else if(pol->polType == POL_TYPE_LIN){
					pol->unitPol[0] = linUnits[0];
				
					for(int iMono=1; iMono<pol->polSize; iMono++)
						pol->unitPol[iMono] = 0;
					pol->unitPol[pol->polSize] = 0xf;
					
					int coor = TUV2Coor(t,u,v,L);
					for(int iMono=0; iMono<=pol->polSize; iMono++){
						pol->coorPol[iMono] = coor;
						coor = cs->AddUnitToCoor(pol->unitPol[iMono], coor, L, &OOB);
					}
					
					int bondOcc0 = 1<<(pol->unitPol[0]^0xf);
					int bondOcc1 = 1<<(pol->unitPol[0]);
					
					cs->bondOcc[pol->coorPol[0]] |= bondOcc0;
					cs->bondOcc[pol->coorPol[1]] |= bondOcc1;
					
					iPol++;
				}
			}
		}
	}
	free(polOrder);
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

double*** GenerateBoundedDistanceMatrix(int L){
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
	
	double tuv[3];
	double xyz[3];
	double xyzZero[3]={0,0,0};
	
	for(int t=-L; t<L; t++){
		for(int u=-L; u<L; u++){
			for(int v=-L; v<L; v++){
				tuv[0]=t; tuv[1]=u; tuv[2]=v;
				DTUV2XYZ(tuv,xyz);
				distances[t][u][v] = DistanceSq(xyz, xyzZero);
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
	
	return distances;
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



