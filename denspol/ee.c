#include "ee.h"



int main(int argc, char** argv){
	EEInit(20);
	DirectMutateTopo(&ee);
	WriteSuperTopo(&ee, "ee_topo.dat");
	return 0;
}

void EEInit(int maxDepth){
	ee.nConfig=0;
	ee.nNodes=0;
	ee.maxTopo=-1;
	ee.maxDepth = maxDepth;
	ee.mt = NewMoveTable();
	ee.lattice = NewLattice(5);
	LatticeSetHole(ee.lattice, ee.mt);
	SetCornerMoves(&ee);
	GenerateMutators(&ee);
// 	PrintMutators(&ee);
	EETopoInit(&ee);
}

void LatticeSetHole(Lattice* lattice, MoveTable* mt){
	for(int i=0; i<lattice->LS; i++) SetLattice(i, 1, lattice);
	for(int i=0; i<lattice->LS; i++) lattice->labels[i]=-1;
	
	SetMidLattice(0, 0, lattice);
	lattice->mLabels[0] = 0;
	
	for(int iUnit=0; iUnit<12; iUnit++){
		int unit = mt->validUnits[iUnit];
		int pos  = AddUnitToCoor(unit, 0, lattice);
		SetMidLattice(pos, 0, lattice);
		lattice->mLabels[pos] = iUnit+1;
	}
}

void SetCornerMoves(ExactEnum* ee){
	for(int unitA=1; unitA<0xf; unitA++){
		for(int unitB=1; unitB<0xf; unitB++){
			int unitC = ee->mt->backMoves[unitA][unitB];
			if(!unitC) continue;
			int posA = AddUnitToCoor((~unitA)&0xf, 0, ee->lattice);
			int posB = AddUnitToCoor(unitB, 0, ee->lattice);
			ee->cornerMoves[ee->nCornerMoves  ][0] = posA;
			ee->cornerMoves[ee->nCornerMoves  ][1] = posB;
			ee->cornerMoves[ee->nCornerMoves  ][2] = unitA;
			ee->cornerMoves[ee->nCornerMoves  ][3] = unitB;
			ee->cornerMoves[ee->nCornerMoves++][4] = unitC;
		}
	}
}

void GenerateMutators(ExactEnum* ee){
	///First generate double mutators
	
	ee->nMutator = 0;
	
	for(int i=0; i<ee->nCornerMoves; i++){
		for(int initial=0; initial<2; initial++){
			for(int k=0; k<2; k++)
				ee->allMutator[ee->nMutator].move[k] = ee->cornerMoves[i][k];
			ee->allMutator[ee->nMutator].bond = ee->cornerMoves[i][4];
			ee->allMutator[ee->nMutator++].initial = initial;
		}
	}
}


void EETopoInit(ExactEnum* ee){
	Topo* dstTopo=NULL;
	Config* startCfg = NewConfig();
	startCfg->nPol=0;
	startCfg->topo=ee->allTopo;
	ee->nTopo=0;
	int inserted;
// 	PrintConfig(startCfg);
	ee->allTopo = malloc(sizeof(Topo)*MAX_TOPO);
	ee->allConfigs = InsertConfig(NULL, 0, startCfg, &inserted, &dstTopo);
	dstTopo->allConfigs = InsertConfig(NULL, 0, startCfg, &inserted, &dstTopo);
	DeleteConfig(startCfg);
	if(!ee->allConfigs) printf("???\n");
}

Topo* NewTopo(){
	if(ee.nTopo>=MAX_TOPO){
		printf("Error: Not enough memory allocated for topo\n");
		exit(192);
	}
	Topo* topo = ee.allTopo+ee.nTopo;
	topo->supTopo = malloc(sizeof(SuperTopo));
	for(int i=0; i<NMUTATOR; i++) topo->supTopo->mutTopo[i]=NULL;
	topo->supTopo->topo = topo;
	topo->supTopo->id=-1;
	topo->supTopo->level=-1;
	topo->topoNext = NULL;
	topo->allConfigs=NULL;
	topo->mutated=0;
	ee.nTopo++;
	return topo;
}



void PrintMutators(ExactEnum* ee){
	printf("Number of mutators: %i\n", ee->nMutator);
	
	for(int i=0; i<ee->nMutator; i++){
		Mutator* mut = ee->allMutator+i;
		printf("#[%i] %i -> %i [0x%x]",i, mut->move[0], mut->move[1], mut->bond);
		printf("\n");
	}
	printf("\n");
}

void PrintCornerMoves(ExactEnum* ee){
	printf("Number of corner moves: %i\n", ee->nCornerMoves);
	for(int i=0; i<ee->nCornerMoves; i++){
		for(int j=0; j<3; j++) printf("%i ", ee->cornerMoves[i][j]);
		printf("\n");
	}
}

void CopyConfig(Config* src, Config* dst){
	dst->nPol = src->nPol;
	for(int iPol=0; iPol<src->nPol; iPol++){
		dst->nBonds[iPol] = src->nBonds[iPol];
		for(int iBond=0; iBond<src->nBonds[iPol]; iBond++){
			dst->bonds[iPol][iBond] = src->bonds[iPol][iBond];
			dst->pos[iPol][iBond] = src->pos[iPol][iBond];
		}
		dst->pos[iPol][src->nBonds[iPol]] = src->pos[iPol][src->nBonds[iPol]];
		dst->topo = src->topo;
	}
}

Config* NewFromConfig(Config* src){
	Config* dst = NewConfig();
	CopyConfig(src,dst);
	return dst;
}

Node* NewNode(int val){
	Node* node = malloc(sizeof(Node));
	node->val = val;
	node->topo = NULL;
	node->next = NULL;
	node->child = NULL;
	ee.nNodes++;
	return node;
}

void DeleteNode(Node* node){
	if(node){
		free(node);
		ee.nNodes--;
	}
}

void DeleteNodeTree(Node* node){
	if(!node) return;
	DeleteNodeTree(node->next);
	DeleteNodeTree(node->child);
	DeleteNode(node);
}

Config* NewConfig(){
	ee.nConfig++;
	Config* cfg = malloc(sizeof(Config));
	cfg->topo = NULL;
	cfg->nPol = 0;
	for(int i=0; i<6; i++) cfg->nBonds[i]=-1;
	return cfg;
}

void DeleteConfig(Config* cfg){
	free(cfg);
	ee.nConfig--;
}

int CheckKey(int* key, int nKey){
	for(int i=0; i<nKey; i++){
		for(int j=0; j<nKey; j++){
			if(i!=j && key[i] != END && key[i] == key[j]){
				printf("Error setting key, double occupancy\n");
				return 1;
			}
		}
	}
	return 0;
}

void SetConfigKey(Config* cfg){
	cfg->nKey=0;
	int previous = -129386401;
	for(int iPol=0; iPol<cfg->nPol; iPol++){
		int lowest = 13984619;
		int bestPol = -1;
		for(int jPol=0; jPol<cfg->nPol; jPol++){
			if(cfg->pos[jPol][0] > previous && cfg->pos[jPol][0] < lowest){
				bestPol = jPol;
				lowest = cfg->pos[jPol][0];
			}
		}
		previous = lowest;
		if(bestPol<0){ printf("Error setting key\n"); exit(192); }
		
		for(int iMono=0; iMono<=cfg->nBonds[bestPol]; iMono++)
			cfg->key[cfg->nKey++] = cfg->pos[bestPol][iMono];
		cfg->key[cfg->nKey++] = END;
	}
	if(!cfg->nPol) cfg->key[cfg->nKey++] = END;
// 	if(CheckKey(cfg->key, cfg->nKey)){
// 		PrintConfig(cfg);
// 		printf("Key:\n\n");
// 		for(int i=0; i<cfg->nKey; i++)
// 			printf("%i ", cfg->key[i]);
// 		printf("\n\n");
// 		exit(0);
// 	}
}

Node* InsertConfig(Node* node, int iKey, Config* cfg, int* inserted, Topo** destTopo){
	if(!iKey){
		SetConfigKey(cfg);
		*inserted=0;
	}
	
	int val = cfg->key[iKey];
	if(!node){ node = NewNode(val); }
	
	if(node->val == val){
		if(iKey == cfg->nKey-1){
			if(! node->topo){
				if(!*destTopo){
					node->topo = NewTopo();
				}
				else
					node->topo = *destTopo;
				*inserted=1;
			}
			*destTopo = node->topo;
		}
		else{
			node->child = InsertConfig(node->child, iKey+1, cfg, inserted, destTopo);
		}
	}
	else{
		node->next = InsertConfig(node->next, iKey, cfg, inserted, destTopo);
	}
	return node;
}

// Topo* TopoGetRoot(Topo* topo){
// 	if(topo->mergedWith)
// 		return TopoGetRoot(topo->mergedWith);
// 	else
// 		return topo;
// }

void MergeSuperTopo(SuperTopo* sTopoA, SuperTopo* sTopoB){
	Topo* topo = sTopoA->topo;
	
	topo->supTopo = sTopoB;
	while(topo->topoNext){
		topo = topo->topoNext;
		topo->supTopo = sTopoB;
	}
	topo->topoNext = sTopoB->topo;
	sTopoB->topo = sTopoA->topo;

	if(topo->topoNext == topo){
		printf("?????\n");
		exit(0);
	}
	
	for(int iMut=0; iMut<NMUTATOR; iMut++){
		if(!sTopoA->mutTopo[iMut]) continue;
		if(!sTopoB->mutTopo[iMut]){ sTopoB->mutTopo[iMut] = sTopoA->mutTopo[iMut]; continue; }
		SuperTopo* destTopoA = sTopoA->mutTopo[iMut]->supTopo;
		SuperTopo* destTopoB = sTopoB->mutTopo[iMut]->supTopo;
		if(destTopoA != destTopoB)
			MergeSuperTopo(destTopoA, destTopoB);
	}
	
	free(sTopoA);
}

int FindInternalCfg(ExactEnum* ee, Config* cfg, Topo** topo){
	MoveTable* mt = ee->mt;
	int totInsert=0;
	int dummy;
	ee->allConfigs = InsertConfig(ee->allConfigs, 0, cfg, &totInsert, topo);
	
	if(!totInsert){ /*printf("config already exists\n");*/ return 0; }
	if(!(*topo)){printf("ooops\n"); exit(192);}
	(*topo)->allConfigs = InsertConfig((*topo)->allConfigs, 0, cfg, &dummy, topo);
	
	for(int iPol=0; iPol<cfg->nPol; iPol++){
		for(int iBond=0; iBond<cfg->nBonds[iPol]; iBond++){
			int bond = cfg->bonds[iPol][iBond];
			if(!IsValid(bond)){ printf("Woooops\n"); PrintConfig(cfg); exit(0); }
			/// Try forward moves.
			for(int iForw=0; iForw<4; iForw++){
				int newBond1 = mt->forwMoves[bond][iForw][0];
				int newBond2 = mt->forwMoves[bond][iForw][1];
				int newPos = AddUnitToCoor(newBond1, cfg->pos[iPol][iBond], ee->lattice);
				
				if(ConfigOccupied(newPos, cfg) || GetMidLattice(newPos, ee->lattice)) continue;
				Config* newCfg = NewFromConfig(cfg);
				InsertValue(newPos, iBond+1, cfg->nBonds[iPol]+1, newCfg->pos[iPol]);
				InsertValue(newBond2, iBond+1, cfg->nBonds[iPol], newCfg->bonds[iPol]);
				newCfg->nBonds[iPol]++;
				newCfg->bonds[iPol][iBond] = newBond1;
				
				int newInsert = FindInternalCfg(ee, newCfg, topo);
				DeleteConfig(newCfg);
				totInsert+=newInsert;
			}
			/// Try backward moves.
			if((iBond>=cfg->nBonds[iPol]-1)) continue;
			int newBond = mt->backMoves[bond][cfg->bonds[iPol][iBond+1]];
			if(!newBond) continue;
			
			Config* newCfg = NewFromConfig(cfg);
			for(int i=iBond+1; i<cfg->nBonds[iPol]-1; i++){
				newCfg->bonds[iPol][i] = cfg->bonds[iPol][i+1];
				newCfg->pos[iPol][i] = cfg->pos[iPol][i+1];
			}
			newCfg->pos[iPol][cfg->nBonds[iPol]-1] = cfg->pos[iPol][cfg->nBonds[iPol]];
			newCfg->bonds[iPol][iBond] = newBond;
			newCfg->nBonds[iPol]--;
			
			int newInsert = FindInternalCfg(ee, newCfg, topo);
			DeleteConfig(newCfg);
			totInsert+= newInsert;
		}
	}
	return totInsert;
}

int CheckIntegrity(Config* cfg){
	for(int iPol=0; iPol<cfg->nPol; iPol++){
		for(int iBond=0; iBond<=cfg->nBonds[iPol]; iBond++){
			if(GetMidLattice(cfg->pos[iPol][iBond], ee.lattice)){
				printf("Error: pos occupied\n");
				exit(0);
			}
			if(iBond<cfg->nBonds[iPol]){
				if(!IsValid(cfg->bonds[iPol][iBond])){
					printf("Error: bond invalid\n");
					exit(0);
				}
				
				if(AddUnitToCoor(cfg->bonds[iPol][iBond], cfg->pos[iPol][iBond], ee.lattice) != cfg->pos[iPol][iBond+1]){
					printf("Error: bond+pos != newPos, (%i,%i)\n", iPol, iBond);
					PrintConfig(cfg);
					exit(0);
				}
			}
			for(int jPol=0; jPol<cfg->nPol; jPol++){
				for(int jBond=0; jBond<=cfg->nBonds[jPol]; jBond++){
					if((iPol != jPol || iBond != jBond) && cfg->pos[iPol][iBond] == cfg->pos[jPol][jBond]){
						printf("Double occupancy!\n");
						exit(0);
					}
				}
			}
		}
	}
	return 0;
}

void InsertBond(int bond, int iBond, int iPol, Config* cfg){
	InsertValue(bond, iBond, cfg->nBonds[iPol], cfg->bonds[iPol]);
}

void InsertPos(int pos, int iPos, int iPol, Config* cfg){
	InsertValue(pos, iPos, cfg->nBonds[iPol]+1, cfg->pos[iPol]);
}

void InsertValue(int val, int iPos, int len, int* arr){
	for(int i=len; i>iPos; i--){
		arr[i]=arr[i-1];
	}
	arr[iPos] = val;
}

void DeleteValue(int iPos, int len, int*arr){
}

int ConfigOccupied(int pos, Config* cfg){
	int occupied=0;
	for(int iPol=0; iPol<cfg->nPol && !occupied; iPol++){
		for(int iMono=0; iMono<=cfg->nBonds[iPol] && !occupied; iMono++){
			if(pos == cfg->pos[iPol][iMono]) occupied = 1;
		}
	}
	return occupied;
}

Config* MutateConfig(Config* cfg, Mutator* mt){
	Config* newCfg = NULL;
	int occ[2];
	occ[0] = ConfigOccupied(mt->move[0], cfg);
	occ[1] = ConfigOccupied(mt->move[1], cfg);
	
	if(mt->initial){
		if( ! (occ[0] || occ[1]) ){
			newCfg = NewFromConfig(cfg);
			newCfg->pos[newCfg->nPol][0] = mt->move[0];
			newCfg->pos[newCfg->nPol][1] = mt->move[1];
			newCfg->bonds[newCfg->nPol][0] = mt->bond;
			newCfg->nBonds[newCfg->nPol++] = 1;
			return newCfg;
		}
		else if( occ[0] && occ[1] ) {
			for(int iPol=0; iPol<cfg->nPol; iPol++){
				if(cfg->pos[iPol][0] == mt->move[0] && cfg->pos[iPol][1] == mt->move[1] && cfg->nBonds[iPol] == 1){
					newCfg = NewFromConfig(cfg);
					int jPol = newCfg->nPol-1;
					if(jPol != iPol){
						for(int iMono=0; iMono<=cfg->nBonds[jPol]; iMono++) 
							newCfg->pos[iPol][iMono] = newCfg->pos[jPol][iMono];
						for(int iBond=0; iBond<cfg->nBonds[jPol]; iBond++) 
							newCfg->bonds[iPol][iBond] = newCfg->bonds[jPol][iBond];
						newCfg->nBonds[iPol] = newCfg->nBonds[jPol];
					}
					newCfg->nPol--;
					return newCfg;
				}
			}
		}
	}
	else{
		if ( occ[0] && ! occ[1] ){
			for(int iPol=0; iPol<cfg->nPol; iPol++){
				if(cfg->pos[iPol][0] == mt->move[0]){
					newCfg = NewFromConfig(cfg);
					InsertValue(mt->move[1], 0, newCfg->nBonds[iPol]+1, newCfg->pos[iPol]);
					InsertValue((~mt->bond)&0xf, 0, newCfg->nBonds[iPol], newCfg->bonds[iPol]);
					newCfg->nBonds[iPol]++;
					return newCfg;
				}
				if(cfg->pos[iPol][cfg->nBonds[iPol]] == mt->move[0]){
					newCfg = NewFromConfig(cfg);
					InsertValue(mt->move[1], newCfg->nBonds[iPol]+1, newCfg->nBonds[iPol]+1, newCfg->pos[iPol]);
					InsertValue(mt->bond, newCfg->nBonds[iPol], newCfg->nBonds[iPol], newCfg->bonds[iPol]);
					newCfg->nBonds[iPol]++;
					return newCfg;
				}
			}
		}
		if( occ[0] && occ[1] ){
			for(int iPol=0; iPol<cfg->nPol; iPol++){
				if(cfg->pos[iPol][0] == mt->move[0] && cfg->nBonds[iPol] != 1 && cfg->pos[iPol][1] == mt->move[1]){
					
					newCfg = NewFromConfig(cfg);
					for(int iBond=newCfg->nBonds[iPol]; iBond>0; iBond--)
						newCfg->pos[iPol][iBond-1] = cfg->pos[iPol][iBond];
					for(int iBond=newCfg->nBonds[iPol]-1; iBond>0; iBond--)
						newCfg->bonds[iPol][iBond-1] = cfg->bonds[iPol][iBond];
					newCfg->nBonds[iPol]--;
// 					printf("---------------------\n");
// 					PrintConfig(cfg);
// 					PrintConfig(newCfg);
					return newCfg;
				}
				else if(cfg->pos[iPol][cfg->nBonds[iPol]] == mt->move[0] && cfg->nBonds[iPol] != 1 && cfg->pos[iPol][cfg->nBonds[iPol]-1] == mt->move[1]){
					newCfg = NewFromConfig(cfg);
					newCfg->nBonds[iPol]--;
					return newCfg;
				}
			}
		}
	}
	return NULL;
}

void PrintConfig(Config* cfg){
	if(cfg->nPol == 0){
		printf("< empty config >\n\n");
	}
	else {
		for(int iPol=0; iPol<cfg->nPol; iPol++){
			for(int iBond=0; iBond<cfg->nBonds[iPol]; iBond++){
				printf("%i -> [%x] -> ", cfg->pos[iPol][iBond], cfg->bonds[iPol][iBond]);
			}
			printf("%i\n", cfg->pos[iPol][cfg->nBonds[iPol]]);
		}
		printf("\n");
	}
}

void PrintTree(Node* node, int depth){
	if(!node) return;
	if(node->topo){
// 		PrintConfig(node->cfg);
	}
// 	for(int i=0; i<depth; i++) printf("-");
// 	printf("[%i]\n", node->val);
	PrintTree(node->child, depth+1); 
	PrintTree(node->next, depth);
}

void UpdateConfig(Config* cfg, int val){
	if(val==END){
		if(cfg->nBonds[cfg->nPol]>0)
			cfg->nPol++;
		return;
	}
	int iPol = cfg->nPol;
	int iPos = cfg->nBonds[iPol];
	cfg->pos[iPol][iPos+1] = val;
	if(iPos>=0){
		int prevPos = cfg->pos[iPol][iPos];
		int fBond=0;
		for(int bond=1; bond<0xf; bond++){
			if(!IsValid(bond)) continue;
			if(AddUnitToCoor(bond, prevPos, ee.lattice) == val){
				cfg->bonds[iPol][iPos] = bond;
				fBond=1;
				break;
			}
		}
		if(!fBond) { printf("Error finding bond: %i => %i...\n", prevPos, val); exit(192); }
	}
	cfg->nBonds[iPol]++;
}

///Two kinds of mutations: Add a new two-step, or extend an existing 

int MutateFindConfig(ExactEnum* ee, Topo* topo, Node* node, Config curConfig){
	Config* newCfg;
	int ret=0;
	if(!node) return 0;
	Config oldConfig = curConfig;
	UpdateConfig(&curConfig, node->val);
	if(node->topo){
		for(int i=0; i<ee->nMutator; i++){
// 			if(topo->mutTopo[i]) continue;
			newCfg = MutateConfig(&curConfig, ee->allMutator+i);
			if(newCfg==NULL) continue;
// 			CheckIntegrity(newCfg);
			Topo* newTopo=NULL;
			int nConfig = FindInternalCfg(ee, newCfg, &newTopo);
			
			if(!topo->supTopo->mutTopo[i]){
				topo->supTopo->mutTopo[i] = newTopo;
			}
			else if(topo->supTopo->mutTopo[i]->supTopo != newTopo->supTopo){
				MergeSuperTopo(topo->supTopo->mutTopo[i]->supTopo, newTopo->supTopo);
			}
			DeleteConfig(newCfg);
		}
		ret++;
	}
	ret += MutateFindConfig(ee, topo, node->next, oldConfig);
	ret += MutateFindConfig(ee, topo, node->child, curConfig);
	return ret;
}
/*
void MutateTopo(ExactEnum* ee, Topo* topo){
	Config* startCfg = NewConfig();
	MutateFindConfig(ee, topo, topo->allConfigs, *startCfg);
	DeleteConfig(startCfg);
	topo->mutated=1;
	
	for(int i=0; i<ee->nMutator; i++){
		if(topo->mutTopo[i] && !topo->mutTopo[i]->mutated){
			MutateTopo(ee, topo->mutTopo[i]);
		}
	}
}*/

void DirectMutateTopo(ExactEnum* ee){
	int iTopo=0;
	int nDepth=0;
	int nextDepthStart=1;
	int maxDepth=ee->maxDepth;
	
	while(iTopo<ee->nTopo){

		Config* startCfg = NewConfig();
		int ret= MutateFindConfig(ee, ee->allTopo+iTopo, ee->allTopo[iTopo].allConfigs, *startCfg);
// 		printf("++++++++++++++++++++++ %i configs found ++++++++++++++++\n", ret); 
// 		PrintConfig(startCfg)
		DeleteConfig(startCfg);
		printf("\b \b\b \b\b \b\b \rnNodes = %i: %.2lf GB, nTopo=%i/%i, depth=%i/%i, next=%i, %.1lf%% done, %i left", ee->nNodes, (double)(ee->nNodes*sizeof(Node))*1e-9, iTopo+1, ee->nTopo, nDepth, maxDepth, nextDepthStart, 100*(iTopo+1)/(double)ee->nTopo, ee->nTopo-(iTopo+1)); 
		fflush(NULL);
// 		free(startCfg);
		DeleteNodeTree(ee->allTopo[iTopo].allConfigs);
// 		if(nDepth==1) exit(0);
		iTopo++;
// 		if(iTopo == nextDepthStart){
// 			nDepth++;
// 			nextDepthStart=ee->nTopo;
// 			if(nDepth > maxDepth){
// 				ee->maxTopo = iTopo;
// 				break;
// 			}
// 		}

	}
	printf("\n");
	ee->nTopo = iTopo;
}

void WriteSuperTopo(ExactEnum* ee, char* file){
// 	int** sTopoMap;
	
// 	sTopoMap = malloc(sizeof(int*)*MAX_TOPO));
// 	for(int i=0; i<MAX_TOPO; i++){
// 		sTopoMap[i] = malloc(sizeof(int)*NMUTATOR);
// 		for(int j=0; j<NMUTATOR; j++)
// 			sTopoMap[i][j] = -1;
// 	}
	
	SuperTopo** stack = malloc(sizeof(SuperTopo*)*ee->nTopo);
	int id=0;
	int stackSize=1;
	stack[0] = ee->allTopo[0].supTopo;
	stack[0]->id = 0;
	stack[0]->level = 0;
	
	while(id<stackSize && stack[id]->level < ee->maxDepth){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(!stack[id]->mutTopo[iMut]) continue;
// 			if(stack[id]->mutTopo[iMut]-ee->allTopo >= ee->maxTopo) continue;
			SuperTopo* sNew = stack[id]->mutTopo[iMut]->supTopo;
			if(sNew->id == -1 && stack[id]->level < ee->maxDepth){
				stack[stackSize] = sNew;
				sNew->id = stackSize++;
				sNew->level = stack[id]->level+1;
			}
		}
		id++;
	}
	
	FILE* pFile = fopen(file, "w");
	
	for(int iSupTopo=0; iSupTopo<stackSize; iSupTopo++){
		int srcId = stack[iSupTopo]->id;
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(stack[iSupTopo]->mutTopo[iMut]){
				int destId = stack[iSupTopo]->mutTopo[iMut]->supTopo->id;
				if(destId >= 0)
					fprintf(pFile, "%i %i %i\n", srcId, iMut, destId);
			}
		}
	}
	fclose(pFile);
	free(stack);
}

/*
void WriteTopo(ExactEnum* ee, char* file){
	FILE* pFile = fopen(file, "w");
	
	for(int iTopo=0; iTopo<ee->nTopo; iTopo++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(ee->allTopo[iTopo].mutTopo[iMut] && ee->allTopo[iTopo].mutTopo[iMut]-ee->allTopo < ee->nTopo){
				fprintf(pFile, "%i %i %li\n", iTopo, iMut, ee->allTopo[iTopo].mutTopo[iMut]-ee->allTopo);
			}
		}
	}
	fclose(pFile);
}
*/
