#include "ee.h"



int main(int argc, char** argv){
	EEInit();
	
	
// 	MutateFindConfig(&ee, ee.allTopo, ee.allConfigs, *NewConfig());
	DirectMutateTopo(&ee);
	WriteTopo(&ee, "ee_topo.dat");
// 	PrintTree(ee.allConfigs, 0);
	return 0;
}

void EEInit(){
	ee.nConfig=0;
	ee.nNodes=0;
	ee.mt = NewMoveTable();
	ee.lattice = NewLattice(5);
	LatticeSetHole(ee.lattice, ee.mt);
	SetCornerMoves(&ee);
	GenerateMutators(&ee);
	EETopoInit(&ee);
// 	ee.recStep=0;
// 	ee.nTopo=0;
	
// 	PrintTree(ee.allConfigs, 0);
// 	PrintMoveTable(ee.mt);
// 	PrintLattice(ee.lattice);
// 	PrintCornerMoves(&ee);
// 	PrintMutators(&ee);
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
			ee->cornerMoves[ee->nCornerMoves++][2] = unitC;
		}
	}
}

void EETopoInit(ExactEnum* ee){
	Topo* dstTopo=NULL;
// 	ee->allTopo = NewTopo();
	Config* startCfg = NewConfig();
	startCfg->nPol=0;
	startCfg->topo=ee->allTopo;
	ee->nTopo=0;
	int inserted;
	PrintConfig(startCfg);
	ee->allTopo = malloc(sizeof(Topo)*MAX_TOPO);
	ee->allConfigs = InsertConfig(NULL, 0, startCfg, &inserted, &dstTopo);
	dstTopo->allConfigs = InsertConfig(NULL, 0, startCfg, &inserted, &dstTopo);
	DeleteConfig(startCfg);
	if(!ee->allConfigs) printf("???\n");
}

// Topo* NewTopo(){
// 	Topo* topo = malloc(sizeof(Topo));
// 	for(int i=0; i<48; i++) topo->mutTopo[i]=NULL;
// // 	topo->nMutations=0;
// 	topo->allConfigs=NULL;
// 	topo->mutated=0;
// 	return topo;
// }

Topo* NewTopo(){
	if(ee.nTopo>=MAX_TOPO){
		printf("Error: Not enough memory allocated for topo\n");
		exit(192);
	}
	Topo* topo = ee.allTopo+ee.nTopo;
	for(int i=0; i<48; i++) topo->mutTopo[i]=NULL;
	topo->allConfigs=NULL;
	topo->mutated=0;
	ee.nTopo++;
	return topo;
}

void GenerateMutators(ExactEnum* ee){
	///First generate double mutators
	
	ee->nMutator = 0;
	
	for(int i=0; i<ee->nCornerMoves; i++){
		for(int k=0; k<2; k++)
			ee->allMutator[ee->nMutator].move[k] = ee->cornerMoves[i][k];
		ee->allMutator[ee->nMutator++].bond = ee->cornerMoves[i][2];
// 		ee->allMutator[ee->nMutator++].initial = 0;
	}
	
// 	for(int i=0; i<ee->nCornerMoves; i++){
// 		for(int j=0; j<ee->nCornerMoves; j++){
// 			if(ee->cornerMoves[i][1] != ee->cornerMoves[j][0]) continue;
// 			if(ee->cornerMoves[i][0] == ee->cornerMoves[j][1]) continue;
// 			int loop=0;
// 			for(int k=0; k<ee->nCornerMoves; k++){
// 				if(ee->cornerMoves[i][0] == ee->cornerMoves[k][0] && ee->cornerMoves[j][1] == ee->cornerMoves[k][1]){
// 					loop=1;
// 					break;
// 				}
// 			}
// 			if(loop) continue;
// 			ee->allMutator[ee->nMutator  ].move[0] = ee->cornerMoves[i][0];
// 			ee->allMutator[ee->nMutator  ].move[1] = ee->cornerMoves[i][1];
// 			ee->allMutator[ee->nMutator  ].move[2] = ee->cornerMoves[j][1];
// 			ee->allMutator[ee->nMutator++].initial = 1;
// 		}
// 	}
}

void PrintMutators(ExactEnum* ee){
	printf("Number of mutators: %i\n", ee->nMutator);
	
	for(int i=0; i<ee->nMutator; i++){
		Mutator* mut = ee->allMutator+i;
		printf("%i -> %i [0x%x]", mut->move[0], mut->move[1], mut->bond);
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
	if(CheckKey(cfg->key, cfg->nKey)){
		PrintConfig(cfg);
		printf("Key:\n\n");
		for(int i=0; i<cfg->nKey; i++)
			printf("%i ", cfg->key[i]);
		printf("\n\n");
		exit(0);
	}
}

Node* InsertConfig(Node* node, int iKey, Config* cfg, int* inserted, Topo** destTopo){
	if(!iKey){
		SetConfigKey(cfg);
		*inserted=0;
// 		*destTopo=NULL;
	}
	
	int val = cfg->key[iKey];
	if(!node){ node = NewNode(val); }
	
	if(node->val == val){
		if(iKey == cfg->nKey-1){
			if(! node->topo){
// 				node->cfg = cfg;
				if(!*destTopo){
// 					printf("Created new Topo\n");
					node->topo = NewTopo();
// 					ee.nTopo++;
				}
				else
					node->topo = *destTopo;
				*inserted=1;
// 				ee.nConfig++;
			}
// 			else{
// 				free(cfg);
// 			}
			*destTopo = node->topo;
// 			printf("destTopo = %lx\n", *destTopo);
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

int FindInternalCfg(ExactEnum* ee, Config* cfg, Topo** topo){
	MoveTable* mt = ee->mt;
// 	Topo* topo = ee->allTopo+ee->nTopo;
// 	Config* newCfg = malloc(sizeof(Config));
	int totInsert=0;
	int dummy;
// 	CheckIntegrity(cfg);
// 	Topo* destTopo = cfg->topo;
// 	printf("Start: \n");
// 	PrintConfig(cfg);
	ee->allConfigs = InsertConfig(ee->allConfigs, 0, cfg, &totInsert, topo);
	if(!totInsert){ /*printf("config already exists\n");*/ return 0; }
	if(!(*topo)){printf("ooops\n"); exit(192);}
	(*topo)->allConfigs = InsertConfig((*topo)->allConfigs, 0, cfg, &dummy, topo);
// 	printf("Insert complete\n");
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
// 				printf("%i -> %i + %i [pos=%i]\n", bond, newBond1, newBond2, newPos);
// 				PrintLattice(ee->lattice);
				Config* newCfg = NewFromConfig(cfg);
				InsertValue(newPos, iBond+1, cfg->nBonds[iPol]+1, newCfg->pos[iPol]);
				InsertValue(newBond2, iBond+1, cfg->nBonds[iPol], newCfg->bonds[iPol]);
				newCfg->nBonds[iPol]++;
				newCfg->bonds[iPol][iBond] = newBond1;
				
// 				SetMidLattice(newPos, 1, ee->lattice);
				int newInsert = FindInternalCfg(ee, newCfg, topo);
				DeleteConfig(newCfg);
// 				free(newCfg);
// 				printf("Finish forw\n");
// 				SetMidLattice(newPos, 0, ee->lattice);
				totInsert+=newInsert;
			}
			/// Try backward moves.
			if((iBond>=cfg->nBonds[iPol]-1)) continue;
			int newBond = mt->backMoves[bond][cfg->bonds[iPol][iBond+1]];
			if(!newBond) continue;
// 			int oldPos = cfg->pos[iPol][iBond];
			
			Config* newCfg = NewFromConfig(cfg);
			for(int i=iBond+1; i<cfg->nBonds[iPol]-1; i++){
				newCfg->bonds[iPol][i] = cfg->bonds[iPol][i+1];
				newCfg->pos[iPol][i] = cfg->pos[iPol][i+1];
			}
			newCfg->pos[iPol][cfg->nBonds[iPol]-1] = cfg->pos[iPol][cfg->nBonds[iPol]];
			newCfg->bonds[iPol][iBond] = newBond;
			newCfg->nBonds[iPol]--;
			
// 			SetMidLattice(oldPos, 0, ee->lattice);
			int newInsert = FindInternalCfg(ee, newCfg, topo);
			DeleteConfig(newCfg);
// 			free(newCfg);
// 			SetMidLattice(oldPos, 1, ee->lattice);
// 			if(newInsert) newCfg = malloc(sizeof(Config));
			totInsert+=newInsert;
// 			printf("Finish back\n");
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

// Config* MutateConfig(Config* cfg, Mutator* mt){
// 	Config* newCfg = NULL;
// 	if(mt->initial){
// 		int available=1;
// 		for(int iPol=0; iPol<cfg->nPol && available; iPol++){
// 			for(int iMono=0; iMono<=cfg->nBonds[iPol] && available; iMono++){
// 				int pos = cfg->pos[pos][iMono];
// 				for(int k=0; k<3 && available; k++){
// 					if(pos == mt->move[k])
// 						available=0;
// 				}
// 			}
// 		}
// 		if(available){
// 			newCfg = malloc(sizeof(Config));
// 			CopyConfig(cfg, newCfg);
// 			for(int k=0; k<3; k++){
// 				newCfg->bonds[newCfg->nPol][k] = mt->move[k];
// 			}
// 		}
// 	}
// 	
// }

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

// void InsertPos(int pos, int iPos, int iPol, Config* cfg){
// 	
// }

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
// 	CheckIntegrity(cfg);
	int occ[2];
	occ[0] = ConfigOccupied(mt->move[0], cfg);
	occ[1] = ConfigOccupied(mt->move[1], cfg);
// 	printf("occ = {%i,%i}\n", occ[0], occ[1]);
	
	if( ! (occ[0] || occ[1]) ){
		newCfg = NewFromConfig(cfg);
		newCfg->pos[newCfg->nPol][0] = mt->move[0];
		newCfg->pos[newCfg->nPol][1] = mt->move[1];
		newCfg->bonds[newCfg->nPol][0] = mt->bond;
		newCfg->nBonds[newCfg->nPol++] = 1;
// 		printf("hi!\n");
// 		PrintConfig(cfg);
// 		PrintConfig(newCfg);
// 		printf("-------------\n");
		return newCfg;
	}
	else if ( occ[0] && ! occ[1] ){
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
	return NULL;
}

void PrintConfig(Config* cfg){
	if(cfg->nPol == 0){
		printf("< empty config >\n\n");
	}
	else {
		for(int iPol=0; iPol<cfg->nPol; iPol++){
			for(int iBond=0; iBond<cfg->nBonds[iPol]; iBond++){
				printf("%i -> [%i] -> ", cfg->pos[iPol][iBond], cfg->bonds[iPol][iBond]);
			}
			printf("%i\n", cfg->pos[iPol][cfg->nBonds[iPol]]);
		}
		printf("\n");
	}
}

// void PrintTree(Node* node, int depth){
// 	if(!node) return;
// 	for(int i=0; i<depth; i++) printf("-");
// 	printf("[%i]\n", node->val);
// 	PrintTree(node->child, depth+1); 
// 	PrintTree(node->next, depth);
// }

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

void MutateFindConfig(ExactEnum* ee, Topo* topo, Node* node, Config curConfig){
	Config* newCfg;
	if(!node) return;
// 	printf("Value=%i\n", node->val);
	Config oldConfig = curConfig;
	UpdateConfig(&curConfig, node->val);
// 	if(ee->nConfig >= 10000){
// 		printf("\n");
// 		PrintTree(ee->allConfigs, 0);
// 		exit(0);
// 	}
	if(node->topo){
// 		PrintConfig(&curConfig);
// 		printf("1\n");
// 		CheckIntegrity(&curConfig);
// 		PrintTree(ee->allConfigs, 0);
// 		printf("nMutator = %i\n", ee->nMutator);
	// 	if(node->cfg){
		for(int i=0; i<ee->nMutator; i++){
			if(topo->mutTopo[i]) continue;
			newCfg = MutateConfig(&curConfig, ee->allMutator+i);
			if(newCfg==NULL) continue;
	// 			int inserted;
	// 			ee->allConfigs = InsertConfig(ee->allConfigs, 0, cfg, &inserted);
	// 			if(!inserted) continue;
// 			PrintConfig(newCfg);
// 			printf("2\n");
// 			CheckIntegrity(newCfg);
// 			printf("3\n");
			Topo* newTopo=NULL;
			FindInternalCfg(ee, newCfg, &newTopo);
			DeleteConfig(newCfg);
// 			free(newCfg);
// 				PrintTree(newTopo->allConfigs, 0);
// 				printf("Adding new topo: %lx\n", newTopo->allConfigs);
// 				PrintTree(newTopo->allConfigs, 0);
// 				printf("++++++++++++++++++\n");
// 				Config* startCfg = NewConfig();
// 				MutateFindConfig(ee, newTopo, newTopo->allConfigs, *startCfg);
// 				free(startCfg);
// 			}
// 			else{
			topo->mutTopo[i] = newTopo;
// 			}
	// 		PrintTree(ee->allConfigs, 0);
	// 		exit(0);
		}
	}
// 	Node* curNode = node->next;
// 	while(curNode){
// 		MutateFindConfig(ee, topo, curNode, oldConfig);
// 		curNode = curNode->next;
// 	}
	MutateFindConfig(ee, topo, node->next, oldConfig);
	MutateFindConfig(ee, topo, node->child, curConfig);
// 	}
}

void MutateTopo(ExactEnum* ee, Topo* topo){
	Config* startCfg = NewConfig();
	MutateFindConfig(ee, topo, topo->allConfigs, *startCfg);
	DeleteConfig(startCfg);
// 	free(startCfg);
	topo->mutated=1;
	
	for(int i=0; i<ee->nMutator; i++){
		if(topo->mutTopo[i] && !topo->mutTopo[i]->mutated){
			MutateTopo(ee, topo->mutTopo[i]);
		}
	}
}

void DirectMutateTopo(ExactEnum* ee){
	int iTopo=0;
	
	while(iTopo<ee->nTopo){
		Config* startCfg = NewConfig();
		MutateFindConfig(ee, ee->allTopo+iTopo, ee->allTopo[iTopo].allConfigs, *startCfg);
		DeleteConfig(startCfg);
		printf("\b \b\b \b\b \b\b \rnNodes = %i: %.2lf GB, nConfig=%i: %.2lf GB, nTopo=%i: %.2lf GB, %i/%i=%.1lf%% done, %i left", ee->nNodes, (double)(ee->nNodes*sizeof(Node))*1e-9, ee->nConfig, ee->nConfig*sizeof(Config)*1e-9, ee->nTopo, ee->nTopo*sizeof(Topo)*1e-9, iTopo, ee->nTopo, 100*iTopo/(double)ee->nTopo, ee->nTopo-iTopo); 
		fflush(NULL);
// 		free(startCfg);
		DeleteNodeTree(ee->allTopo[iTopo].allConfigs);
		iTopo++;
	}
}

void WriteTopo(ExactEnum* ee, char* file){
	FILE* pFile = fopen(file, "w");
	
	for(int iTopo=0; iTopo<ee->nTopo; iTopo++){
		for(int iMut=0; iMut<48; iMut++){
			if(ee->allTopo[iTopo].mutTopo[iMut]){
				fprintf(pFile, "%i %i %li\n", iTopo, iMut, ee->allTopo[iTopo].mutTopo[iMut]-ee->allTopo);
			}
		}
	}
	fclose(pFile);
}

