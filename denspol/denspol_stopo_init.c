#include "denspol.h"

SuperTopo* NewSuperTopo(Topo* topo){
	SuperTopo* sTopo = malloc(sizeof(SuperTopo));
	sTopo->topo = topo;
	sTopo->id = NON_EXISTING;
	sTopo->permBond = NON_EXISTING;
	sTopo->keyInserted = 0;
	sTopo->topoSameBonds = NULL;
	for(int i=0; i<NMUTATOR; i++) sTopo->mutTopo[i] = NULL;
	return sTopo;
}

void ClearTopo(Topo* topo){
	topo->nBonds=0;
	for(int i=0; i<NMUTATOR; i++) topo->mutTopo[i]=NULL;
	topo->next = NULL;
	topo->supTopo = NewSuperTopo(topo);
	topo->set=0;
}

Topo* NewTopo(){
	Topo* topo = malloc(sizeof(Topo));
	topo->nBonds=0;
	for(int i=0; i<NMUTATOR; i++) topo->mutTopo[i]=NULL;
	topo->next = NULL;
	topo->set = 0;
	topo->supTopo = NewSuperTopo(topo);
	return topo;
}

void CopyTopo(Topo* src, Topo* dst){
	
	dst->nBonds=src->nBonds;
	for(int i=0; i<src->nBonds; i++){
		for(int k=0; k<2; k++){
			dst->bonds[i][k] = src->bonds[i][k];
		}
		dst->bondFree[i] = src->bondFree[i];
	}
	
}

void MergeTopo(Topo* topoA, Topo* topoB){
	if(!(topoA && topoB)) return;
	
	SuperTopo *sTopoA = topoA->supTopo;
	SuperTopo *sTopoB = topoB->supTopo;
	
	if(sTopoA == sTopoB) return;
	
	Topo* curTopo=sTopoA->topo;
	
	while(curTopo->next){
		curTopo->supTopo = sTopoB;
		curTopo = curTopo->next;
	}
	
	curTopo->supTopo = sTopoB;
	curTopo->next = sTopoB->topo;
	
	sTopoB->topo = sTopoA->topo;
	
	free(sTopoA);
	nSupTopo--;
}

void PrintTopo(Topo* topo){
	if(!topo->nBonds) printf("[empty]\n");
	else{
		for(int iBond=0; iBond<topo->nBonds; iBond++){
			printf("["); PrintUnit(topo->bonds[iBond][0]); printf(" => "); PrintUnit(topo->bonds[iBond][1]); printf("]\n");
		}
	}
	printf("\n");
}

void PrintSTopo(SuperTopo* sTopo){
	
	int minBonds=999;
	Topo* topoMin = NULL;
	
	Topo* topo = sTopo->topo;
	
	while(topo){
		if(topo->nBonds <minBonds){
			topoMin = topo;
			minBonds = topo->nBonds;
		}
		topo = topo->next;
	}
	PrintTopo(topoMin);
}

void MutateTopo(Topo* topo, int iMut, Topo* destTopo, LookupTables* lt){
	
	int *newBond = lt->revMutTable[iMut/2];
	CopyTopo(topo, destTopo); 
	
	if(iMut%2 == 1){
		int exists=0;
		for(int iBond=0; iBond<topo->nBonds; iBond++){
			if( destTopo->bonds[iBond][0] == newBond[0] && destTopo->bonds[iBond][1] == newBond[1] ){
				destTopo->nBonds--;
				for(int k=0; k<2; k++)
					destTopo->bonds[iBond][k] = destTopo->bonds[destTopo->nBonds][k];
				exists = 1;
			}
		}
		
		if(!exists){
			for(int k=0; k<2; k++)
				destTopo->bonds[destTopo->nBonds][k] = newBond[k];
			destTopo->nBonds++;
		}
	}
	else{
		int inserted=0;
		for(int iBond=0; iBond<topo->nBonds && ! inserted; iBond++){
			for(int k=0; k<2 && !inserted; k++){
				if(destTopo->bonds[iBond][k] == (newBond[0]^(k*0xf)) ){
					destTopo->bonds[iBond][k] = newBond[1]^((0x1^k)*0xf);
					inserted = 1;
				}
			}
		}
		if(!inserted){
			printf("Error: mutation failed to insert bond\n");
			PrintTopo(topo); 
			PrintTopo(destTopo);
			PrintUnit(newBond[0]); printf(" --> "); PrintUnit(newBond[1]); printf("\n");
			exit(0);
		}
	}
}

KeyNode* NewKey(int val){
	KeyNode* kn = malloc(sizeof(KeyNode));
	
	kn->val = val;
	kn->next = NULL;
	kn->child = NULL;
	kn->sTopo = NULL;
	return kn;
}

KeyNode* InsertSTKey(KeyNode* kn, int* key, int nKey, SuperTopo* sTopo){
	KeyNode* curKey=NULL;
	
	if(!kn)
		kn = NewKey(key[0]);
	
	curKey = kn;
	while(curKey && curKey->val != key[0])
		curKey = curKey->next;
	
	if(!curKey){
		curKey = NewKey(key[0]);
		curKey->next = kn;
		kn = curKey;
	}
	
	if(nKey == 1){
		if(curKey->sTopo){
			sTopo->topoSameBonds = curKey->sTopo->topoSameBonds;
			curKey->sTopo->topoSameBonds = sTopo;
// 			printf("Double key found\n");
// 			exit(0);
		}
		else{
			curKey->sTopo = sTopo;
			sTopo->topoSameBonds = sTopo;
		}
	}
	else
		curKey->child = InsertSTKey(curKey->child, key+1, nKey-1, sTopo);
	
	return kn;
}

void SortKey(int key[12], int nKey){
	int newKey[12];
	
	for(int iKey=0; iKey<nKey/2; iKey++){
		int minBond=999;
		int iMin=-1;
		for(int jKey=0; jKey<nKey/2; jKey++){
			if(key[jKey*2]<minBond){
				minBond=key[jKey*2];
				iMin=jKey;
			}
		}
		newKey[2*iKey  ] = key[2*iMin];
		newKey[2*iKey+1] = key[2*iMin+1];
		key[2*iMin  ] = 999;
		key[2*iMin+1] = 999;
	}
	for(int i=0; i<nKey; i++) key[i] = newKey[i];
// 	exit(0);
}

void UpdateTopoMutations(Topo* topo, LookupTables* lt){
	
	int permBond=0;
	int key[12];
	int nKey=0;
	for(int iBond=0; iBond<topo->nBonds; iBond++){
		int mut = lt->mutIdTableDouble[topo->bonds[iBond][0]][topo->bonds[iBond][1]];
		if(mut>=0 && topo->mutTopo[2*mut+1]) continue;
		key[nKey  ] = topo->bonds[iBond][0];
		key[nKey+1] = topo->bonds[iBond][1];
		nKey+=2;
		permBond |= 1 << topo->bonds[iBond][0];
		permBond |= 1 << (topo->bonds[iBond][1]^0xf);
	}
	
	if(topo->supTopo->permBond >= 0 && topo->supTopo->permBond != permBond){
		printf("Making a mistake while figuring out the perma bonds\n");
		printf("%x vs %x\n", topo->supTopo->permBond, permBond);
		PrintTopo(topo);
		exit(0);
	}
	
	if(!topo->supTopo->keyInserted){
		if(nKey == 0){
			key[nKey++] = NON_EXISTING;
			key[nKey++] = NON_EXISTING;
		}
		SortKey(key, nKey);
// 		for(int i=0; i<nKey; i++) printf("%i ", key[i]);
// 		printf("\n");
		lt->sameTopoTree = InsertSTKey(lt->sameTopoTree, key, nKey, topo->supTopo);
		topo->supTopo->keyInserted = 1;
	}
	
	topo->supTopo->permBond = permBond;
	
	for(int iMut=0; iMut<NMUTATOR; iMut++){
		if(!topo->mutTopo[iMut]) continue;
		
		Topo* destTopo = topo->mutTopo[iMut];
		if(!destTopo->set){
			MutateTopo(topo, iMut, destTopo, lt);
		}
		
		if(iMut%2 == 1) /// if it only adds/removes a link, it does not alter the super topo.
			continue;
		
		int start = lt->revMutTable[iMut/2][0];
		int end = lt->revMutTable[iMut/2][1];
		
		int foundBond=0;
		for(int iBond=0; iBond<topo->nBonds && !foundBond; iBond++){
			int oldMut = lt->mutIdTableDouble[topo->bonds[iBond][0]][topo->bonds[iBond][1]];
			for(int k=0; k<2 && !foundBond; k++){
				if((topo->bonds[iBond][k]^(k*0xf)) == start){
					int mutation=-1;
					if(oldMut>=0 && topo->mutTopo[2*oldMut+1]){
						if(k==0)
							mutation = lt->mutIdTableTriple[end^0xf][start^0xf][topo->bonds[iBond][1]];
						if(k==1)
							mutation = lt->mutIdTableTriple[topo->bonds[iBond][0]][start^0xf][end];
					}
					else{
						mutation = iMut/2; 
					}
					if(mutation == NON_EXISTING){
						printf("Uh oh! Not a valid mutation!\n");
						PrintUnit(topo->bonds[iBond][0]); printf(" -> "); PrintUnit(topo->bonds[iBond][1]); printf(" + "); PrintUnit(start); printf(" -> "); PrintUnit(end);
						printf("\nk=%i\n", k);
						exit(0);
					}
					else if (mutation != SAME_TOPO){
						if(topo->supTopo->mutTopo[mutation] && topo->supTopo->mutTopo[mutation]->supTopo != destTopo->supTopo){
							printf("Overwriting mutation!\n");
							exit(192);
						}
						topo->supTopo->mutTopo[mutation] = destTopo;
					}
					foundBond=1;
				}
			}
		}
	}
}

void PrintPermBond(int permBond){
	printf("(");
	for(int iB=0; iB<16; iB++){
		if(permBond & (1<<iB)){
			PrintUnit(iB); printf(" ");
		}
	}
	printf(")");
}

void PrintStraightTopo(LookupTables* lt){
	FILE* pFile = fopen("straight_topo.dat", "w");
	for(int i=0; i<lt->nTopoComp; i++){
		int permBond = lt->topComp[i].permBond;
		int bonds[2], nBonds=0;
		for(int iB=0; iB<16; iB++){
			if(permBond & (1<<iB)){
				if(nBonds>=2){
					nBonds=3; break;
				}
				else
					bonds[nBonds++] = iB;
			}
		}
		
		if(nBonds==2 && bonds[0] == 15-bonds[1]){
// 			fprintf(pFile, "%x %x %i\n", bonds[0], bonds[1], i);
			
			int found=0;
			for(int iMut=0; iMut<NMUTATOR && !found; iMut++){
				int topo1 = lt->topComp[i].mutators[iMut];
				if(topo1 < 0 ) continue;
				
				for(int jMut=0; jMut<NMUTATOR && !found; jMut++){
					int topo2 = lt->topComp[topo1].mutators[jMut];
					if(topo2 != 0) continue;
					
// 					PrintPermBond(lt->topComp[topo1].permBond);
// 					PrintUnit(bonds[0]); printf(" "); PrintUnit(bonds[1]); PrintPermBond(lt->topComp[i].permBond); printf(" => ");
// 					for(int k=0; k<2; k++){
// 						PrintUnit(lt->revMutTable[iMut][k]); printf(" ");
// 					}
// 					PrintPermBond(lt->topComp[topo1].permBond);
// 					printf("\n");
// 					for(int k=0; k<2; k++){
// 						PrintUnit(lt->revMutTable[jMut][k]); printf(" ");
// 					}
// 					printf("\n");
// 					printf("%i (%i)=> %i (%i)=> %i\n", i, iMut, topo1, jMut, topo2);
					
					int tripleMut=-1;
					for(int kMut=0; kMut<NMUTATOR; kMut++){
						if(lt->topComp[0].mutators[kMut] == topo1){
// 							for(int k=0; k<3; k++){
// 								PrintUnit(lt->revMutTableTriple[kMut][k]); printf(" ");
// 							}
							tripleMut = kMut;
						}
					}
// 					printf("\n\n");
					
					for(int k=0; k<2; k++){
						for(int l=0; l<3; l+=2){
							if(lt->revMutTableTriple[tripleMut][l] == bonds[k]){
								fprintf(pFile, "%x %i\n", bonds[k], i);
// 								printf("Correct bond = %i\n", bonds[k]);
								found=1;
							}
						}
					}
				}
			}
// 			printf("---------------------------\n");
// 			PrintUnit(bonds[0]); printf(" -> "); PrintUnit(bonds[1]); printf(" = %i\n", i); 
		}
	}
	exit(0);
}

void SuperTopoInit(LookupTables* lt){
	nSupTopo = lt->nTopo;
	Topo* allTopo = malloc(sizeof(Topo)*lt->nTopo);
	
	for(int i=0; i<lt->nTopo; i++) ClearTopo(allTopo+i);
	
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(lt->mutTopo[iTopo][iMut] == -1) continue;
			allTopo[iTopo].mutTopo[iMut] = allTopo + lt->mutTopo[iTopo][iMut];
		}
	}
	
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(iMut%2 == 1){
				MergeTopo(allTopo+iTopo, allTopo[iTopo].mutTopo[iMut]);
			}
		}
	}
	
	allTopo[0].set=1;
	for(int iTopo=0; iTopo<lt->nTopo; iTopo++){
		UpdateTopoMutations(allTopo+iTopo, lt);
	}
// 	exit(0);
// 	for(int i=0; i<NMUTATOR; i++){
// 		if(allTopo[0].supTopo->mutSupTopo[i]){
// 			printf("%i\n", i);
// 		}
// 		else{
// 			printf("N/A\n");
// 		}
// 	}
	
	SuperTopo** sTopoArray = malloc(sizeof(SuperTopo*)*nSupTopo);
	lt->topComp = malloc(sizeof(SuperTopo)*nSupTopo);
	sTopoArray[0] = allTopo[0].supTopo;
	sTopoArray[0]->id = 0;
	for(int i=0; i<nSupTopo; i++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			lt->topComp[i].mutators[iMut] = NON_EXISTING;
		}
	}
	
	int nIds=1;
	for(int i=0; i<nSupTopo; i++){
		lt->topComp[i].permBond = sTopoArray[i]->permBond;
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			if(!sTopoArray[i]->mutTopo[iMut]) continue;
			SuperTopo* newSTopo= sTopoArray[i]->mutTopo[iMut]->supTopo;
			if(newSTopo->id <0 ){
				newSTopo->id = nIds;
				sTopoArray[nIds++] = newSTopo;
			}
			lt->topComp[i].mutators[iMut] = newSTopo->id;
		}
	}
	
	lt->nTopoComp = nIds;
	
	for(int i=0; i<nSupTopo; i++){
		SuperTopo* sTopo = sTopoArray[i];
		lt->topComp[i].sameTopo = sTopo->topoSameBonds->id;
	}
	
	for(int i=0; i<nSupTopo; i++){
		for(int iMut=0; iMut<NMUTATOR; iMut++){
			int newTopo = lt->topComp[i].mutators[iMut];
			if(newTopo < 0 ) continue;
			int foundBack=0;
			for(int jMut=0; jMut<NMUTATOR && !foundBack; jMut++){
				if(lt->topComp[newTopo].mutators[jMut] == i)
					foundBack=1;
			}
			if(!foundBack){
				PrintTopo(sTopoArray[i]->topo);
				for(int iMut=0; iMut<NMUTATOR; iMut++){
					if(sTopoArray[i]->mutTopo[iMut]){
						PrintSTopo(sTopoArray[i]->mutTopo[iMut]->supTopo);
					}
				}
				printf("==============================\n");
				PrintTopo(sTopoArray[newTopo]->topo);
				for(int iMut=0; iMut<NMUTATOR; iMut++){
					if(sTopoArray[newTopo]->mutTopo[iMut]){
						PrintSTopo(sTopoArray[newTopo]->mutTopo[iMut]->supTopo);
					}
				}
				printf("Error contructing topComp, no detailed balance\n");
				exit(192);
			}
		}
	}
	free(allTopo);
	free(sTopoArray);
}
