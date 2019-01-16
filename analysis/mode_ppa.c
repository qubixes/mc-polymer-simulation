typedef struct BinHeapElem{
	int mono[2];
	int gap;
	int dist;
	int priority;
	int heapId;
	struct BinHeapElem* next[2];
}BinHeapElem;

BinHeapElem* NewBinHeapElem(int gapSq, int distSq, int iMono, int jMono){
	BinHeapElem* bh = malloc(sizeof(BinHeapElem));
	bh->mono[0] = iMono;
	bh->mono[1] = jMono;
	bh->gapSq = gapSq;
	bh->distSq = distSq;
	for(int k=0; k<2; k++) bh->next[k]=NULL;
	return bh;
}
	

void ComputeRingPPA(SimProperties* sp, PolyConfig* pcfg){
	///First make the list of monomers without stored length:
	
	int nRealBond=0;
	int** xyz = pcfg->xyzRealBond;
	for(int iMono=0; iMono<pcfg->polSize; iMono++){
		int jMono = (iMono+1)%pcfg->polSize;
		if(pcfg->x[iMono] != pcfg->x[jMono] || pcfg->y[iMono] != pcfg->y[jMono] || pcfg->z[iMono] != pcfg->z[jMono]){
			xyz[nRealBond][0] = pcfg->x[iMono];
			xyz[nRealBond][1] = pcfg->y[iMono];
			xyz[nRealBond][2] = pcfg->z[iMono];
			nRealBond++;
		}
	}
	
	BinHeapElem** heap = pcfg->binHeap;
	BinHeapElem** heapStart = pcfg->heapStart;
	int nHeap=0;
	
	for(int i=0; i<nRealBond; i++){
		heapStart[i] = NULL;
	}
	for(int i=0; i<nRealBond; i++){
		int j= (i+1)%pcfg->polSize;
		heap[i] = NewBinHeapElem(1,1,i,j);
		heap[i]->next[0] = heapStart[i];
		heap[i]->next[1] = heapStart[j];
		heapStart[i] = heap[i];
		heapStart[j] = heap[j];
	}
	nHeap = nRealBond;
	
	int gap=1;
	
	while(1){
		BinHeapElem* jump = HeapExtract(heap);
		gap = jump->gap;
		///Merge left
		for(int dir=0; dir<2; dir++){
			pElem = heapStart[jump->mono[dir]];
			while(pElem){
				if(pElem->mono[dir^0x1] == jump->mono[dir] && pElem->gap <= gap){
					int newMono[2];
					newMono[0] = dir?jump->mono[0]:pElem->mono[0];
					newMono[1] = dir?pElem->mono[1]:pElem->mono[1];
					pNewElem = FindInList(newMono[0], newMono[1]);
					if(!pNewElem){
						double newDist = DRXYZ(xyz[newMono[0]], xyz[newMono[1]]);
						pNewElem = NewBinHeapElem(newDist, newDist, newMono[0], newMono[1]);
						InsertInList(pNewElem, heapStart);
					}
					newDist = MAX(MAX(pElem->gap, jump->gap), pNewElem->dist);
					if(newDist < pNewElem->gap){
						if(pNewElem->heapId >= 0)
							HeapUpdate(heap, pNewElem->heapId);
						else
							HeapInsert(heap, pNewElem);
					}
				}
				pElem = (pElem->mono[dir^0x1] == jump->mono[dir])?pElem->next[dir^0x1]:pElem->next[dir];
			}
		}
	}
	pcfg->gap = gap;
}