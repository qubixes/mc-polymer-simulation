#include "efpol.h"

uint LargeToSmallUnit(uint uLarge, Constants* con){
	uint uSmall;
	
	uSmall = (uLarge&0x1) | ((uLarge>>(con->BL))&0x2) | ((uLarge>>(2*con->BL))&0x4) | ((uLarge>>(3*con->BL))&0x8);
	return uSmall;
}

uint SmallToLargeUnit(uint uSmall, Constants* con){
	uint uLarge;
	
	uLarge = (uSmall&0x1) | ((uSmall&0x2)<<(con->BL)) | ((uSmall&0x4)<<(2*con->BL)) | ((uSmall&0x8)<<(3*con->BL));
	return uLarge;
}

int TUVtoCoor(int t, int u, int v, Constants* con){
	return (t+u*con->ALLOC_T+v*con->ALLOC_T*con->ALLOC_U);
}

uint ValidateAddUnitVectors(uint a, uint b, uint* c){
	uint r, valid;
	if((a|b) != 0xf && (a&b))
		return 0;
	r = (((a|b)==0xf)?(a&b):(a|b));
	valid = IsValid(r);
	*c = r;
	return valid;
}

uint AddUnitVectors(uint a, uint b){
	return (((a|b)==0xf)?(a&b):(a|b));
}


uint AddUnitToCoor(uint addUnit, uint coor, Constants* con){
	uint addLarge = SmallToLargeUnit(addUnit, con);
// 	printf("addUnit=%x, addLarge=%x\n", addUnit, addLarge);
	int newCoor = coor+(((con->LARGE_W_UNIT&addLarge)>>(3*con->BL+3))*con->BOUND_MASK)+(addLarge&con->LARGE_TUVW);
// 	if(addUnit&0x8)
// 	printf("%x+%x+%x (%x&%x)=>%x\n", coor, (((con->LARGE_W_UNIT&addLarge)>>(3*con->BL+3))*con->BOUND_MASK), (addLarge&con->LARGE_TUVW), addLarge, con->LARGE_TUVW, newCoor);
	return (newCoor&con->BOUND_MASK);
}


uint IsValid(uint a){
	return (a!=0 && a!=0x3 && a!=0xc && a!=0xf);
}

int GetNPol(double density, int polSize, int LAT_SIZE){
	double dPol;
	
	dPol = density*LAT_SIZE/(double)polSize;
	
	return ((int)(dPol+0.5));
}


uint TCoor(unsigned int coor, Constants* con){
	return (coor&con->T_MASK);
}

uint UCoor(unsigned int coor, Constants* con){
	return ((coor&con->U_MASK)/con->ALLOC_T);
}

uint VCoor(unsigned int coor, Constants* con){
	return ((coor&con->V_MASK)/(con->ALLOC_T*con->ALLOC_U));
}

int X(int t, int u, int v){
	return (t+u-v);
}

int Y(int t, int u, int v){
	return (t-u);
}

int Z(int t, int u, int v){
	return (v);
}

int XUnit(int unit){
	return ((unit&0x1)+((unit>>1)&0x1)-((unit>>2)&0x1)-(unit>>3));
}

int YUnit(int unit){
	return ((unit&0x1)-((unit>>1)&0x1));
}

int ZUnit(int unit){
	return (((unit>>2)&0x1)-(unit>>3));
}

void SetLattice(int coor, CurState* cs){
	uint bit, word;
	
	bit = coor%32;
	word = coor/32;
	
	cs->lattice[word] |= (uint)0x1<<bit;
}

void UnsetLattice(int coor, CurState* cs){
	uint bit, word;
	
	bit = coor%32;
	word = coor/32;
	
	cs->lattice[word] &= ~((uint)1<<bit);
}

uint OccupLattice(int coor, CurState* cs){
	uint bit, word;
	
	bit = coor%32;
	word = coor/32;
	
	return (cs->lattice[word]& ((uint)1<<bit));
}

long DoMCStep(CurState* cs, long nTime){
	long nSmall = nTime*cs->polSize*((cs->nPol-1)/BLOCK_SELECT+1);
	long realStep=0;
	for(long j=0; j<nSmall; j++){
		realStep += DoStep(cs);
	}
	return realStep;
}


int DoStep(CurState* cs){
	int mono, prevMono;
	unsigned int src1, src2, dst1, dst2, move;
	uint option, randOpt;
	int newCoor=-1;
	int polStart, polEnd, nPol;
	
	int nChange=0;
	
	mono = DRng(cfg.rngState)*cs->polSize;
#if POL_TYPE == POL_LIN
	prevMono = (mono==0)?(cs->polSize-1):(mono-1);
#else
	prevMono = (mono-1+cs->polSize)%cs->polSize;
#endif
	polStart = (int)(DRng(cfg.rngState)*((cs->nPol-1)/BLOCK_SELECT+1))*BLOCK_SELECT;
	polEnd = MIN(cs->nPol,polStart+BLOCK_SELECT);
	nPol = polEnd-polStart;
	
	for(int i=polStart/8; i<(polEnd-1)/8+1; i++){
		cs->unit1Cache[i-polStart/8] = cs->unitPol[i+prevMono*cs->intsPerMono];
	}
	
	for(int i=polStart/8; i<(polEnd-1)/8+1; i++){
		cs->unit2Cache[i-polStart/8] = cs->unitPol[i+mono*cs->intsPerMono];
	}
	
	for(int i=polStart; i<polEnd; i++){
		cs->coorCache[i-polStart] = cs->coorPol[i+mono*cs->nPol];
	}
// 	printf("--------------------------------------\n");
	for(int iPol=0; iPol<nPol; iPol++){
		int uShift=(iPol%8)*4;
		if(uShift==0) 
			randOpt = Rng(cfg.rngState);
		src1 = (cs->unit1Cache[iPol/8]>>uShift)&0xf;
		src2 = (cs->unit2Cache[iPol/8]>>uShift)&0xf;
		option = src1*16*NOPT+src2*NOPT+((randOpt>>uShift)&0xf);
		
		dst1 = tab.transTable[option]&0xf;
		dst2 = (tab.transTable[option]>>4)&0xf;
		move = tab.transTable[option]>>8;
		
		if(dst1==0 && dst2 == 0) continue;
// 		printf("%x+%x => %x+%x\n", src1, src2, dst1, dst2); 
// 		if(iPol+polStart == cfg.efPolId && tab.accRat[option] < DRng(cfg.rngState)) continue;
		newCoor = AddUnitToCoor(move, cs->coorCache[iPol], &cs->con);
// 		printf("coor: %x->%x\n", cs->coorCache[iPol], newCoor);
		if(tab.addOccup[option]){
			int occup = OccupLattice(newCoor, cs);
			if(occup) continue;
			SetLattice(newCoor, cs);
		}
		if(tab.remOccup[option])
			UnsetLattice(cs->coorCache[iPol], cs);
		cs->coorCache[iPol] = newCoor;
		cs->unit1Cache[iPol/8] &= ~(0xf<<uShift); cs->unit1Cache[iPol/8] |= dst1<<uShift;
		cs->unit2Cache[iPol/8] &= ~(0xf<<uShift); cs->unit2Cache[iPol/8] |= dst2<<uShift;
		nChange++;
	}
	
	for(int i=polStart/8; i<(polEnd-1)/8+1; i++){
		cs->unitPol[i+prevMono*cs->intsPerMono] = cs->unit1Cache[i-polStart/8];
	}
	
	for(int i=polStart/8; i<(polEnd-1)/8+1; i++){
		cs->unitPol[i+mono*cs->intsPerMono] = cs->unit2Cache[i-polStart/8];
	}
	
	for(int i=polStart; i<polEnd; i++){
		cs->coorPol[i+mono*cs->nPol] = cs->coorCache[i-polStart];
	}
// 	if(CheckIntegrity(cs))
// 		exit(0);
	return nChange;
}


void CopyState(CurState* srcCS, CurState* dstCS){
	if(!dstCS->allocated)
		CSInit(dstCS, srcCS->con.BL, srcCS->polSize, srcCS->nPol, NULL);
	for(int coor=0; coor<cs->con.LAT_ALLOC; coor++)
		dstCS->lattice[coor] = srcCS->lattice[coor];
	for(int i=0; i<cs->polSize*cs->nPol; i++)
		dstCS->coorPol[i] = srcCS->coorPol[i];
	for(int i=0; i<cs->intsPerMono*cs->polSize; i++)
		dstCS->unitPol[i] = srcCS->unitPol[i];
}

int CharToHex(char c){
	int hex;
	
	hex = c-'0';
	if(hex>=10) hex = 10+(int)(c-'a');
	return hex;
}

double TRelax(int polSize){
	double tRouse, tInter, tRept, tau;
	
	if(POL_TYPE == POL_RING){
		tRouse = 0.08  *pow(polSize, 2.2);
		tInter = 0.015 *pow(polSize, 2.5);
		tRept  = 0.5e-3*pow(polSize, 2.97);
	}
	else{
		tRouse = 0.22    *pow(polSize, 2.3);
		tInter = 0.1    *pow(polSize, 2.5);
		tRept  = 0.4e-2*pow(polSize, 3.1);
	}
	
	tau = MAX(1e2, tRouse);
	tau = MAX(tau, tInter);
	tau = MAX(tau, tRept);
	return tau;
}


