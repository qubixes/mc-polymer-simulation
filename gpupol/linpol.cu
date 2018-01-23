#define TUD 2
#define BPW 3
#define WST 8
#define WSU 8
#define WSV 8

#define WS (WST*WSU*WSV)
#define CELL_LENGTH 4
#define CELL_SIZE (4*4*4)
#define BLOCK_SIZE (WS*CELL_SIZE)
#define BS_NUINT (BLOCK_SIZE/4)
#define WLT (WST*CELL_LENGTH)
#define WLU (WSU*CELL_LENGTH)
#define WLV (WSV*CELL_LENGTH)
#define TUV_MASK 0x421

// #define NWT 3
// #define NWU 3
// #define NWV 3
#define NWG (NWT*NWU*NWV)

#define T_MASK 0x0079
#define U_MASK 0x0782
#define V_MASK 0x7804

#define CELL_MASK 0x88f
#include <stdint.h>
#define uint uint32_t

// #define SLAB_MASK 0x3001

__device__ uint NextPrev(uint site, uint* lattice);
__device__ uint GlobalSite(uint lSite, uint lid);
__device__ uint EdgeOutside(uint site, uint next);
__device__ uint OutOfBounds(uint site, uint next);
__device__ uint OnEdge(uint site);
__device__ uint RelPos(uint unit);
__device__ uint AddUnitToSite(uint unit, uint site);
__device__ uint GPUValidateAddUnitVectors(uint a, uint b, uint* c);
__device__ uint GPUAddUnitVectors(uint a, uint b);
__device__ void TransMove(uint* lattice, uint site, uint* slab, uint* trans, uint4* rngState);
__device__ void DiffuseSL(uint* lattice, uint site, uint* slab);
__device__ int LocalSiteToSite(uint localSite);
__device__ void SetEdges(uint* slab, uint* lattice, uint lidSiteStart);
__device__ uint IsValidUnit(uint unit);


__device__ uint Rng4(uint4* state){
	uint b;
	b = ((((*state).x <<  6)^(*state).x) >> 13);
	(*state).x = ((((*state).x & 4294967294) << 18) ^b);
	b = ((((*state).y <<  2)^(*state).y) >> 27);
	(*state).y = ((((*state).y & 4294967288) <<  2) ^b);
	b = ((((*state).z << 13)^(*state).z) >> 21);
	(*state).z = ((((*state).z & 4294967280) <<  7) ^b);
	b = ((((*state).w <<  3)^(*state).w) >> 12);
	(*state).w = ((((*state).w & 4294967168) << 13)^b);
	return ((*state).x^(*state).y^(*state).z^(*state).w);
}

__device__ uint NextPrev(uint site, uint* lattice){
	return (lattice[site>>2]>>(8*(site&0x3)))&0xff;
}

__device__ uint AddUnitToSite(uint unit, uint site){
	uint flip;
	uint add;
	
	flip = (unit&0x8)?((unit&0x7)^0x7):unit;
	
	add = (unit&0x8)?((flip^site)&flip):(flip&site);
	add = ((add&0x1)<<3)|((add&0x2)<<6)|((add&0x4)<<9);
	site += (unit&0x8)?(-add):(add);
	site ^= flip;
	return site;
}

__device__ uint GPUValidateAddUnitVectors(uint a, uint b, uint* c){
	uint r, valid;
	if((a|b) != 0xf && (a&b))
		return 0;
	r = (((a|b)==0xf)?(a&b):(a|b));
	valid = (r==0x3||r==0xc)?0:1;
		
	*c = r;
	return valid;
}

__device__ uint GPUAddUnitVectors(uint a, uint b){
	return (((a|b)==0xf)?(a&b):(a|b));
}

__device__ uint IsValidUnit(uint unit){
	return (!(unit==0x0 || unit==0x3 || unit==0xc || unit==0xf));
}

/** Data structure:
 * The polymer's next and prev vectors, are in the "lattice[]" array.
 * lattice: [prev][next] x4 in blocks of 8 bit each.
 * The polymer's label and stored length data are combined in the "slab[]" array
 * A label is two bit, since it can be on either side of stored length.
 * The 
 * slab: [s][l [p][n]], with 8 labels (position 0-15), then 8 stored length (from position 16 and higher). 
 * 
**/
__device__ void TransMove(uint* lattice, uint site, uint* slab, uint* trans, uint4* rngState){
// 	uint next = NextPrev(site, lattice);
	uint next;
	uint prev;
	uint bond, temp;
	uint oldNextSite, oldPrevSite, bondSite, newSite;
	uint nNext, nPrev;
	uint prevBond, rand;
	uint latSiteComplete, nextPrev, slabSite, latNew, latBond, latPrev, latNext, newBond, slabNew, insBond;
	uint slabBond, label;
	
	latSiteComplete = lattice[site/4];
	nextPrev = (latSiteComplete>>(8*(site&0x3)))&0xff;
	prev = nextPrev>>4; next=nextPrev&0xf;
	slabSite = slab[site/8];
	
	if(!nextPrev) return;
	rand = Rng4(rngState);
	if((slabSite>>(16|(site&0x7)))&0x1){
		///Forward: SL-> no SL
		prevBond = (rand&0x4)>>2;
		bond = (nextPrev>>(prevBond*4))&0xf;
		if(bond==0){
			prevBond ^= 0x1;
			newBond = (rand>>3)&0xf;
			
			insBond = newBond^(0xf*(prevBond^0x1));
			
			if(!IsValidUnit(newBond)) return;
			if(EdgeOutside(site, newBond)) return;
			newSite = AddUnitToSite(newBond, site);
			latNew = lattice[newSite/4];
			if(latNew&(0xff<<((newSite&0x3)*8))) return;
			
			label = (slabSite>>((site&0x7)*2))&0x3;
			
			///Displace the label from the prev to the next site.
			slabSite &= ~(0x3<<((site&0x7)*2));
			
			slabSite |= ((label>>prevBond)&0x1)<<((site&0x7)*2);
			
			
			
			///Remove the stored length bit.
			slabSite &= ~(1<<(16|(site&0x7)));
			
			slab[site/8] = slabSite;
			
			
			slabNew = slab[newSite/8];
			slabNew |= ((label>>(prevBond^0x1))&0x1)<<((newSite&0x7)*2);
			
			slab[newSite/8] = slabNew;
			
			///Move end to new place.
			
			latSiteComplete |= insBond<<((site&0x3)*8+4*(prevBond^0x1));
			lattice[site/4] = latSiteComplete;
			
			if(site/4 == newSite/4)
				latNew = lattice[newSite/4];
			
			latNew |= insBond<<((newSite&0x3)*8+4*prevBond);
			lattice[newSite/4] = latNew;
			
			return;
		}
		
		
		if(EdgeOutside(site, bond^(0xf*prevBond)))
			return;
		
		nPrev = (trans[bond/4]>>(4*(2*(bond%4)+(rand&0x1))))&0xf;
		nNext = GPUAddUnitVectors((~nPrev)&0xf, bond);
		temp = nPrev;
		nPrev = (rand&0x2)?nPrev:nNext;
		nNext = (rand&0x2)?nNext:temp;
		
		bondSite = AddUnitToSite(bond^(0xf*prevBond), site);
		
		oldNextSite = (prevBond)?site:bondSite;
		oldPrevSite = (prevBond)?bondSite:site;
		
		newSite = AddUnitToSite(nPrev, oldPrevSite);
		
		if(EdgeOutside(oldPrevSite, nPrev))
			return;
		
		latNew = lattice[newSite/4];
		if(latNew&(0xff<<((newSite&0x3)*8)))
			return;
		
		latNew |= (nNext|(nPrev<<4))<<((newSite&0x3)*8);
		
		latSiteComplete &= ~(0xf<<((site&0x3)*8+4*prevBond));
		latSiteComplete |=  ((prevBond)?(nNext<<4):nPrev)<<((site&0x3)*8);
		
		if(newSite/4 == site/4){
			lattice[site/4] = (latSiteComplete&(0xff<<((site&0x3)*8))) | (latNew&(~(0xff<<((site&0x3)*8))));
		}
		else{
			lattice[site/4] = latSiteComplete;
			lattice[newSite/4] = latNew;
		}
		latBond = lattice[bondSite/4];
		
		latBond &= ~(0xf<<((bondSite&0x3)*8+4*(prevBond^0x1)));
		latBond |= ((prevBond)?(nPrev):(nNext<<4))<<((bondSite&0x3)*8);
		
		lattice[bondSite/4] = latBond;
		label = (slabSite >> (((site&0x7)*2)|prevBond))&0x1;
		slabSite &= ~((1<<((site&0x7)|16))|(1<<(((site&0x7)*2)|prevBond)));
		if(!prevBond){
			slabSite |= (slabSite&(1<<((site&0x7)*2+1)))>>1;
			slabSite &= ~(1<<((site&0x7)*2+1));
		}
		
		slab[site/8] = slabSite;
		slab[newSite/8] |= label<<((newSite&0x7)*2);
	}
	else{
		///Backward: no SL -> SL
		
		if(((rand>>4)&0x1) != 0) return;
		
		if(next == 0 || prev == 0){
			if(next==0){
				bond = prev^0xf;
				prevBond=1;
			}
			else{
				bond = next;
				prevBond=0;
			}
			
			if(EdgeOutside(site, bond)) return;
			newSite = AddUnitToSite(bond, site);
			slabBond = slab[newSite/8];
			if(slabBond & (1<<(16|(newSite&0x7))) ) return;
			
			label = (slabSite >> (((site&0x7)*2)))&0x1;
			slabSite &= ~(1<<((site&0x7)*2));
			
			slab[site/8] = slabSite;
			if(newSite/8 == site/8) slabBond = slab[newSite/8];
			
			label <<= prevBond^0x1;
			label |= ((slabBond >> ((newSite&0x7)*2))&0x1)<<prevBond;
			slabBond &= ~(0x3<<((newSite&0x7)*2));
			slabBond |= label <<((newSite&0x7)*2);
			slabBond |= (1<<(16|(newSite&0x7)));
			
			slab[newSite/8] = slabBond;
			
			latSiteComplete &= ~(0xf<<(8*(site&0x3)+4*prevBond));
			lattice[site/4]  = latSiteComplete;
			
			latNew = lattice[newSite/4];
			
			latNew &= ~(0xf<<((newSite&0x3)*8+4*(prevBond^0x1)));
			
			lattice[newSite/4] = latNew;
			return;
		}
		prevBond = rand&0x1;
		
		if(EdgeOutside(site, (~prev)&0xf)|EdgeOutside(site,next))
			return;
		
		if(!GPUValidateAddUnitVectors(next,prev,&bond))
			return;
		oldNextSite = AddUnitToSite(next, site);
		oldPrevSite = AddUnitToSite((~prev)&0xf, site);
		bondSite = (prevBond)?oldPrevSite:oldNextSite;
		
		if(bondSite == site) return;
		
		slabBond = slab[bondSite/8];
		if(slabBond & (1<<((bondSite&0x7)|16)))
			return;
		slab[bondSite/8] = slabBond | (1<<((bondSite&0x7)|16));
		
		slabSite = slab[site/8];
		label = (slabSite>>(2*(site&0x7)))&0x1;
		slab[site/8] = slabSite&(~(1<<(2*(site&0x7))));
		slabBond = slab[bondSite/8];
		if(prevBond){
			slabBond |= (slabBond&(1<<(2*(bondSite&0x7))))<<1;
			slabBond &= ~(1<<(2*(bondSite&0x7)));
		}
		slab[bondSite/8] = slabBond | (label<<((2*(bondSite&0x7))|(prevBond^0x1)));
		
		lattice[site/4] = latSiteComplete &(~(0xff<<((site&0x3)*8)));
		
		latPrev = lattice[oldPrevSite/4];
		latPrev &= ~(0xf<<((oldPrevSite%4)*8));
		latPrev |= (bond<<((oldPrevSite%4)*8));
		lattice[oldPrevSite/4] = latPrev;
		
		latNext = lattice[oldNextSite/4];
		latNext &= ~(0xf0<<((oldNextSite%4)*8));
		latNext |= (bond<<((oldNextSite%4)*8+4));
		lattice[oldNextSite/4] = latNext;
	}
}


__device__ void DiffuseSL(uint* lattice, uint site, uint* slab){
	uint next, prev, nextSite, tsl;
	uint bit, nBit, lBit, lnBit;
	uint slabSite, slabNext, label;
	
	next = NextPrev(site, lattice);
	prev = next>>4;
	next &= 0xf;
	
	if(!next || EdgeOutside(site, next))
		return;
	
	nextSite = AddUnitToSite(next,site);
	slabSite = slab[site/8]; slabNext = slab[nextSite/8];
	
	
	bit = 16|(site&0x7); nBit = 16|(nextSite&0x7);
	lBit = (site&0x7)*2; lnBit = (nextSite&0x7)*2;
	
	if((slabSite&(1<<bit)) && !(slabNext&(1<<nBit))){
		label = slabSite&(0x3<<lBit);
		slabSite &= ~(0x3<<lBit);
		slabSite |= (label&(0x2<<lBit))>>1;
		label = (label>>lBit)&0x1;
		slabNext |= label<<(lnBit+1);
	}
	else if(!(slabSite&(1<<bit)) && (slabNext&(1<<nBit))){
		label = (slabNext>>(lnBit+1))&0x1;
		slabNext &= ~(1<<(lnBit+1));
		slabSite |= (slabSite&(1<<lBit))<<1;
		slabSite &= ~(1<<lBit);
		slabSite |= label<<lBit;
	}
	else
		return;
	
	tsl = ((slabSite>>bit) & 0x1)<<nBit;
	slabSite &= ~(1<<bit);
	slabSite |= ((slabNext>>nBit)&0x1)<<bit;
	slabNext &= ~(1<<nBit);
	slabNext |= tsl;
	
	if(site/8 == nextSite/8){
		slab[site/8] = (slabSite&((0x3<<lBit)|(0x1<<bit)))|(slabNext&(~((0x3<<lBit)|(0x1<<bit))));
	}
	else{
		slab[site/8] = slabSite;
		slab[nextSite/8] = slabNext;
	}
}

__device__ uint EdgeOutside(uint site, uint next){
	int t,u,v,w;
	
	w = next>>3;
	t =  (next&0x1)    -w;
	u = ((next&0x2)>>1)-w;
	v = ((next&0x4)>>2)-w;
	
	if((site&T_MASK) == 0 && t<0)
		return 1;
	if((site&T_MASK) == T_MASK && t>0)
		return 1;
	if((site&U_MASK) == 0 && u<0)
		return 1;
	if((site&U_MASK) == U_MASK && u>0)
		return 1;
	if((site&V_MASK) == 0 && v<0)
		return 1;
	if((site&V_MASK) == V_MASK && v>0)
		return 1;
	return 0;
}

__device__ int LocalSiteToSite(uint localSite){
	
	localSite |= (localSite&0x10)<<3;
	localSite |= (localSite&0x20)<<6;
	localSite &= CELL_MASK;
	return localSite;
}

__device__ void SetEdges(uint* sledge, uint* lattice, uint lidSiteStart){
	uint site, lSite;
	uint next,prev;
	
// 	for(lSite=0; lSite<CELL_SIZE/16; lSite++)
// 		edged[lSite]=0;
	
	for(lSite=0; lSite<CELL_SIZE; lSite++){
// 		site = GlobalSite(lSite, lid);
		site = lidSiteStart|LocalSiteToSite(lSite);
		next = NextPrev(site, lattice);
		prev = (~(next>>4))&0xf;
		next &= 0xf;
		sledge[site>>3] |= EdgeOutside(site, next)<<(2*(site&0x7));
		sledge[site>>3] |= EdgeOutside(site, prev)<<(2*(site&0x7)+1);
	}
}


__global__ void polmove(uint nStep, uint4* seeds, uint* gLattice, uint* gTrans, uint tuvOffset, uint NWT, uint NWU, uint NWV){
	__shared__ uint lattice[BLOCK_SIZE/4];
	__shared__ uint slab[BLOCK_SIZE/8];
	uint trans[4];
	
	uint locId, localOffset,t,u,v, destId, memAddr, unit;
	uint tOffset, uOffset, vOffset;
	uint lid = threadIdx.x;
	uint wid = blockIdx.x;
	uint gid = wid * blockDim.x + lid;
	
	uint i;
	uint lidSiteStart;
	
	uint4 rngl;
	uint4 rngp;
	
	uint site;
	for(i=0; i<4; i++) trans[i] = gTrans[i];
	
	tOffset = tuvOffset%WLT;
	uOffset = (tuvOffset/WLT)%(WLU);
	vOffset = tuvOffset/(WLT*WLU);

	
	tOffset += ((lid&0x7)<<2)|((wid%NWT)*CELL_LENGTH*WST);
	uOffset += ((lid&0x38)>>1)|(((wid/NWT)%NWU)*CELL_LENGTH*WSU);
	vOffset += ((lid&0x1c0)>>4)|((wid/NWU/NWT)*CELL_LENGTH*WSV);
	
	localOffset = ((lid&0x7)<<4)|((lid&0x38)<<5)|((lid&0x1c0)<<6);
	
	for(locId=0; locId<CELL_SIZE/4; locId++){
		destId=(locId&0x3)|((locId&0x4)<<3)|((locId&0x8)<<6);
		destId |= localOffset>>2;
		lattice[destId]=0;
	}
	
	for(locId=0; locId<CELL_SIZE/8; locId++){
		destId=(locId&0x1)|((locId&0x2)<<3)|((locId&0x4)<<6);
		destId |= localOffset>>3;
		slab[destId]=0;
	}
	
	for(locId=0; locId<CELL_SIZE; locId++){
		t = tOffset+(locId&0x3);
		u = uOffset+((locId>>2)&0x3);
		v = vOffset+((locId>>4)&0x3);
		
		t %= NWT*CELL_LENGTH*WST;
		u %= NWU*CELL_LENGTH*WSU;
		v %= NWV*CELL_LENGTH*WSV;
		memAddr = ((t&0x3)<<9)|((t&0x1c)>>2) | ((u&0x3)<<11) | ((u&0x1c)<<1) | ((v&0x3)<<13) | ((v&0x1c)<<4);
		memAddr |= ((t>>5)+(u>>5)*NWT+(v>>5)*NWT*NWU)<<15;
		unit = gLattice[memAddr];
		destId = (locId&0x1)|((locId&0x2)<<2)|((locId&0x4)>>1)|((locId&0x8)<<4)|((locId&0x10)>>2)|((locId&0x20)<<6);
		
		destId |= localOffset;
		lattice[destId/4] |= (unit&0xff)<<(8*(destId%4));
		
		slab[destId/8] |= ((unit&0x100)>>8)<<(16|(destId&0x7));
		slab[destId/8] |= ((unit&0x600)>>9)<<((destId&0x7)*2);
	}
	
	lidSiteStart = ((lid&0x7)<<4)|((lid&0x38)<<5)|((lid&0x1c0)<<6);
// 	SetEdges(sledge, lattice, lidSiteStart);
	rngp = seeds[gid*2];
	rngl = seeds[gid*2+1];
	
	for(i=0; i<nStep; i++){
		__syncthreads();
		site = lidSiteStart | (Rng4(&rngl)&CELL_MASK);
		DiffuseSL(lattice, site, slab);
		__syncthreads();
		site = lidSiteStart | (Rng4(&rngl)&CELL_MASK);
		TransMove(lattice, site, slab, trans, &rngp);
	}
	__syncthreads();
	
	for(locId=0; locId<CELL_SIZE; locId++){
		t = tOffset+(locId&0x3);
		u = uOffset+((locId>>2)&0x3);
		v = vOffset+((locId>>4)&0x3);
		
		t %= NWT*CELL_LENGTH*WST;
		u %= NWU*CELL_LENGTH*WSU;
		v %= NWV*CELL_LENGTH*WSV;
		memAddr = ((t&0x3)<<9)|((t&0x1c)>>2) | ((u&0x3)<<11) | ((u&0x1c)<<1) | ((v&0x3)<<13) | ((v&0x1c)<<4);
		memAddr |= ((t>>5)+(u>>5)*NWT+(v>>5)*NWT*NWU)<<15;
		destId = (locId&0x1)|((locId&0x2)<<2)|((locId&0x4)>>1)|((locId&0x8)<<4)|((locId&0x10)>>2)|((locId&0x20)<<6);
		destId |= localOffset;
		
		unit = (lattice[destId/4]>>(8*(destId%4)))&0xff;
		unit |= ((slab[destId/8]>>(16|(destId&0x7)))&0x1)<<8;
		unit |= ((slab[destId/8]>>(2*(destId&0x7)))&0x3)<<9;
		gLattice[memAddr]=unit;
	}
	seeds[gid*2]=rngp;
	seeds[gid*2+1]=rngl;
	__syncthreads();
}
