#include <stdint.h>
#define uint uint32_t

#define WST 8
#define WSU 8
#define WSV 8

#define WS (WST*WSU*WSV)
#define CELL_LENGTH 3
#define CELL_SIZE (CELL_LENGTH*CELL_LENGTH*CELL_LENGTH)
#define BLOCK_SIZE (WS*CELL_SIZE)
#define WLT (WST*CELL_LENGTH)
#define WLU (WSU*CELL_LENGTH)
#define WLV (WSV*CELL_LENGTH)

#define WS_MASK  (WS-1)
#define TID_MASK (WST-1)
#define UID_MASK (WSU-1)
#define VID_MASK (WSV-1)

#include <stdint.h>
#define uint uint32_t

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


__device__ int TUVToIndex(int t, int u, int v){
	int index=0; 
	index += (t%3)<<9;
	index += ((u%3)*3)<<9;
	index += ((v%3)*9)<<9;
	
	index += (t/3)+(u/3)*WST+(v/3)*WST*WSU;
	return index;
}


__device__ int AddUnitToIndex(int unit, int index, int& OOB){
	int dw = (unit>>3)&0x1;
	int dt = ((unit>>0)&0x1)-dw;
	int du = ((unit>>1)&0x1)-dw;
	int dv = ((unit>>2)&0x1)-dw;
	
	int tr=(index>>9)%3;
	int ur=((index>>9)/3)%3;
	int vr=((index>>9)/9);
	
	tr += dt;
	ur += du;
	vr += dv;
	
	// And now for a new bag of tricks! The idea is to find out whether it is out of bounds after 
	// adding a combination of +/- t,u,v. We don't check for overflow, since that is actually useful for 
	// the subtraction itself (without the out-of-bounds). We know that we either go 
	// from 111 -> 000 or 111 -> 000 in case of out-of bounds.
	
	int mask  = (tr+1)/4*(0x001); //+t
	mask += (3-tr)/4*(0x007); //-t
	mask += (ur+1)/4*(0x008); //+u
	mask += (3-ur)/4*(0x038); //-u
	mask += (vr+1)/4*(0x040); //+v
	mask += (3-vr)/4*(0x1c0); //-v
	
	int newIndex = ((index&WS_MASK)+mask)&WS_MASK;
	int oldIndex = index&WS_MASK;
	
	int a = oldIndex|(oldIndex>>1)|(oldIndex>>2);
	int b = newIndex|(newIndex>>1)|(newIndex>>2);
	
	int c = oldIndex&(oldIndex>>1)&(oldIndex>>2);
	int d = newIndex&(newIndex>>1)&(newIndex>>2);
	
	OOB= (((~a)&d)|((~b)&c))&0x49;
	
	newIndex += (tr+ur*3+vr*9)*WS;
	
	return newIndex;
}

__device__ uint GPUValidateAddUnitVectors(int a, int b, int& c){
	int valid;
	if((a|b) != 0xf && (a&b))
		return 0;
	c = (((a|b)==0xf)?(a&b):(a|b));
	valid = (c==0x3||c==0xc)?0:1;
	
	return valid;
}

__device__ uint GPUAddUnitVectors(uint a, uint b){
	return (((a|b)==0xf)?(a&b):(a|b));
}

__device__ void TransForw(char* lattice, int index, uint* trans, uint4* rngState){
	int OOB;
	
	int latSiteComplete = lattice[index];
	if(!latSiteComplete) return;
	
	int next = latSiteComplete&0xf;
	int label = latSiteComplete&0x30;
	int sl = latSiteComplete&0x40;
	
	int newIndex = AddUnitToIndex(next, index, OOB);
	if(OOB) return;
	int newSiteComp = lattice[newIndex];
	int newSl=newSiteComp&0x40;
	if(sl+newSl==0) return;
	
	uint rand = Rng4(rngState);
	int newBond1 = (trans[next/4]>>(4*(2*next%4)+(next&0x1)))&0xf;
	int newBond2 = GPUAddUnitVectors((~newBond1)&0xf, next);
	
	int temp = newBond1;
	newBond1 = (rand&0x2)?newBond1:newBond2;
	newBond2 = (rand&0x2)?newBond2:temp;
	
	int destIndex = AddUnitToIndex(index,newBond1,OOB);
	if(OOB) return;
	int destSiteComp = lattice[destIndex];
	if(destSiteComp) return;
	
	int moveFirst;
	if(sl+newSl==0x80){
		moveFirst = (rand&0x4)>>2;
	}
	else if(sl)
		moveFirst = 1;
	else
		moveFirst = 0;
	
	destSiteComp = newBond2;
	if(moveFirst){
		latSiteComplete = newBond1|((label>>1)&0x10);
		destSiteComp |= label&0x10;
	}
	else{
		latSiteComplete = newBond1|label|sl;
		destSiteComp |= (newSiteComp&0x20)>>1;
		newSiteComp = newSiteComp&0x5f;
	}
	
	lattice[index] = latSiteComplete;
	lattice[destIndex] = destSiteComp;
	if(!moveFirst)
		lattice[newIndex] = newSiteComp;
}

__device__ void TransBack(char* lattice, int index, uint* trans, uint4* rngState){
	int OOB;
	int latSiteComplete = lattice[index];
	int next = latSiteComplete&0xf;
	int label = latSiteComplete&0x30;
	int sl = latSiteComplete&0x40;
	if(!latSiteComplete) return;
	int srcIndex = AddUnitToIndex(next, index, OOB);
	if(OOB) return;
	int srcSiteComp = lattice[srcIndex];
	int srcNext = srcSiteComp&0xf;
	int srcLabel= srcSiteComp&0x30;
	int srcSl = srcSiteComp&0x40;
	int newNext;
	if(srcSl) return;
	if(!GPUValidateAddUnitVectors(next, srcNext, newNext)) return;
	
	int newIndex = AddUnitToIndex(newNext, index, OOB);
	if(OOB) return;
	
	int newSiteComp = lattice[newIndex];
	int newSiteSl = newSiteComp&0x40;
	
	if(sl+newSiteSl == 0x80) return;
	
	uint rand = Rng4(rngState);
	
	int moveFirst;
	if(sl+newSiteSl == 0x0){
		moveFirst = rand&0x1;
	}
	else if(sl == 0x40)
		moveFirst = 0;
	else
		moveFirst = 1;
	
	if(moveFirst){
		latSiteComplete = newNext|(label<<1)|0x40;
	}
	else{
		latSiteComplete = newNext|label|sl;
		newSiteComp = (newSiteComp&0x3f)|(srcLabel<<1)|0x40;
	}
	
	lattice[srcIndex]=0;
	lattice[index] = latSiteComplete;
	lattice[newIndex] = newSiteComp;
}

__device__ void DiffuseSL(char* lattice, int index){
	int OOB;
	int latSiteComplete = lattice[index];
	int next = latSiteComplete&0xf;
	int label = latSiteComplete&0x30;
	int sl = latSiteComplete&0x40;
	if(!latSiteComplete) return;
	int newIndex = AddUnitToIndex(next, index, OOB);
	int newSiteComp = lattice[newIndex];
	int newSiteLabel = newSiteComp&0x30;
	int newSiteSl = newSiteComp&0x40;
	
	if(OOB) return;
	if(newSiteSl + sl != 0x40) return;
	
	if(sl){
		newSiteComp = newSiteComp | ((label&0x10)<<1) | 0x40;
		latSiteComplete = next|((label>>1)&0x10);
	}
	else{
		latSiteComplete = next|(label<<1)|((newSiteLabel>>1)&0x10)|0x40;
		newSiteComp = newSiteComp&0x1f;
	}
	
	lattice[index] = latSiteComplete;
	lattice[newIndex] = newSiteComp;
}


__global__ void polmove(int nStep,  uint4* seeds,  char* srcLattice, char* dstLattice, uint* gTrans, int dtuv, int dtuv_next, uint NWT, uint NWU, uint NWV){
	__shared__ char lattice[BLOCK_SIZE];
	uint trans[4];
	
	int lid = threadIdx.x;
	int wid = blockIdx.x;
	int gid = wid * blockDim.x + lid;
	
	
	uint4 rngl;
	uint4 rngp;
	
	uint site;
		
   int dt = dtuv%WLT;
	int du = (dtuv/WLT)%(WLU);
	int dv = dtuv/(WLT*WLU);
	
	int p=0;
	int pSwitchNext=dt*du*dv;
	int dtBlock=dt;
	int duBlock=du;
	int dvBlock=dv;
	int memOffSet=0;
	
	for(int src=wid*BLOCK_SIZE+lid*4; src<(wid+1)*BLOCK_SIZE; src += 4*WS){
		for(int i=0; i<4 && i+src<4*WS; i++){
			while(i+src>=pSwitchNext){
				memOffSet += pSwitchNext;
				dtBlock = (p&0x1)?(WLT-dt):dt;
				duBlock = (p&0x2)?(WLU-du):du;
				dvBlock = (p&0x4)?(WLV-dv):dv;
				pSwitchNext += dtBlock*duBlock*dvBlock;
				p++;
			}
			int offSet = src+i-memOffSet;
			int t = ((p&0x1)?dt:0) + (offSet%dtBlock);
			int u = ((p&0x2)?du:0) + ((offSet/dtBlock)%duBlock);
			int v = ((p&0x4)?dv:0) + (offSet/(dtBlock*duBlock));
			int index = TUVToIndex(t,u,v);
			lattice[index]=srcLattice[src+i];
		}
	}
	
	for(int i=0; i<4; i++) trans[i] = gTrans[i];
		
	int indexStart = ((lid&0x1f)<<2)|((lid&0x60)>>5)|(lid&0x180);
	
	rngp = seeds[gid*2];
	rngl = seeds[gid*2+1];
	
	__syncthreads();
	for(int i=0; i<nStep; i++){
		site = indexStart | ((Rng4(&rngl)%27)<<9);
		DiffuseSL(lattice, site); __syncthreads();
		
		site = indexStart | ((Rng4(&rngl)%27)<<9);
		TransForw(lattice, site, trans, &rngp); __syncthreads();
		
		site = indexStart | ((Rng4(&rngl)%27)<<9);
		DiffuseSL(lattice, site); __syncthreads();
		
		site = indexStart | ((Rng4(&rngl)%27)<<9);
		TransBack(lattice, site, trans, &rngp); __syncthreads();
	}
	
   dt = dtuv_next%WLT;
	du = (dtuv_next/WLT)%(WLU);
	dv = dtuv_next/(WLT*WLU);
	
	memOffSet=0;
	for(int p=0; p<8; p++){
		int dtBlock = (p^0x1)?(WLT-dt):dt;
		int duBlock = (p^0x2)?(WLU-du):du;
		int dvBlock = (p^0x4)?(WLV-dv):dv;
		int dstWid  =  ((wid%NWT)+NWT-(((p>>0)&0x1)^0x1))%NWT;
		    dstWid += (((wid%NWU)+NWU-(((p>>1)&0x1)^0x1))%NWU)*NWT;
		    dstWid += (((wid%NWU)+NWU-(((p>>2)&0x1)^0x1))%NWU)*NWT*NWU;
		for(int i=lid; i<dtBlock*duBlock*dvBlock; i+=WS){
			int t = i%dtBlock;
			int u = (i/dtBlock)%duBlock;
			int v = (i/(dtBlock*duBlock));
			
			int dst = dstWid*BLOCK_SIZE+memOffSet+i;
			int index = TUVToIndex(t, u, v);
			dstLattice[dst] = lattice[index];
		}
		memOffSet += dtBlock*duBlock*dvBlock;
	}
	__syncthreads();
}
