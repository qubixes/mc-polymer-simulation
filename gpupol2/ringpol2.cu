#include <stdint.h>
#include <stdio.h>
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

__device__ uint Rng4(uint4& state){
	uint t=state.w;
	t^= t << 11;
	t^= t >> 8;
	state.w=state.z; state.z=state.y; state.y=state.x;
	t ^= state.x;
	t ^= state.x>>19;
	state.x=t;
	return t;
}

__device__ int TUVToIndex(int t, int u, int v){
	int index=0; 
	
	index += (t%3)<<9;
	index += ((u%3)*3)<<9;
	index += ((v%3)*9)<<9;
	
	index += (t/3)+(u/3)*WST+(v/3)*WST*WSU;
	return index;
}

__device__ void IndexToTUV(int index, int& t, int& u, int& v){
	t=(index>>9)%3;
	u=((index>>9)/3)%3;
	v=((index>>9)/9);
	
	t+= (index&0x7)*3;
	u+= ((index>>3)&0x7)*3;
	v+= ((index>>6)&0x7)*3;
}

__device__ int AddUnitToIndex(int unit, int index, int& OOB){
	int dt = ((unit>>0)&0x1);
	int du = ((unit>>1)&0x1);
	int dv = ((unit>>2)&0x1);
	int dw =  (unit>>3)&0x1;
	
	int t,u,v;
	IndexToTUV(index, t,u,v);
	
	t+= dt-dw;
	u+= du-dw;
	v+= dv-dw;
	
	OOB  = ((t+WLT)/WLT-1);
	OOB |= ((u+WLU)/WLU-1);
	OOB |= ((v+WLV)/WLV-1);
	int newIndex = TUVToIndex(t,u,v);
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

__device__ void TransForw(char* lattice, int index, uint* trans, uint4& rngState){
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
	int newBond1 = (trans[next/4]>>(4*(2*(next%4)+(rand&0x1))))&0xf;
	int newBond2 = GPUAddUnitVectors((~newBond1)&0xf, next);
	
	int temp = newBond1;
	newBond1 = (rand&0x2)?newBond1:newBond2;
	newBond2 = (rand&0x2)?newBond2:temp;
	
	int destIndex = AddUnitToIndex(newBond1,index, OOB);
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
	
// 	int t,u,v, tn, un, vn,td,ud,vd;
// 	IndexToTUV(index, t,u,v);
// 	IndexToTUV(newIndex,tn,un,vn);
// 	IndexToTUV(destIndex,td,ud,vd);

// 	printf("index=%i (%i,%i,%i) [%x], newIndex=%i (%i %i %i) [%x],  destIndex=%i (%i %i %i) [%x]\n", index, t,u,v, latSiteComplete, newIndex, tn,un,vn, newSiteComp, destIndex, td,ud,vd, destSiteComp);

	
	destSiteComp = newBond2;
	if(moveFirst){
		latSiteComplete = newBond1|((label>>1)&0x10);
		destSiteComp |= label&0x10;
	}
	else{
		latSiteComplete = newBond1|label|sl;
		destSiteComp |= (newSiteComp&0x20)>>1;
		newSiteComp = newSiteComp&0x1f;
	}
	
// 	printf("index=%i (%i,%i,%i) [%x], newIndex=%i (%i %i %i) [%x],  destIndex=%i (%i %i %i) [%x]: %x=>%x+%x, %x\n", index, t,u,v, latSiteComplete, newIndex, tn,un,vn, newSiteComp, destIndex, td,ud,vd, destSiteComp, next, newBond1, newBond2, (trans[next/4]>>(4*(2*(next%4))))&0xff);

	lattice[index] = latSiteComplete;
	lattice[destIndex] = destSiteComp;
	if(!moveFirst)
		lattice[newIndex] = newSiteComp;
}

__device__ void TransBack(char* lattice, int index, uint* trans, uint4& rngState){
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
		latSiteComplete = newNext|(label<<1)|srcLabel|0x40;
	}
	else{
		latSiteComplete = newNext|label|sl;
		newSiteComp = (newSiteComp&0x3f)|(srcLabel<<1)|0x40;
	}
// 	printf("%x + %x -> %x\n", next, srcNext, newNext);
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
	if(OOB) return;
	int newSiteComp = lattice[newIndex];
	int newSiteLabel = newSiteComp&0x30;
	int newSiteSl = newSiteComp&0x40;
	if(newSiteSl + sl != 0x40) return;
	
	if(sl){
		newSiteComp = newSiteComp | ((label&0x10)<<1) | 0x40;
		latSiteComplete = next|((label>>1)&0x10);
	}
	else{
		latSiteComplete = next|(label<<1)|((newSiteLabel>>1)&0x10)|0x40;
		newSiteComp = newSiteComp&0x1f;
	}
// 	if(!sl){
// 	int t,u,v, tn, un, vn;
// 	IndexToTUV(index, t,u,v);
// 	IndexToTUV(newIndex,tn,un,vn);

// 	printf("sl=%i, next=%x, index=%i (%i,%i,%i), newIndex=%i (%i %i %i)\n", sl, next, index, t,u,v, newIndex, tn,un,vn);
	lattice[index] = latSiteComplete;
	lattice[newIndex] = newSiteComp;
// 	}
}


__global__ void polmove(int nStep,  uint4* seeds,  char* srcLattice, char* dstLattice, uint* gTrans, int dtuv, int dtuv_next, uint NWT, uint NWU, uint NWV){
	__shared__ char lattice[BLOCK_SIZE];
	uint trans[4];
	
	int lid = threadIdx.x;
	int wid = blockIdx.x;
	int gid = wid * blockDim.x + lid;
	int widt = wid%NWT;
	int widu = (wid/NWT)%NWU;
	int widv = wid/(NWU*NWT);
	
	uint4 rngl;
	uint4 rngp;
	
	uint site;
		
   int dt = dtuv%WLT;
	int du = (dtuv/WLT)%(WLU);
	int dv = dtuv/(WLT*WLU);
	
	int p=0;
	int dtBlock=WLT-dt;
	int duBlock=WLU-du;
	int dvBlock=WLV-dv;
	int pSwitchNext=dtBlock*duBlock*dvBlock;
	int memOffSet=0;
	
// 	printf("pSwitchNext=%i\n", pSwitchNext);
	
	int src;
	for(src=lid*4; src<BLOCK_SIZE; src += 4*WS){
		for(int i=0; i<4 && i+src<BLOCK_SIZE; i++){
			while(i+src>=pSwitchNext){
				memOffSet = pSwitchNext;
				p++;
				dtBlock = (p&0x1)?dt:(WLT-dt);
				duBlock = (p&0x2)?du:(WLU-du);
				dvBlock = (p&0x4)?dv:(WLV-dv);
				pSwitchNext += dtBlock*duBlock*dvBlock;
			}
			int offSet = src+i-memOffSet;
			int t = ((p&0x1)?(WLT-dt):0) + (offSet%dtBlock);
			int u = ((p&0x2)?(WLU-du):0) + ((offSet/dtBlock)%duBlock);
			int v = ((p&0x4)?(WLV-dv):0) + (offSet/(dtBlock*duBlock));
			int index = TUVToIndex(t,u,v);
			lattice[index]=srcLattice[src+i+wid*BLOCK_SIZE];
		}
	}
	
	for(int i=0; i<4; i++) trans[i] = gTrans[i];
		
	int indexStart = ((lid&0x1f)<<2)|((lid&0x60)>>5)|(lid&0x180);
	
	rngp = seeds[gid*2];
	rngl = seeds[gid*2+1];
	
	__syncthreads();
	for(int i=0; i<nStep; i++){
// 		site = indexStart | ((Rng4(rngl)%27)<<9);
// 		DiffuseSL(lattice, site); __syncthreads();
		
		uint randLoc;
		do {
			randLoc = Rng4(rngl);
		}while(randLoc>=4294574721); /// 8081*27^4, so that we are getting good random numbers.
		
// 		site = indexStart | ((Rng4(rngl)%27)<<9);
// 		TransForw(lattice, site, trans, rngp); __syncthreads();
// 		
// 		site = indexStart | ((Rng4(rngl)%27)<<9);
// 		DiffuseSL(lattice, site); __syncthreads();
// 		
// 		site = indexStart | ((Rng4(rngl)%27)<<9);
// 		TransForw(lattice, site, trans, rngp); __syncthreads();
// 		
// 		site = indexStart | ((Rng4(rngl)%27)<<9);
// 		TransBack(lattice, site, trans, rngp); __syncthreads();
		
		site = indexStart | ((randLoc%27)<<9);
		TransForw(lattice, site, trans, rngp); __syncthreads();
		randLoc /= 27;
		site = indexStart | ((randLoc%27)<<9);
		DiffuseSL(lattice, site); __syncthreads();
		randLoc /= 27;
		site = indexStart | ((randLoc%27)<<9);
		TransForw(lattice, site, trans, rngp); __syncthreads();
		randLoc /= 27;
		site = indexStart | ((randLoc%27)<<9);
		TransBack(lattice, site, trans, rngp); __syncthreads();
	}
	
   dt = dtuv_next%WLT;
	du = (dtuv_next/WLT)%(WLU);
	dv = dtuv_next/(WLT*WLU);
	
	memOffSet=0;
// 	printf("????\n");
	for(int p=0; p<8; p++){
		int dtBlock = (p&0x1)?dt:(WLT-dt);
		int duBlock = (p&0x2)?du:(WLU-du);
		int dvBlock = (p&0x4)?dv:(WLV-dv);
		int dstWid  =  (widt+NWT-(((p>>0)&0x1)))%NWT;
		    dstWid += ((widu+NWU-(((p>>1)&0x1)))%NWU)*NWT;
		    dstWid += ((widv+NWV-(((p>>2)&0x1)))%NWV)*NWT*NWU;
// 			if(lid==0)
// 				printf("p=%i, wid=(%i,%i,%i), dstWid=(%i,%i,%i)=%i\n", p,widt,widu,widv,(widt+NWT-(((p>>0)&0x1)))%NWT, ((widu+NWU-(((p>>1)&0x1)))%NWU), ((widv+NWV-(((p>>2)&0x1)))%NWV), dstWid);
// 		if(lid==0 && wid==0)
// 			printf("block=(%i,%i,%i), p=%i\n", dtBlock, duBlock, dvBlock, p);
		for(int i=lid; i<dtBlock*duBlock*dvBlock; i+=WS){
			int t = i%dtBlock + ((p&0x1)?0:dt);
			int u = (i/dtBlock)%duBlock + ((p&0x2)?0:du);
			int v = i/(dtBlock*duBlock) + ((p&0x4)?0:dv);
			
			int dst = dstWid*BLOCK_SIZE+memOffSet+i;
			int index = TUVToIndex(t, u, v);
// 			if(lid%55==0)
// 				printf("dstWid=%i,%i (p=%i), memOffSet=%i, i=%i, (%i,%i,%i)\n", dstWid, dst, p, memOffSet, i, t,u,v);
			dstLattice[dst] = lattice[index];
		}
		memOffSet += dtBlock*duBlock*dvBlock;
	}
	seeds[gid*2]=rngp;
	seeds[gid*2+1]=rngl;
	__syncthreads();
}
