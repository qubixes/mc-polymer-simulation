#include "oclpol.h"

uint GetGpuSite(uint t, uint u, uint v, SimProperties* sp){
	uint site,widt,widu,widv,wid;
	
	widt=t/(LCELL*WST);
	widu=u/(LCELL*WSU);
	widv=v/(LCELL*WSV);
	
	wid = widt+widu*sp->nwt+widv*sp->nwu*sp->nwt;
	
	site = ((t&0x3)<<9)|((u&0x3)<<11)|((v&0x3)<<13);
	site |= ((t&0x1c)>>2)|((u&0x1c)<<1)|((v&0x1c)<<4);
	site |= wid<<15;
	return site;
}


void CopyCPUToGPULattice(){
	uint gSite, cSlot, cSite;
	uint t,u,v, iDev;
	SimState* curState;
	
	for(iDev=0; iDev<sp.nDevices; iDev++){
		curState = ss+iDev;
		for(gSite=0; gSite<sp.gpuLatSize; gSite++) 
			curState->gpuLattice[gSite]=0;
		for(t=0; t<sp.LT; t++){
			for(u=0; u<sp.LU; u++){
				for(v=0; v<sp.LV; v++){
					gSite = GetGpuSite(t,u,v,&sp);
					GetSiteSlot(t,u,v,&cSite, &cSlot, &sp, MSPIN);
					curState->gpuLattice[gSite] |= (curState->lattice[cSite]>>(cSlot*8))&0xff;
					GetSiteSlot(t,u,v,&cSite, &cSlot, &sp, SL_MSPIN);
					curState->gpuLattice[gSite] |= ((curState->slLattice[cSite]>>cSlot)&0x1)<<8;
					GetSiteSlot(t,u,v,&cSite, &cSlot, &sp, LAB_MSPIN);
					curState->gpuLattice[gSite] |= ((curState->labLattice[cSite]>>(2*cSlot))&0x3)<<9;
				}
			}
		}
	}
}

void CopyGPUToCPULattice(){
	uint site, cSlot, cSite, gSite;
	uint t,u,v,iDev;
	SimState* curState;
	
	
	for(iDev=0; iDev<sp.nDevices; iDev++){
		curState = ss + iDev;
		
		for(site=0; site<sp.latSize; site++) curState->lattice[site]=0;
		for(site=0; site<sp.slSize; site++) curState->slLattice[site]=0;
		for(site=0; site<sp.labSize; site++) curState->labLattice[site]=0;
		for(t=0; t<sp.LT; t++){
			for(u=0; u<sp.LU; u++){
				for(v=0; v<sp.LV; v++){
					gSite = GetGpuSite(t,u,v,&sp);
					GetSiteSlot(t,u,v,&cSite, &cSlot, &sp, MSPIN);
					curState->lattice[cSite] |= (curState->gpuLattice[gSite]&0xff)<<(cSlot*8);
					GetSiteSlot(t,u,v,&cSite, &cSlot, &sp, SL_MSPIN);
					curState->slLattice[cSite] |= ((curState->gpuLattice[gSite]&0x100)>>8)<<cSlot;
					GetSiteSlot(t,u,v,&cSite, &cSlot, &sp, LAB_MSPIN);
					curState->labLattice[cSite] |= ((curState->gpuLattice[gSite]&0x600)>>9)<<(cSlot*2);
				}
			}
		}
	}
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

uint IsValid(uint a){
	return (a!=0 && a!=0x3 && a!=0xc && a!=0xf);
}
	

uint GetSl(uint* sl, uint t, uint u, uint v){
	uint site,slot;
	
	GetSiteSlot(t,u,v,&site,&slot,&sp,SL_MSPIN);
	return (sl[site]>>slot)&0x1;
}

uint GetLabel(uint* labelLat, uint t, uint u, uint v){
	uint site, slot;
	GetSiteSlot(t,u,v, &site, &slot, &sp, LAB_MSPIN);
	return ((labelLat[site]>>(2*slot))&0x3);
}
	

void AddUnitTUV(uint unit, uint* t, uint* u, uint* v){
	
	uint tu = unit&0x1;
	uint uu = (unit>>1)&0x1;
	uint vu = (unit>>2)&0x1;
	uint wu = (unit>>3)&0x1;
	
	*t += tu-wu;
	*u += uu-wu;
	*v += vu-wu;
	
	*t += sp.LT;
	*u += sp.LU;
	*v += sp.LV;
	
	*t %= sp.LT;
	*u %= sp.LU;
	*v %= sp.LV;
}

void SetBondVecs(SimState* ss, uint t, uint u, uint v, uint* bonds, uint polSize, uint* labels){
	int i;
	uint site,slot, slSite, slSlot, labSite, labSlot;
	GetSiteSlot(t,u,v,&site,&slot, &sp, MSPIN);
	
	for(i=0; i<polSize; i++){
		if(bonds[i]==0x0){
			GetSiteSlot(t,u,v,&slSite, &slSlot, &sp, SL_MSPIN);
			ss->slLattice[slSite] |= 0x1<<slSlot;
			GetSiteSlot(t,u,v,&labSite, &labSlot, &sp, LAB_MSPIN);
			ss->labLattice[labSite] |= (labels[i]*0x2)<<(2*labSlot);
// 			if(labels[i]) printf("%i ", i);
		}
		else if(bonds[i]==0xf)
			break;
		else{
			if((ss->lattice[site]>>(slot*8))&0x0f){
				printf("Error setting bond vec: site already occupied\n");
				printf("(%i %i %i)\n",t,u,v);
				exit(0);
			}
			ss->lattice[site] |= bonds[i]<<(slot*8);
			GetSiteSlot(t,u,v,&labSite,&labSlot, &sp, LAB_MSPIN);
			ss->labLattice[labSite] |= (labels[i]*0x1)<<(2*labSlot);
			AddUnitTUV(bonds[i],&t,&u,&v);
			GetSiteSlot(t,u,v,&site,&slot, &sp, MSPIN);
			if((ss->lattice[site]>>(slot*8))&0xf0){
				printf("Error setting bond vec: site already occupied[2]\n");
				printf("(%i %i %i)\n",t,u,v);
				exit(0);
			}
			ss->lattice[site] |= bonds[i]<<(4+slot*8);
		}
// 		if(labels[i]) printf("%i ", i);
	}
// 	printf("\n");
}

void SetPrevNextVecs(uint* gLattice, uint* gSl, uint* t, uint* u, uint* v, uint unit, uint sl){
	uint curId, slot;
	
	GetSiteSlot(*t, *u, *v, &curId, &slot, &sp, MSPIN);
	gLattice[curId] |= unit<<(slot*8);
	GetSiteSlot(*t, *u, *v, &curId, &slot, &sp, SL_MSPIN);
	gSl[curId] |= sl<<slot;
	
	AddUnitTUV(unit,t,u,v);
	GetSiteSlot(*t,*u,*v, &curId, &slot, &sp, MSPIN);
	gLattice[curId] |= (unit)<<(4+slot*8);
}

void GetSiteSlot(uint t, uint u, uint v, uint* site, uint* slot, SimProperties* sp, uint mSpin){
	if(t>=sp->LT || u>=sp->LU || v>=sp->LV){
		printf("Internal error: coordinates are out of range: %i %i %i\n", t,u,v);
		exit(1);
	}
	*site = t+u*sp->LT+v*sp->LT*sp->LU;
	*slot = (*site)%mSpin;
	*site /= mSpin;
}

void GetAllRingPolymers(){
	
	uint* latticeCopy;
	uint* slLattice, *labLattice;
	uint site, slot;
	uint t,u,v;
	uint n=0;
	int length;
	int iDev;
	SimState* curState;
	Polymer* pol;
	latticeCopy = (uint*) malloc(sizeof(uint)*sp.latSize);
	
	for(iDev=0; iDev<sp.nDevices; iDev++){
		curState = ss+iDev;
		memcpy(latticeCopy, curState->lattice, sp.latSize*sizeof(uint));
		slLattice = curState->slLattice;
		labLattice = curState->labLattice;
		n=0;
		for(t=0; t<sp.LT; t++){
			for(u=0; u<sp.LU; u++){
				for(v=0; v<sp.LV; v++){
					GetSiteSlot(t,u,v, &site, &slot, &sp, MSPIN);
					if((latticeCopy[site]>>(slot*8))&0xff){
						if(n==curState->nPol){
							printf("Too many polymers\n");
							return;
						}
						pol = curState->pol+n;
// 						printf("Pol number %i\n", n);
						length=GetRingPolymer(curState, t,u,v, pol, latticeCopy, slLattice, labLattice);
						pol->length=length;
// 						pol->startTUV.t = t;
// 						pol->startTUV.u = u;
// 						pol->startTUV.v = v;
						n++;
// 						if(length!= 4096) printf("Error in length: len=%i!\n", length);
					}
				}
			}
		}
// 		printf("found %i polymers\n", n);
	}
	free(latticeCopy);
}

int GetRingPolymer(SimState* ss, uint t, uint u, uint v, Polymer* pol, uint* lattice, uint* slLattice, uint* labLattice){
	uint tp, up, vp;
	
	FindPolymerStart(ss, &t,&u,&v, pol, lattice, slLattice, labLattice);
// 	printf("start at: %i %i %i\n", t,u,v);
	
	pol->startTUV.t = t; pol->startTUV.u = u; pol->startTUV.v = v;
	tp=t; up=u; vp=v;
	int iMono=0;
	Coor curCoor;
	uint unit, label;
	uint direction=DIR_FORWARD;
	curCoor = GetCoor(t,u,v);
// 	int i;
	label=0;
	int nSl=0;
	memset(pol->label,0,sp.maxPolLength*sizeof(uint));
	int leftSL=pol->startLeftSL;
	do{
		pol->polCoor[iMono] = curCoor;
		label = GetLabel(labLattice, t,u,v);
		if(leftSL==1){
			if(GetSl(slLattice, t,u,v)){
				pol->label[iMono] = label>>1;
				pol->bonds[iMono] = 0;
				nSl++;
			}
			else iMono--;
		}
		else{
			pol->label[iMono] = label&0x1;
			unit=GetNextSite(lattice, &t, &u, &v, direction);
			if(unit==0x0){
				pol->bonds[iMono] = 0xf;
				SetSite(lattice, tp,up,vp, 0);
				return iMono+1;
			}
			curCoor = AddUnitToCoor(unit,curCoor);
			SetSite(lattice, tp,up,vp, 0);
			pol->bonds[iMono] = unit;
			tp=t; up=u; vp=v;
		}
		iMono++; leftSL ^= 0x1;
		if(iMono-sp.maxPolLength <1) { printf("polymer too long\n"); return 0;}
	} while(!((pol->polCoor[0].x == curCoor.x) && (pol->polCoor[0].y==curCoor.y) && (pol->polCoor[0].z == curCoor.z) && leftSL == pol->startLeftSL));
// 	if(GetSl(slLattice,t,u,v) && !pol->startLeftSL){
// 		label = GetLabel(labLattice, t,u,v);
// 		pol->polCoor[iMono] = curCoor;
// 		pol->bonds[iMono] = 0x0;
// 		pol->label[iMono] = (label>>1)&0x1;
// 		iMono++;
// 		nSl++;
// 	}
	return iMono;
}
/*
int GetRingPolymer(SimState* ss, uint t, uint u, uint v, TUVCoor* sTUV, Coor* polCoor, uint* bonds, uint* labels, uint* lattice, uint* slLattice, uint* labLattice){
	uint tp, up, vp;
	
	FindPolymerStart(ss, &t,&u,&v,lattice,labLattice);
// 	printf("start at: %i %i %i\n", t,u,v);
	
	sTUV->t = t; sTUV->u = u; sTUV->v = v;
	tp=t; up=u; vp=v;
	int iMono=0;
	Coor curCoor;
	uint unit, label;
	uint direction=DIR_FORWARD;
	curCoor = GetCoor(t,u,v);
// 	int i;
	label=0;
	int nSl=0;
	memset(labels,0,sp.maxPolLength*sizeof(uint));
	do{
		polCoor[iMono] = curCoor;
		label = GetLabel(labLattice, t,u,v);
		if(GetSl(slLattice, t,u,v)){
			labels[iMono] = label>>1;
			polCoor[iMono+1] = curCoor;
			bonds[iMono] = 0;
			iMono++;
			nSl++;
		}
		else if(label&0x2) printf("!! ");
		labels[iMono] = label&0x1;
		unit=GetNextSite(lattice, &t, &u, &v, direction);
		
		if(unit==0x0){
			bonds[iMono] = 0xf;
			SetSite(lattice, tp,up,vp, 0);
			return iMono+1;
		}
		
		curCoor = AddUnitToCoor(unit,curCoor);
		SetSite(lattice, tp,up,vp, 0);
		bonds[iMono] = unit;
		tp=t; up=u; vp=v;
		iMono++;
		if(iMono-sp.maxPolLength <1) { printf("polymer too long\n"); return 0;}
	} while(!((polCoor[0].x == curCoor.x) && (polCoor[0].y==curCoor.y) && (polCoor[0].z == curCoor.z)));
	return iMono;
}*/


// void FindPolymerStart(SimState* ss, uint *t, uint *u, uint *v, uint* lattice, uint* labLattice){
// 	
// 	uint direction = DIR_BACKWARD;
// 	
// 	uint maxDLastLabel;
// 	uint dLastLabel=0;
// 	uint foundLabel=0;
// 	uint tLabel, uLabel, vLabel, label, unit;
// 	int i=0;
// 	
// 	
// 	for(maxDLastLabel=2; (1<<(maxDLastLabel/2))<ss->nPol; maxDLastLabel+=2);
// 		
// 	
// 	do{
// 		label = GetLabel(labLattice, *t,*u,*v);
// 		if(label&0x3){
// 			foundLabel=1;
// 			dLastLabel=0;
// 			tLabel = *t; uLabel = *u; vLabel=*v;
// 		}
// 		
// 		unit = GetNextSite(lattice, t, u, v, direction);
// 		
// 		if(unit==0x0) {
// // 			printf("Linear polymer detected\n");
// 			break;
// 		}
// 		
// 		if(foundLabel && dLastLabel > maxDLastLabel){
// 			*t = tLabel; *u = uLabel; *v = vLabel;
// 			break;
// 		}
// 		i++;
// 		dLastLabel++;
// 	}while(i<100000);
// 	if(i==100000) {
// 		printf("\nWoops: found=%i, dLast=%i, maxDLast=%i?!\n", foundLabel, dLastLabel, maxDLastLabel);
// 		exit(193);
// 	}
// }


void FindPolymerStart(SimState* ss, uint *t, uint *u, uint *v, Polymer* pol, uint* lattice, uint* slLattice, uint* labLattice){
	
	uint direction = DIR_BACKWARD;
	
	uint label, unit;
	int i=0;
	int curLabelHead=0;
	int sl;
	do{
		label = GetLabel(labLattice, *t,*u,*v);
		sl = GetSl(slLattice, *t, *u, *v);
		pol->startLeftSL=0;
		if(sl){
			curLabelHead>>=1;
			curLabelHead |= (label&0x1)<<(ss->headBit-1);
			i++;
			if(curLabelHead == ss->labelHead && i >=ss->headBit) break;
			label >>= 1;
			pol->startLeftSL=1;
		}
		curLabelHead>>=1;
		curLabelHead |= (label&0x1)<<(ss->headBit-1);
		i++;
		if(curLabelHead == ss->labelHead && i >=ss->headBit) break;
		
		unit = GetNextSite(lattice, t, u, v, direction);
		
		if(unit==0x0) {
// 			printf("Linear polymer detected\n");
			break;
		}
	}while(i<100000);
	if(i==100000) {
		printf("Error: Did not find label head\n");
// 		printf("\nWoops: found=%i, dLast=%i, maxDLast=%i?!\n", foundLabel, dLastLabel, maxDLastLabel);
		exit(193);
	}
}
Coor GetCoor(uint t, uint u, uint v){
	Coor c;
	
	c.x=(int)t-(int)u;
	c.y=(int)(t+u)-(int)v;
	c.z=(int)v;
	return c;
}

Coor DGetCoor(double t, double u, double v){
	Coor c;
	
	c.x=t-u;
	c.y=t+u-v;
	c.z=v;
	return c;
}



void SetSite(uint* gLattice, uint t, uint u, uint v, uint val){
	
	uint site, slot;
	
	if(val > 0xff){
		printf("Error error\n");
		exit(192);
	}
	GetSiteSlot(t,u,v, &site, &slot, &sp, MSPIN);
	
	gLattice[site] &= ~(0xff<<(8*slot));
	gLattice[site] |= (val<<(8*slot));
}


uint GetNextSite(uint* gLattice, uint* t, uint* u, uint* v, uint dirMode){
	uint curId, slot;
	uint next=0x0;
	
	GetSiteSlot(*t, *u, *v, &curId, &slot, &sp, MSPIN);
	if(dirMode == DIR_FORWARD){
		next = (gLattice[curId]>>(8*slot))&0xf;
		if(!next) {
// 			printf("0x%x\n", gLattice[curId]);
			return 0;
		}
	}
	else if(dirMode == DIR_BACKWARD){
		next =  ((~gLattice[curId])>>(4+8*slot))&0xf;
		if(next == 0xf) return 0;
	}
	
	AddUnitTUV(next, t,u,v);
	return next;
}




int RSQA(int* d){
	int dx, dy, dz;
	
	dx = d[0]-d[1];
	dy = d[0]+d[1]-d[2];
	dz = d[2];
	
	return (dx*dx+dy*dy+dz*dz);
}
	

int PolDistance(uint* start, uint* end, uint nMono){
	
	int i;
	int rsq=0, rsq2=0;
	
		
	for(i=0; i<nMono; i++){
		rsq += MinDistance(start+3*i, end+3*i);
		rsq2 += MinDistance(start+3*i, end+3*(nMono-i-1));
	}
	
	return (rsq<rsq2)?rsq:rsq2;
}

int MinDistance(uint* start, uint* end){
	int minRsq=1049230198;
	int q;
	int d[3];
	int l[3] = {3*WST*sp.nwt, 3*WSU*sp.nwu, 3*WSV*sp.nwv};
	for(q=0; q<27; q++){
		d[0] = start[0] - (int)end[0] + ((q%3)-1)*l[0];
		d[1] = start[1] - (int)end[1] + ((q%9)/3-1)*l[1];
		d[2] = start[2] - (int)end[2] + ((q/9)-1)*l[2];

		if(RSQA(d) <minRsq) minRsq = RSQA(d);
	}
	return minRsq;
}

double DRSQ(Coor* a, Coor* b){
	double dx,dy,dz;
	
	dx = a->x-b->x;
	dy = a->y-b->y;
	dz = a->z-b->z;
	
	return (dx*dx+dy*dy+dz*dz);
}

double Distance(Coor* a, Coor* b){
	int q;
	Coor aPrime;
	double d, dMin=1e10;
	int dt,du,dv;
	
	for(q=0; q<27; q++){
		dt=((q%3)-1)*sp.LT;
		du=(((q/3)%3)-1)*sp.LU;
		dv=((q/9)-1)*sp.LV;
		aPrime.x = a->x+dt-du;
		aPrime.y = a->y+dt+du-dv;
		aPrime.z = a->z+dv;
		
		d = DRSQ(&aPrime, b);
		dMin = (dMin<d)?dMin:d;
	}
	return dMin;
}

Coor AddUnitToCoor(uint unit, Coor coor){
	int t,u,v,w;
	
	w=(int)(unit&0x8)>>3;
	t=(int)(unit&0x1)-w;
	u=(int)((unit&0x2)>>1)-w;
	v=(int)((unit&0x4)>>2)-w;
	
	coor.x += t-u;
	coor.y += t+u-v;
	coor.z += v;
	return coor;
}

Coor CoorToBaseBlock(Coor coor){
	double t,u,v;
	
	t = 0.5*( coor.x+coor.y+coor.z);
	u = 0.5*(-coor.x+coor.y+coor.z);
	v = coor.z;
	
	while(t>sp.LT) t -= sp.LT;
	while(t<0) t += sp.LT;
	while(u>sp.LU) u -= sp.LU;
	while(u<0) u += sp.LU;
	while(v>sp.LV) v -= sp.LV;
	while(v<0) v += sp.LV;
	return DGetCoor(t,u,v);
}

Coor GetCMS(Coor* coor, uint polSize){
	int i;
	Coor cms;
	
	cms.x=0; cms.y=0; cms.z=0;
	
	for(i=0; i<polSize; i++){
		cms.x += coor[i].x;
		cms.y += coor[i].y;
		cms.z += coor[i].z;
	}
	cms.x /= polSize; cms.y /= polSize; cms.z /= polSize;
	cms = CoorToBaseBlock(cms);
	return cms;
}

int CEQ(Coor a, Coor b){
	return (a.x==b.x && a.y==b.y && a.z==b.z);
}

void SetSecAttPolymer(uint t, uint u, uint v, Polymer* pol){
	int i;
	
	pol->startTUV.t = t;
	pol->startTUV.u = u;
	pol->startTUV.v = v;
	
	pol->polCoor[0] = GetCoor(t,u,v);
	for(i=1; i<pol->length; i++)
		pol->polCoor[i] = AddUnitToCoor(pol->bonds[i-1], pol->polCoor[i-1]);
#if POL_TYPE == POL_RING
	if(!CEQ(AddUnitToCoor(pol->bonds[i-1], pol->polCoor[i-1]), pol->polCoor[0]))
		printf("Error: not a ring polymer!\n");
#endif
// 	pol->oldModes[0] = GetCMS(pol->polCoor, pol->length);
}

// void UpdatePolymerModes(SimState* ss){
// 	Coor cms;
// 	int* polModesTaken;
// 	double dif, minDif=1e20, minDifTaken=1e20;
// 	int minPol, iOldPol, iPol;
// 	Coor* newModes;
// 	Polymer* pol, *oldPol;
// 	
// 	newModes = (Coor*) malloc(sizeof(Coor)*ss->nPol);
// 	memset(newModes, 0, sizeof(Coor)*ss->nPol);
// 	
// 	polModesTaken = (int*) malloc(sizeof(int)*ss->nPol);
// 	memset(polModesTaken, 0, sizeof(int)*ss->nPol);
// 	
// 	for(iPol=0; iPol<ss->nPol; iPol++){
// 		pol = ss->pol+iPol;
// 		cms=GetCMS(pol->polCoor, pol->length);
// // 		printf("cms=%.1lf %.1lf %.1lf\n", cms.x, cms.y, cms.z);
// 		minPol=0;
// 		minDif=1e20; minDifTaken=1e20;
// 		for(iOldPol=0; iOldPol<ss->nPol; iOldPol++){
// 			oldPol = ss->pol+iOldPol;
// 			dif = Distance(&cms, oldPol->oldModes);
// 			if(dif<minDif){
// 				if(polModesTaken[iOldPol])
// 					minDifTaken=dif;
// 				else{
// 					minPol=iOldPol;
// 					minDif=dif;
// 				}
// 			}
// 		}
// 		if(minDifTaken<minDif)
// 			printf("Warning: polymer assignment is ambiguous!\n");
// // 		printf("d=%le\n", minDif);
// 		polModesTaken[minPol]=1;
// 		newModes[iPol] = cms;
// 		ss->polNumTranslate[iPol]=minPol;
// 	}
// 	for(iPol=0; iPol<ss->nPol; iPol++)
// 		ss->pol[iPol].oldModes[0] = newModes[iPol];
// 	free(newModes); free(polModesTaken);
// }

// void UpdatePolymerWithLabels(SimState* ss){
// 	int fNZIndex;
// 	int curZero, maxZero, maxNZIndex, shift;
// 	uint posEncodings, label, i;
// 	
// 	
// 	for(Polymer* pol = ss->pol; pol-ss->pol < ss->nPol; pol++){
// 		
// 		fNZIndex=0;
// 		for(i=0; !pol->label[i]; i++);
// 		fNZIndex = i;
// 		maxZero = fNZIndex;
// 		curZero = 0;
// 		for(i++; i<pol->length; i++){
// 			if(pol->label[i]){
// 				if(curZero>maxZero){
// 					maxNZIndex = i;
// 					maxZero=curZero;
// 				}
// 				curZero=0;
// 			}
// 			else
// 				curZero++;
// 		}
// 		if(fNZIndex+curZero>= maxZero){
// 			maxNZIndex = fNZIndex;
// 		}
// 		posEncodings=1; label=0; shift=0;
// 		for(i=2; posEncodings<=ss->nPol; i+=2, posEncodings<<=1){
// 			label |= (pol->label[(maxNZIndex+i)%pol->length])<<shift;
// 			shift++;
// 		}
// 		pol->labelStart = maxNZIndex;
// 		for(i=0; i<maxNZIndex; i++){
// 			AddUnitTUV(pol->bonds[i], &(pol->startTUV.t), &(pol->startTUV.u), &(pol->startTUV.v));
// 		}
// 		ss->polNumTranslate[pol-ss->pol] = label;
// // 		printf("first non-zero index = %i, label start = %i\n", fNZIndex, pol->labelStart);
// // 		printf("%li <==> %i\n", pol-ss->pol, label);
// 	}
// }	

void UpdatePolymerWithLabels(SimState* ss){
	int* polAccounted = (int*) malloc(sizeof(int)*ss->nPol);
	for(int i=0; i<ss->nPol; i++) polAccounted[i] = 0;
	for(Polymer* pol = ss->pol; pol-ss->pol < ss->nPol; pol++){
		int label=0;
		for(int i=ss->headBit; i<ss->headBit+ss->labelBit; i++){
			label |= pol->label[i]<<(ss->headBit+ss->labelBit-i-1);
		}
		int polId = ss->label2Id[label];
		if(polId<0) {
			printf("\n");
			PrintBonds(pol->label, pol->length);
			PrintBonds(pol->bonds, pol->length);
			printf("start=[%u %u %u]\n", pol->startTUV.t, pol->startTUV.u, pol->startTUV.v);
			printf("length = %i\n", pol->length);
			printf("Ouch...\n");
			exit(123);
		}
		ss->polNumTranslate[pol-ss->pol] = polId;
		polAccounted[polId]++;
	}
	int fail=0;
	for(int i=0; i<ss->nPol; i++){
		if(polAccounted[i] != 1){
			printf("Account failure: %i %i\n", i, polAccounted[i]);
			fail++;
		}
	}
	free(polAccounted);
	if(fail) exit(0);
}	


