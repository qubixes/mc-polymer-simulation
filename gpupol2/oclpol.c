#include "oclpol.h"

uint GetGpuSite(uint t, uint u, uint v, SimProperties* sp){
	uint widt=t/(LCELL*WST);
	uint widu=u/(LCELL*WSU);
	uint widv=v/(LCELL*WSV);
	
	uint lidt=t%(LCELL*WST);
	uint lidu=u%(LCELL*WSU);
	uint lidv=v%(LCELL*WSV);
	
	uint wid = widt+widu*sp->nwt+widv*sp->nwu*sp->nwt;
	uint lid = lidt+lidu*(LCELL*WST)+lidv*(LCELL*WST*LCELL*WSU);
// 	printf("(%u,%u,%u) -> %i (max = %i)\n",t,u,v,lid+wid*(CELL_SIZE*WS), sp->latSize);
	return lid+wid*(CELL_SIZE*WS);
}

void CopyGPUToCPULattice(char* gpuLattice, char* cpuLattice, int tOff, int uOff, int vOff, int dt, int du, int dv, SimProperties* sp){
	
	for(int widt=0; widt<sp->nwt; widt++){
		for(int widu=0; widu<sp->nwu; widu++){
			for(int widv=0; widv<sp->nwv; widv++){
				int wid=widt + widu*sp->nwt + widv*sp->nwt*sp->nwu;
				int memOffSet=0; 
				int curTOff=tOff+widt*LCELL*WST;
				int curUOff=uOff+widu*LCELL*WSU;
				int curVOff=vOff+widv*LCELL*WSV;
				for(int p=0; p<8; p++){
					int dtBlock = (p&0x1)?dt:(WLT-dt);
					int duBlock = (p&0x2)?du:(WLU-du);
					int dvBlock = (p&0x4)?dv:(WLV-dv);
					for(int i=0; i<dtBlock*duBlock*dvBlock; i++){
						int t = ((p&0x1)?(WLT-dt):0) + (i%dtBlock);
						int u = ((p&0x2)?(WLU-du):0) + ((i/dtBlock)%duBlock);
						int v = ((p&0x4)?(WLV-dv):0) + (i/(dtBlock*duBlock));
						t = (t+curTOff)%sp->LT;
						u = (u+curUOff)%sp->LU;
						v = (v+curVOff)%sp->LV;
						uint site = GetGpuSite(t,u,v,sp);
						cpuLattice[site] = gpuLattice[i+memOffSet+wid*WS*CELL_SIZE];
					}
					memOffSet += dtBlock*duBlock*dvBlock;
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
	

char GetSl(char* lattice, uint t, uint u, uint v, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	return GetSlSite(lattice, site);
}

char GetSlSite(char* lattice, uint site){
	return (lattice[site]>>6)&0x1;
}

void SetSl(char* lattice, uint t, uint u, uint v, int sl, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	SetSlSite(lattice, site, sl);
}

void SetSlSite(char* lattice, uint site, int sl){
	lattice[site] &= ~(1<<6);
	lattice[site] |= (sl<<6);
}

uint GetBond(char* lattice, uint t, uint u, uint v, int sl, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	return GetBondSite(lattice, site);
}

uint GetBondSite(char* lattice, uint site){
	return (lattice[site]&0xf);
}

void SetBond(char* lattice, uint t, uint u, uint v, int bond, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	SetBondSite(lattice, site, bond);
}

void SetBondSite(char* lattice, uint site, int bond){
	lattice[site] &= ~(0xf);
	lattice[site] |= bond;
}

uint GetLabel(char* lattice, uint t, uint u, uint v, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	return GetLabelSite(lattice, site);
}

uint GetLabelSite(char* lattice, uint site){
	return (lattice[site]>>4)&0x3;
}

void SetLabel(char* lattice, uint t, uint u, uint v, uint label, uint left, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	SetLabelSite(lattice, site, label, left);
}

void SetLabelSite(char* lattice, uint site, uint label, uint left){
	lattice[site] &= ~(1<<(4+left));
	lattice[site] |= label<<(4+left);
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

void SetLatticeFromSS(SimState* ss, SimProperties* sp){
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		memset(ss[iDev].lattice, 0, sizeof(char)*sp->latSize);
		for(int iPol=0; iPol<ss[iDev].nPol; iPol++){
			Polymer* pol = ss[iDev].pol+iPol;
			SetBondVecs(ss[iDev].lattice, pol->startTUV.t, pol->startTUV.u, pol->startTUV.v, pol->bonds, pol->length, pol->label, sp);
		}
	}
}


void SetBondVecs(char* lattice, uint t, uint u, uint v, uint* bonds, uint polSize, uint* labels, SimProperties* sp){
	uint site = GetGpuSite(t,u,v,sp);
	for(int i=0; i<polSize; i++){
		if(bonds[i]==0x0){
			SetSlSite(lattice, site, 1);
			SetLabelSite(lattice, site, labels[i], 1);
		}
		else if(bonds[i]==0xf){
			printf("Error: Linear polymer detected!\n");
			break;
		}
		else{
			if(GetBondSite(lattice, site)){
				printf("Error setting bond vector %i: site already occupied\n", i);
				printf("(%i %i %i)\n",t,u,v);
				exit(0);
			}
			SetBondSite(lattice, site, bonds[i]);
			SetLabelSite(lattice, site, labels[i], 0);
// 			printf("Setting: %i -> (%x, %x)\n", site, bonds[i], labels[i]);
			AddUnitTUV(bonds[i],&t,&u,&v);
			site = GetGpuSite(t,u,v,sp);
		}
	}
}

void GetAllRingPolymers(){
	
	char* latticeCopy;
	uint t,u,v;
	uint n=0;
	int length;
	int iDev;
	SimState* curState;
	Polymer* pol;
	latticeCopy = (char*) malloc(sizeof(uint)*sp.latSize);
	
	for(iDev=0; iDev<sp.nDevices; iDev++){
		curState = ss+iDev;
		memcpy(latticeCopy, curState->lattice, sp.latSize*sizeof(char));
		n=0;
		for(t=0; t<sp.LT; t++){
			for(u=0; u<sp.LU; u++){
				for(v=0; v<sp.LV; v++){
					uint site = GetGpuSite(t,u,v,&sp);
					if(GetBondSite(latticeCopy, site)){
						if(n==curState->nPol){
							printf("Too many polymers\n");
							return;
						}
						pol = curState->pol+n;
						length=GetRingPolymer(curState, t,u,v, pol, latticeCopy);
						pol->length=length;
						n++;
					}
				}
			}
		}
// 		printf("found %i polymers\n", n);
	}
	free(latticeCopy);
}

int GetRingPolymer(SimState* ss, uint t, uint u, uint v, Polymer* pol, char* lattice){
	FindPolymerStart(ss, &t,&u,&v, pol, lattice);
// 	printf("start at: %i %i %i\n", t,u,v);
	
	pol->startTUV.t = t; pol->startTUV.u = u; pol->startTUV.v = v;
	int iMono=0;
// 	Coor curCoor;
	uint unit;
// 	curCoor = GetCoor(t,u,v);
// 	int i;
	int nSl=0;
	memset(pol->label,0,sp.maxPolLength*sizeof(uint));
	int leftSL=pol->startLeftSL;
	uint site = GetGpuSite(t,u,v,&sp);
	uint label = GetLabelSite(lattice, site);
	do{
		if(leftSL==1){
			if(GetSlSite(lattice, site)){
				pol->label[iMono] = label>>1;
				pol->bonds[iMono] = 0;
				nSl++;
			}
			else iMono--;
		}
		else{
			pol->label[iMono] = label&0x1;
			unit=GetBondSite(lattice, site);
			if(unit==0x0){
				pol->bonds[iMono] = 0xf;
				SetBondSite(lattice, site, 0);
				return iMono+1;
			}
			SetBondSite(lattice, site, 0);
			pol->bonds[iMono] = unit;
			AddUnitTUV(unit, &t, &u, &v);
			site = GetGpuSite(t,u,v, &sp);
			label = GetLabelSite(lattice, site);
		}
		iMono++; leftSL ^= 0x1;
		if(iMono-sp.maxPolLength <1) { printf("polymer too long\n"); exit(0);}
	} while(!((pol->startTUV.t == t) && (pol->startTUV.u == u) && pol->startTUV.v == v && leftSL == pol->startLeftSL));
	return iMono;
}

void PrintLabelHead(int label){
	for(int i=4; i>=0; i--){
		printf("%i", (label>>i)&0x1);
	}
	printf("\t");
}

void PrintLabel(int label){
	for(int i=1; i>=0; i--){
		printf("%i", (label>>i)&0x1);
	}
	printf("\n");
}

void FindPolymerStart(SimState* ss, uint *t, uint *u, uint *v, Polymer* pol, char* lattice){
	uint label, unit;
	int i=0;
	int curLabelHead=0;
	int headMask = (1<<ss->headBit)-1;
	int sl;
	uint site = GetGpuSite(*t,*u,*v,&sp);
	do{
		sl = GetSlSite(lattice, site);
		label = GetLabelSite(lattice, site);
		pol->startLeftSL=sl?1:0;
		
		curLabelHead<<=1;
		curLabelHead &= headMask;
		curLabelHead |= (label>>sl)&0x1;
// 		PrintLabelHead(curLabelHead); PrintLabel(label);
		i++;
		if(curLabelHead == ss->labelHead && i >=ss->headBit) break;
		
		if(sl){
			curLabelHead <<= 1;
			curLabelHead &= headMask;
			curLabelHead |= (label&0x1);
// 			PrintLabelHead(curLabelHead);  PrintLabel(label);
			i++;
			pol->startLeftSL=0;
			if(curLabelHead == ss->labelHead && i >= ss->headBit) break;
			label >>= 1;
		}
		unit = GetBondSite(lattice, site);
		if(unit==0x0) {
			printf("Linear polymer detected: not supported\n");
			exit(0);
			break;
		}
		AddUnitTUV(unit, t,u,v);
		site = GetGpuSite(*t,*u,*v, &sp);
// 		printf("unit=%x, (%i,%i,%i), lattice: %hhx\n", unit, *t,*u,*v, lattice[site]);
	}while(i<2*sp.maxPolLength);
// 	printf("\n\n");
	if(i==2*sp.maxPolLength) {
		printf("Error: Did not find label head\n");
		for(int v=0; v<sp.LV; v++){
			for(int u=0; u<sp.LU; u++){
				for(int t=0; t<sp.LT; t++){
					printf("%hhx ", ss[0].lattice[GetGpuSite(t,u,v,&sp)]);
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("\n");
// 		printf("\nWoops: found=%i, dLast=%i, maxDLast=%i?!\n", foundLabel, dLastLabel, maxDLastLabel);
		exit(193);
	}
	///Advance one more bond.
// 	printf("curLabelHead:"); PrintLabelHead(curLabelHead);printf("\n");
	if(pol->startLeftSL) pol->startLeftSL=0;
	else{
		unit = GetBondSite(lattice,site);
		AddUnitTUV(unit,t,u,v);
		pol->startLeftSL = GetSl(lattice, *t,*u,*v, &sp);
	}
// 	printf("tuv=(%i %i %i)\n", *t, *u, *v);
}


void UpdatePolymerWithLabels(SimState* ss){
	int* polAccounted = (int*) malloc(sizeof(int)*ss->nPol);
	for(int i=0; i<ss->nPol; i++) polAccounted[i] = 0;
	for(Polymer* pol = ss->pol; pol-ss->pol < ss->nPol; pol++){
		int label=0;
		for(int i=0; i<ss->labelBit; i++){
			label |= pol->label[i]<<(ss->labelBit-i-1);
		}
		int polId = ss->label2Id[label];
		if(polId<0) {
			printf("\n");
			PrintBonds(pol->label, pol->length);
			PrintBonds(pol->bonds, pol->length);
			printf("start=[%u %u %u]\n", pol->startTUV.t, pol->startTUV.u, pol->startTUV.v);
			printf("length = %i\n", pol->length);
			printf("label = %i\n", label);
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
