#include <stdio.h>
#include <stdlib.h>
#include "lowm_modes.h"

#define IS_MAIN
#include "rng.h"

typedef struct EnumConfig{
	int bonds[7];
}EnumConfig;

typedef struct Enumeration{
	EnumConfig cfgs[64][32];
	int nCfgs[64];
}Enumeration;

typedef struct PortCombos{
	int combo[16][9][3];
	int nCombo[16];
}PortCombos;

PortCombos* FindPortCombos();
void PrintEnumeration(Enumeration* er);
void DoStep(int lattice[8], int bonds[7], int pos, int startPos, int n, Enumeration* er);
Enumeration* EnumerateBox();
void UpscaleFile(char* file, char* fileOut, Enumeration* er, PortCombos* pc, int polType);
int GetRandomPorts(int bond, int portUsed[2], int ports[2], unsigned int* state, PortCombos* pc);

int main(int argc, char** argv){
	char* file, *fileOut;
	if(argc<3){
		printf("Not enough arguments\n");
		return 192;
	}
	
	file = argv[1];
	fileOut = argv[2];
	Enumeration* er = EnumerateBox();
	PortCombos* pc = FindPortCombos();
	UpscaleFile(file, fileOut, er, pc, POL_TYPE_RING);
}

PortCombos* FindPortCombos(){
	PortCombos* pc = malloc(sizeof(PortCombos));
	
	for(int bond=1; bond<0xf; bond++){
		if(!IsValid(bond)) continue;
		pc->nCombo[bond]=0;
		
		int w = (bond>>3)&0x1;
		int dt = (bond&0x1) - w;
		int du = ((bond>>1)&0x1) - w;
		int dv = ((bond>>2)&0x1) - w;
		
		for(int pos=0; pos<8; pos++){
			int tPos = pos&0x1;
			int uPos = (pos>>1)&0x1;
			int vPos = (pos>>2)&0x1;
			for(int nextBond=1; nextBond<0xf; nextBond++){
				if(!IsValid(nextBond)) continue;
				int ntPos = tPos + (nextBond&0x1)      - ((nextBond>>3)&0x1);
				int nuPos = uPos + ((nextBond>>1)&0x1) - ((nextBond>>3)&0x1);
				int nvPos = vPos + ((nextBond>>2)&0x1) - ((nextBond>>3)&0x1);
				
// 				if(bond == 0x7 && nextBond == 0x7) printf("%i %i %i [%i]\n", ntPos, nuPos, nvPos, pos);
				
				if( (ntPos+2)/2-1 != dt) continue;
				if( (nuPos+2)/2-1 != du) continue;
				if( (nvPos+2)/2-1 != dv) continue;
				
				int newTPos = (ntPos+2)%2;
				int newUPos = (nuPos+2)%2;
				int newVPos = (nvPos+2)%2;
				
				int newPos = newTPos+2*newUPos+4*newVPos;
				
				pc->combo[bond][pc->nCombo[bond]  ][0] = pos;
				pc->combo[bond][pc->nCombo[bond]  ][1] = newPos;
				pc->combo[bond][pc->nCombo[bond]++][2] = nextBond;
			}
		}
	}
	return pc;
}

void PrintEnumeration(Enumeration* er){
	for(int start=0; start<8; start++){
		for(int end=0; end<8; end++){
			int comb = start+8*end;
			printf("\n(%i->%i) [n=%i]\n", start, end, er->nCfgs[comb]);
			for(int iCfg=0; iCfg<er->nCfgs[comb]; iCfg++){
				for(int i=0; i<7; i++){
					printf("%x", er->cfgs[comb][iCfg].bonds[i]);
				}
				printf("\n");
			}
		}
	}
}

void DoStep(int lattice[8], int bonds[7], int pos, int startPos, int n, Enumeration* er){
	if(lattice[pos]) return;
	if(n==7){
		EnumConfig* eCfg = er->cfgs[startPos+8*pos]+er->nCfgs[startPos+8*pos];
		for(int i=0; i<7; i++){
			eCfg->bonds[i] = bonds[i];
		}
		er->nCfgs[startPos+8*pos]++;
		return;
	}
	
	lattice[pos]=1;
	int newPos;
	for(int bond=1; bond<15; bond++){
		if(!IsValid(bond)) continue;
		if(bond&0x8){
			if((~(bond|pos))&0x7) continue;
			newPos = (pos ^ (~bond))&0x7;
		}
		else{
			if(bond&pos) continue;
			newPos = pos|bond;
		}
		bonds[n]=bond;
		DoStep(lattice, bonds, newPos, startPos, n+1, er);
	}
	lattice[pos]=0;
}



Enumeration* EnumerateBox(){
	Enumeration* er = malloc(sizeof(Enumeration));
	int lattice[8];
	int bonds[7];
	
	for(int i=0; i<64; i++)
		er->nCfgs[i]=0;
	
	for(int i=0; i<8; i++)
		lattice[i]=0;
	
	for(int start=0; start<8; start++){
		DoStep(lattice, bonds, start, start, 0, er);
	}
// 	PrintEnumeration(er);
	return er;
}

/// Use bonds and ports with all stored length removed. That is way easier...
/// 

void SetPorts(int* bonds, int* ports, int len, unsigned int* state, PortCombos* pc){
	int portUsed[2], newCombo[3];
	int nSet;
	int totSet=0;
	
	for(int i=0; i<2*len; i++) ports[i]=-1;
	
	do{
		nSet=0;
		for(int prio=1, nPrioSet=0; prio<10 && !nPrioSet; prio++){
			for(int iPort=0; iPort<len; iPort++){
				if(ports[2*iPort]>=0) continue;
				portUsed[0] = ports[ (2*iPort-1+2*len)%(2*len) ];
				portUsed[1] = ports[ (2*iPort+2)%(2*len) ];
				int curPrio = GetRandomPorts(bonds[iPort], portUsed, newCombo, state, pc);
				if(curPrio == prio){ 
					for(int i=0; i<2; i++) ports[2*iPort+i] = newCombo[i]; 
					bonds[iPort] = newCombo[2];
					nSet++; 
					nPrioSet++;
					totSet++;
				}
			}
		}
	} while(nSet);
}

int GetRandomPorts(int bond, int portUsed[2], int combo[3], unsigned int* state, PortCombos* pc){
	int candidates[9][3];
	int nCandid=0;
	
	for(int i=0; i<pc->nCombo[bond]; i++){
		if(portUsed[0] == pc->combo[bond][i][0]) continue;
		if(portUsed[1] == pc->combo[bond][i][1]) continue;
		for(int k=0; k<3; k++){
			candidates[nCandid][k] = pc->combo[bond][i][k];
		}
		nCandid++;
	}
	
	int nDif[2] = {1,1};
	int last[2] = {candidates[0][0],candidates[0][1]};
	for(int i=1; i<nCandid; i++){
		for(int k=0; k<2; k++){
			if(last[k] != candidates[i][k]) nDif[k]++;
		}
	}
	int candidate = (int)(DRng(state)*nCandid);
	for(int k=0; k<3; k++)
		combo[k] = candidates[candidate][k];
	
	int prio;
	if(nCandid==1) prio=1;
	else prio = 1+MIN(nDif[0], nDif[1]);
	
	return prio;
}


void UpscaleFile(char* file, char* fileOut, Enumeration* er, PortCombos* pc, int polType){
	unsigned int state[4];
	Seed(state, 192869123);
	
	FILE* pFile = fopen(file, "r");
	FILE* pFileOut = fopen(fileOut, "w");
	if(!pFile) printf("woops: not opening file %s\n", file);
	if(!pFileOut) printf("Not opening file for writing: %s\n", fileOut);
	
// 	for(int i=0; i<3; i++) fscanf(pFile, "%*s %*s");
	int maxLen, len, t, u, v, nPol;
	fscanf(pFile, "%*s %i", &LT);
	fscanf(pFile, "%*s %i", &LU);
	fscanf(pFile, "%*s %i", &LV);
	fscanf(pFile, "%*s %i", &nPol);
	fscanf(pFile, "%*s %i", &maxLen);
	
	fprintf(pFileOut, "LT= %i\n", 2*LT);
	fprintf(pFileOut, "LU= %i\n", 2*LU);
	fprintf(pFileOut, "LV= %i\n", 2*LV);
	fprintf(pFileOut, "np= %i\n", nPol);
	fprintf(pFileOut, "maxPolLength= %i\n", 8*maxLen);
	fflush(pFileOut);
	char* str     = malloc(sizeof(char)*maxLen);
	int* bonds    = malloc(sizeof(int)*maxLen);
	int* nSl      = malloc(sizeof(int)*maxLen);
	int* newBonds = malloc(sizeof(int)*maxLen*8);
	int* ports    = malloc(sizeof(int)*maxLen*2);
	for(int iPol=0; iPol<nPol; iPol++){
		fscanf(pFile, "%*s %i", &len);
		fscanf(pFile, "%i %i %i", &t, &u, &v);
		fscanf(pFile, "%s", str);
		
		for(int i=0; i<=len; i++)
			nSl[i]=0;
		
		int bondLen=0;
		for(int i=0; i<len; i++){
			int bond= CharToHex(str[i]);
			if(!bond){
				nSl[bondLen]++;
			}
			else{
				bonds[bondLen++] = bond;
			}
		}
		nSl[0] += nSl[bondLen];
		for(int i=0; i<2*bondLen; i++) ports[i]=-1;
		
		SetPorts(bonds, ports, bondLen, state, pc);
		
		int iBondNew=0;
		///Just works for ring polymers. Patch it if you want linear polymers!
		for(int i=0; i<bondLen; i++){
			int start = ports[(i*2-1+2*bondLen)%(2*bondLen)];
			int end = ports[i*2];
			int combi = start+end*8;
			EnumConfig* eCfg = er->cfgs[combi] + (int)(DRng(state)*er->nCfgs[combi]);
			
			for(int ib=0; ib<7; ib++){
				newBonds[iBondNew++] = eCfg->bonds[ib];
				for(int iSl=0; iSl<nSl[i]; iSl++)
					newBonds[iBondNew++]=0;
			}
			
			newBonds[iBondNew++] = bonds[i];
			for(int iSl=0; iSl<nSl[i]; iSl++) newBonds[iBondNew++]=0;
		}
		
		fprintf(pFileOut, "len= %i\n", iBondNew);
		int tNew, uNew, vNew;
		tNew = 2*t+(ports[bondLen*2-1]&0x1);
		uNew = 2*u+((ports[bondLen*2-1]>>1)&0x1);
		vNew = 2*v+((ports[bondLen*2-1]>>2)&0x1);
		fprintf(pFileOut, "%i %i %i\n", tNew, uNew, vNew);
		for(int i=0; i<iBondNew; i++)
			fprintf(pFileOut, "%x", newBonds[i]);
		fprintf(pFileOut, "\n");
		fflush(pFileOut);
	}
	fclose(pFile); fclose(pFileOut);
}