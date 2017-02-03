#include "efpol.h"

void LoadPolymer(PolyConfig* pcfg, CurState* cs, int iPol){
	CoorI* tuv = pcfg->tuv;
	CoorI* xyz = pcfg->xyz;
	
	tuv[0].x[0] = TCoor(cs->coorPol[iPol], &cs->con);
	tuv[0].x[1] = UCoor(cs->coorPol[iPol], &cs->con);
	tuv[0].x[2] = VCoor(cs->coorPol[iPol], &cs->con);
	
	xyz[0].x[0] = X(tuv[0].x[0], tuv[0].x[1], tuv[0].x[2]);
	xyz[0].x[1] = Y(tuv[0].x[0], tuv[0].x[1], tuv[0].x[2]);
	xyz[0].x[2] = Z(tuv[0].x[0], tuv[0].x[1], tuv[0].x[2]);
	
	
// 	printf("iPol=%i\n", iPol);
	for(int iMono=1; iMono<cs->polSize; iMono++){
		int unit = (cs->unitPol[iPol/8+iMono*cs->intsPerMono]>>((iPol%8)*4))&0xf;
		int t = tuv[iMono-1].x[0]+ (unit&0x1)     - ((unit>>3)&0x1);
		int u = tuv[iMono-1].x[1]+((unit>>1)&0x1) - ((unit>>3)&0x1);
		int v = tuv[iMono-1].x[2]+((unit>>2)&0x1) - ((unit>>3)&0x1);
		
		int x = X(t,u,v);
		int y = Y(t,u,v);
		int z = Z(t,u,v);
		
		tuv[iMono].x[0] = t;
		tuv[iMono].x[1] = u;
		tuv[iMono].x[2] = v;
		
		xyz[iMono].x[0] = x;
		xyz[iMono].x[1] = y;
		xyz[iMono].x[2] = z;
		
// 		printf("%i %i %i\n", x,y,z);
		
	}
	pcfg->polSize = cs->polSize;
}

void PrintPCFG(PolyConfig* pcfg){
	for(int iMono=0; iMono<pcfg->polSize; iMono++){
		printf("[%i %i %i] / (%i %i %i)\n", pcfg->xyz[iMono].x[0], pcfg->xyz[iMono].x[1],pcfg->xyz[iMono].x[2], pcfg->tuv[iMono].x[0],pcfg->tuv[iMono].x[1],pcfg->tuv[iMono].x[2]);
	}
	printf("\n\n");
}

void MeasVars(CurState* cs, Result* res){
	for(int iPol=0; iPol<cs->nPol; iPol++){
		LoadPolymer(&res->pcfg, cs, iPol);
// 		if(cs->polSize==32)
// 			PrintPCFG(&res->pcfg);
		AddGenom(res);
	}
}



void AddGenom(Result* res){
	int di=20;
	int dx, dr;
	int gMax;
	PolyConfig* pcfg = &res->pcfg;
#if POL_TYPE == POL_LIN
	gMax=pcfg->polSize-1;
#else
	gMax = pcfg->polSize/2;
#endif
	for(int g=1, dg=1; g<=gMax; g+=dg){
		di = MAX(1,MAX(pcfg->polSize/100, dg/10));
		for(int i=0; i<pcfg->polSize; i+=di){
#if POL_TYPE == POL_LIN
			if(i+g >= pcfg->polSize) continue;
#endif
			int j = (i+g)%pcfg->polSize;
			
			dr=0;
			for(int k=0; k<3; k++){
				dx = pcfg->xyz[i].x[k]-pcfg->xyz[j].x[k];
// 				if(pcfg->polSize==4)
// 					printf("%i - %i\n", pcfg->xyz[i].x[k],pcfg->xyz[j].x[k]);
				
				dr += dx*dx;
			}
// 			if(pcfg->polSize==32){
// 				printf("i=%i, j=%i, g=%i, dr=%i\n", i,j,g,dr);
// 			}
			res->gen[g] += dr/2.0;
			res->genCounts[g]++;
		}
		dg = MAX(1,MIN(gMax/2-g, dg/10));
	}
}

