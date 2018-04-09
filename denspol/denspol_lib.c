#include "denspol_lib.h"

int ValidateAddUnitVectors(int a, int b, int* c){
	int r, valid;
	if((a|b) != 0xf && (a&b))
		return 0;
	r = (((a|b)==0xf)?(a&b):(a|b));
	valid = IsValid(r);
	*c = r;
	return valid;
}

int AddUnitVectors(int a, int b){
	return (((a|b)==0xf)?(a&b):(a|b));
}

int IsValid(int a){
	return (a!=0 && a!=0x3 && a!=0xc && a!=0xf);
}

int TUV2Coor(int t, int u, int v, int L){
	return (t+u*L+v*L*L);
}

int TCoor(int coor, int L){
	return (coor%L);
}

int UCoor(int coor, int L){
	return ((coor/L)%L);
}

int VCoor(int coor, int L){
	return (coor/(L*L));
}

int AddUnitToCoorPeriod(int unit, int coor, int L, int* OOB){
	*OOB = 0; /// Never out of bounds
	int dw =   unit>>3;
	int dt = ( unit    &0x1) - dw;
	int du = ((unit>>1)&0x1) - dw;
	int dv = ((unit>>2)&0x1) - dw;
	
	int t = TCoor(coor, L);
	int u = UCoor(coor, L);
	int v = VCoor(coor, L);
	
	return TUV2Coor((t+dt+L)%L, (u+du+L)%L, (v+dv+L)%L, L);
}


int AddUnitToCoorWBounds(int unit, int coor, int L, int* OOB){
	int dw =   unit>>3;
	int dt = ( unit    &0x1) - dw;
	int du = ((unit>>1)&0x1) - dw;
	int dv = ((unit>>2)&0x1) - dw;
	
	int t = TCoor(coor, L)+dt;
	int u = UCoor(coor, L)+du;
	int v = VCoor(coor, L)+dv;
	
	*OOB  = ((t+L)/L-1);
	*OOB |= ((u+L)/L-1);
	*OOB |= ((v+L)/L-1);
	
	if(*OOB) return -1;
	
	return TUV2Coor(t,u,v,L);
}

void PrintCoor(int coor, int L){
	int t = TCoor(coor, L);
	int u = UCoor(coor, L);
	int v = VCoor(coor, L);
	
	printf("(%i %i %i)", t,u,v);
}

void Coor2TUV(int coor, int tuv[3], int L){
	tuv[0] = TCoor(coor, L);
	tuv[1] = UCoor(coor, L);
	tuv[2] = VCoor(coor, L);
}

void TUV2XYZ(int tuv[3], int xyz[3]){
	xyz[0] = tuv[0]+tuv[1]-tuv[2];
	xyz[1] = tuv[0]-tuv[1]       ;
	xyz[2] =               tuv[2];
}

void DTUV2XYZ(double tuv[3], double xyz[3]){
	double invSqrt2=1./sqrt(2);
	xyz[0] = (tuv[0]+tuv[1]-tuv[2])*invSqrt2;
	xyz[1] = (tuv[0]-tuv[1]       )*invSqrt2;
	xyz[2] = (              tuv[2])*invSqrt2;
}

double DistanceSq(double xyz1[3], double xyz2[3]){
	double rsq=0;
	for(int k=0; k<3; k++)
		rsq += (xyz1[k]-xyz2[k])*(xyz1[k]-xyz2[k]);
	return (rsq/2.0);
}


double Distance(double xyz1[3], double xyz2[3]){
	return sqrt(DistanceSq(xyz1,xyz2));
}

void UnitToXYZ(int unit, int xyz[3]){
	int w = unit>>3;
	int t = (unit&0x1)      - w;
	int u = ((unit>>1)&0x1) - w;
	int v = ((unit>>2)&0x1) - w;
	
	int tuv[3] = {t,u,v};
	
	TUV2XYZ(tuv, xyz);
}

int UnitInProd(int unit1, int unit2){
	int xyz1[3];
	int xyz2[3];
	
	UnitToXYZ(unit1, xyz1);
	UnitToXYZ(unit2, xyz2);
	
	int dot=0;
	for(int k=0; k<3; k++) dot += xyz1[k]*xyz2[k];
	return dot;
}

int CharToHex(char c){
	int hex;
	
	hex = c-'0';
	if(hex>=10) hex = 10+(int)(c-'a');
	return hex;
}

double MagnificationRatio(int nMono, int nMonoOrig, int polType){
	if(polType == POL_TYPE_LIN)
		return ((nMono-1)/(double)(nMonoOrig-1));
	else
		return (nMono/(double)nMonoOrig);
}

int SetLatticeSphere(int* topo, int* bondOcc, int L, int* nBondUsed){
	int nLatticeUsed=0;
	(*nBondUsed)=0;
	double r = L/sqrt(12);
	
	double tuvMid[3] = {0.5*L, 0.5*L, 0.5*L};
	double xyzMid[3];
	DTUV2XYZ(tuvMid, xyzMid);
	
	for(int t=0; t<L; t++){
		for(int u=0; u<L; u++){
			for(int v=0; v<L; v++){
				int site = t+u*L+v*L*L;
				double newTuv[3];
				
				int inside=1;
				for(int dt=0; dt<2 && inside; dt++){
					for(int du=0; du<2 && inside; du++){
						for(int dv=0; dv<2 && inside; dv++){
							newTuv[0] = t+dt;
							newTuv[1] = u+du;
							newTuv[2] = v+dv;
							
							double dxyz[3];
							DTUV2XYZ(newTuv, dxyz);
							
							if(Distance(dxyz, xyzMid) > r)
								inside=0;
						}
					}
				}
				
				if(inside){
					topo[site] = 0;
					bondOcc[site] = 0;
					nLatticeUsed++;
				}
				else{
					topo[site]      = -1;
					bondOcc[site]   = 0xffffffff;
				}
			}
		}
	}
	
	for(int coor=0; coor<L*L*L; coor++){
		if(topo[coor]<0) continue;
		for(int unit=0x1; unit<0xf; unit++){
			if(!IsValid(unit)) continue;
			int OOB;
			int newCoor = AddUnitToCoorWBounds(unit, coor, L, &OOB);
			if(!OOB && topo[newCoor] >=0)
				(*nBondUsed)++;
		}
	}
	(*nBondUsed) /= 2;
	return nLatticeUsed;
}

int SetLatticeEmpty(int* topo, int* bondOcc, int L, int* nBondUsed){
	int LSize = L*L*L;
	
	for(int site=0; site<LSize; site++){
		topo[site] = 0;
		bondOcc[site] = 0;
	}
	*nBondUsed = LSize*6;
	return LSize;
}

int* ComputePolLengths(char* file, int nBondUsed, double density){
	FILE* pFile= fopen(file, "r");
	if(!pFile) printf("Error opening file %s for reading\n", file);
	
	int nPol;
	int nTotMonoOrig = 0;
	
	fscanf(pFile, "%*s %i", &nPol);
	int* origLengths = malloc(sizeof(int)*nPol);
	int* newLengths = malloc(sizeof(int)*nPol);
	int* polTypes = malloc(sizeof(int)*nPol);
	char* pt = malloc(sizeof(char)*300);
	
	for(int i=0; i<2; i++) fscanf(pFile, "%*s %*s");
	for(int iPol=0; iPol<nPol; iPol++){
		fscanf(pFile, "%s %i", pt, &origLengths[iPol]);
		if(!strcmp(pt, "lin"))
			polTypes[iPol] = POL_TYPE_LIN;
		else
			polTypes[iPol] = POL_TYPE_RING;
		nTotMonoOrig += origLengths[iPol];
	}
	fclose(pFile);
	
	int nTotMonoTarget = (int)((density/6.0)*nBondUsed+0.5);
	
	double leftover=0;
	for(int iPol=0; iPol<nPol; iPol++){
		double dblNewNMono = origLengths[iPol]*(nTotMonoTarget/(double)nTotMonoOrig)+leftover;
		int newNMono = (int)(dblNewNMono+0.5);
		leftover = dblNewNMono-newNMono;
		if(polTypes[iPol] == POL_TYPE_LIN)
			newLengths[iPol] = MAX(3,newNMono);
		else
			newLengths[iPol] = MAX(4,newNMono);
	}
	free(origLengths);
	free(pt); free(polTypes);
	return newLengths;
}


/// Returns 1 if the move is accepted on the grounds of the harmonic potential.
/// Otherwise returns 0.

int TestMoveHP(CurState* cs, LookupTables* lt, int iMono, Polymer* pol, int newUnit){
	int iPol = pol-cs->pol;
	int ituv[3], newituv[3];
// 	printf("iMono = %i, iPol = %i\n", iMono, iPol);
	int OOB;
	int newCoor = cs->AddUnitToCoor(newUnit, pol->coorPol[iMono], cs->L, &OOB);
	if(OOB) return 0;
	
	Coor2TUV(pol->coorPol[iMono], ituv, cs->L);
	Coor2TUV(newCoor, newituv, cs->L);
	double dE=0;
	
	for(int iInter=0; iInter<lt->hp.nInter[iPol][iMono]; iInter++){
		int jMono       = lt->hp.inter[iPol][iMono][iInter].iMono;
		int jPol        = lt->hp.inter[iPol][iMono][iInter].iPol;
		double strength = lt->hp.inter[iPol][iMono][iInter].strength;
		
		Polymer* polJ = cs->pol+jPol;
		int jtuv[3], dtuv[3];
		Coor2TUV(polJ->coorPol[jMono], jtuv, cs->L);
		
		for(int k=0; k<3; k++) dtuv[k] = jtuv[k] - ituv[k];
		double distance = lt->hp.distance[dtuv[0]][dtuv[1]][dtuv[2]];
		dE -= strength*distance;
		
		for(int k=0; k<3; k++) dtuv[k] = jtuv[k] - newituv[k];
		distance = lt->hp.distance[dtuv[0]][dtuv[1]][dtuv[2]];
		dE += strength*distance;
	}
	
	double prob = (dE<=0)?1:(exp(-dE));
	double comp = DRng(cs->rngState);
	if(comp<prob){
		return 1;
	}
	else{
// 		printf("comp=%lf, prob=%lf\n", comp, prob);
		return 0;
	}
}

