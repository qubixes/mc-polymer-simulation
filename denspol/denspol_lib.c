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
	xyz[0] = tuv[0]+tuv[1]-tuv[2];
	xyz[1] = tuv[0]-tuv[1]       ;
	xyz[2] =               tuv[2];
}

double Distance(double xyz1[3], double xyz2[3]){
	double rsq=0;
	for(int k=0; k<3; k++)
		rsq += (xyz1[k]-xyz2[k])*(xyz1[k]-xyz2[k]);
	return sqrt(rsq);
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

