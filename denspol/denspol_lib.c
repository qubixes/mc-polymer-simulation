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

int AddUnitToCoor(int unit, int coor, int L){
	int dw = unit>>3;
	int dt = (unit&0x1) - dw;
	int du = ((unit>>1)&0x1) - dw;
	int dv = ((unit>>2)&0x1) - dw;
	
	int t = TCoor(coor, L);
	int u = UCoor(coor, L);
	int v = VCoor(coor, L);
	
	return TUV2Coor((t+dt+L)%L, (u+du+L)%L, (v+dv+L)%L, L);
}

void TUV2XYZ(int tuv[3], int xyz[3]){
	xyz[0] = tuv[0]+tuv[1]-tuv[2];
	xyz[1] = tuv[0]-tuv[1]       ;
	xyz[2] =               tuv[2];
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