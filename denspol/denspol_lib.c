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