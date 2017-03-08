#include "lowm_modes.h"

double RSQ(Coor c){
	return (c.x*c.x+c.y*c.y+c.z*c.z);
}

double DRXYZ(double xyz[3], double xyz2[3]){
	double dr=0;
	for(int k=0; k<3; k++)
		dr += (xyz[k]-xyz2[k])*(xyz[k]-xyz2[k]);
	return dr;
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

int TUV2Coor(int tuv[3]){
	return (tuv[0]+tuv[1]*LT+tuv[2]*LT*LU);
}

void Coor2TUV(int coor, int tuv[3]){
	tuv[0] = coor%LT;
	tuv[1] = (coor/LT)%LU;
	tuv[2] = coor/(LT*LU);
}

int AddCoor2TUV(int coor1, int tuv[3]){
	int tuv2[3];
	Coor2TUV(coor1, tuv2);
	for(int i=0; i<3; i++){
		tuv2[i]  = (tuv2[i]+tuv[i]+LT)%LT;
	}
	return TUV2Coor(tuv2);
}

double ShortestDistanceSQ(int coor1[3], int coor2[3]){
	
	int xyz1[3], xyz2[3];
	int newCoor[3];
	int minRSq=2*MAX_LT*MAX_LU*MAX_LV;
	
	for(int i=0; i<3; i++){
		while(coor1[i]<0)   coor1[i] += LT;
		while(coor1[i]>=LT) coor1[i] -= LT;
		while(coor2[i]<0)   coor2[i] += LT;
		while(coor2[i]>=LT) coor2[i] -= LT;
	}
	
	TUV2XYZ(coor2, xyz2);
	for(int tAdd=-LT; tAdd<=LT; tAdd += LT){
		newCoor[0] = coor1[0]+tAdd;
		for(int uAdd=-LU; uAdd<=LT; uAdd += LU){
			newCoor[1] = coor1[1]+uAdd;
			for(int vAdd=-LV; vAdd<=LT; vAdd += LV){
				newCoor[2] = coor1[2]+vAdd;
				
				TUV2XYZ(newCoor, xyz1);
				int newRsq=0; 
				for(int i=0; i<3;i++) newRsq += (xyz1[i]-xyz2[i])*(xyz1[i]-xyz2[i]);
				minRSq = MIN(newRsq, minRSq);
			}
		}
	}
	if(minRSq == 2*MAX_LT*MAX_LU*MAX_LV){
		printf("Error, no shortest distance found\n");
		exit(192);
	}
	
	return minRSq/(double)2;
}

int CompareIDouble(const void* id1, const void* id2){
	IDouble *i1, *i2;
	i1 = (IDouble*)id1; i2 = (IDouble*)id2;
	if(i1->val <  i2->val) return -1;
	if(i1->val == i2->val) return  0;
	if(i1->val >  i2->val) return  1;
	return 0;
}

void RemPeriod(int coor[3], int coorRet[3]){
	
	for(int i=0; i<3; i++){
		coorRet[i] = coor[i];
		while(coorRet[i]<0  ) coorRet[i] += LT;
		while(coorRet[i]>=LT) coorRet[i] -= LT;
	}
}

int IsValid(int a){
	return (a!=0 && a!=0x3 && a!=0xc && a!=0xf);
}

