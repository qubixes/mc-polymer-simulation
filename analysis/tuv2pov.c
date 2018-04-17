#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#define IS_MAIN
#include "lowm_modes.h"
#include "rng.h"
#define ALL_POLYMERS -1

// #define PI 3.141592653589793

typedef struct Coor3{
	double x, y, z;
}Coor3;

typedef struct lattice{
	int*** tuvLattice;
	int minTUV[3];
}Lattice;

typedef struct Data{
	int*** tuv;
	Coor3** pos;
	int* nMono;
	int* polTypes;
	double (*TUV2Distance) (double*, double*, int);
	int nPol;
	int maxNMono;
	int nTotMono;
	int L;
	Lattice lat;
}Data;




typedef struct Camera{
	Coor3 camPos;
	Coor3 pointAt;
	Coor3 lightPos;
	Coor3 lightAreaVec1;
	Coor3 lightAreaVec2;
	double lightAngle;
	double angle;
}Camera;
	
Data* ReadData(char* file);
int WriteToFile(FILE* pFile, int length);
int WriteToPovFile(Data* dat, char* file, int polId);
Data* ReadPolymerFromFile(char* file, int column);
int TrashLine(FILE *pFile);
int GetNColumns(FILE* pFile);
void FWriteMono(FILE* pFile, Coor3 coor, char* color, double size);
void FWriteLink(FILE* pFile, Coor3 coor1, Coor3 coor2, char* color, double size);

int main(int argc, char** argv){
	char* tuvFile;
	char* povFile;
	int polId = ALL_POLYMERS;
	
	if(argc < 3){
		printf("Bad number of arguments\n");
		return 192;
	}
	
	tuvFile = argv[1];
	povFile = argv[2];
	if(argc >= 4)
		polId = atoi(argv[3]);
	
	Data* data = ReadData(tuvFile);
	
	WriteToPovFile(data, povFile, polId);
	return 0;	
}


void TUV2DblXYZ(int tuv[3], Coor3* coor){
	coor->x = (tuv[0]+tuv[1]-tuv[2])/sqrt(2);
	coor->y = (tuv[0]-tuv[1]       )/sqrt(2);
	coor->z = (              tuv[2])/sqrt(2);
}

Data* NewData(int maxNMono, int nPol, int L, int boundaryCond){
	Data* data    = malloc(sizeof(Data));
	data->nPol    = nPol;
	data->maxNMono= maxNMono;
	data->L       = L;
	
	data->tuv = malloc(sizeof(int**)*nPol);
	data->pos = malloc(sizeof(Coor3*)*nPol);
	for(int i=0; i<nPol; i++){
		data->tuv[i] = malloc(sizeof(int*)*maxNMono);
		data->pos[i] = malloc(sizeof(Coor3)*maxNMono);
		for(int j=0; j<maxNMono; j++)
			data->tuv[i][j] = malloc(sizeof(int)*3);
	}
	
	data->nMono    = malloc(sizeof(int)*nPol);
	data->polTypes = malloc(sizeof(int)*nPol);
	
	return data;
}

void LatticeInit(Data* data){
	Lattice* lat = &data->lat;
	int maxTUV[] = {-9999999,-9999999,-9999999};
	for(int k=0; k<3; k++){
		lat->minTUV[k] = 9999999;
	}
	for(int iPol=0; iPol<data->nPol; iPol++){
		for(int iMono=0; iMono<data->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++){
				maxTUV[k] = MAX(data->tuv[iPol][iMono][k], maxTUV[k]);
				lat->minTUV[k] = MIN(data->tuv[iPol][iMono][k], lat->minTUV[k]);
			}
		}
	}
	
// 	printf("dtuv=(%i,%i,%i)\n", maxTUV[0]-lat->minTUV[0]+1, maxTUV[1]-lat->minTUV[1]+1, maxTUV[2]-lat->minTUV[2]+1);
	
	lat->tuvLattice = malloc(sizeof(int**)*(maxTUV[0]-lat->minTUV[0]+1));
	for(int tr=0; tr<maxTUV[0]-lat->minTUV[0]+1; tr++){
		lat->tuvLattice[tr] = malloc(sizeof(int*)*(maxTUV[1]-lat->minTUV[1]+1));
		for(int ur=0; ur<maxTUV[1]-lat->minTUV[1]+1; ur++){
			lat->tuvLattice[tr][ur] = malloc(sizeof(int)*(maxTUV[2]-lat->minTUV[2]+1));
			for(int vr=0; vr<maxTUV[2]-lat->minTUV[2]+1; vr++){
				lat->tuvLattice[tr][ur][vr] = 0;
			}
		}
	}
}

void AddCharToTUV(char c, int tuv[3]){
	int bond = CharToHex(c);
	
	for(int i=0; i<3; i++){
		tuv[i] += (bond>>i)&0x1;
		tuv[i] -= (bond>>3)&0x1;
	}
}

Data* ReadData(char* file){
	int nPol, maxNMono, L, nOrigMono;
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int i=0; i<3; i++)
		fscanf(pFile, "%*s %i", &L);
	
	fscanf(pFile, "%*s %i", &nPol);
	
	fscanf(pFile, "%*s %i", &maxNMono);
	Data* data = NewData(maxNMono, nPol, L, 0);
	char* str  = malloc(sizeof(char)*(maxNMono+1));
	
	data->nTotMono=0;
	for(int iPol=0; iPol<nPol; iPol++){
		int tuv[3];
		fscanf(pFile, "%*s %i %i %i %i %s", &nOrigMono, tuv, tuv+1, tuv+2, str);
		data->nMono[iPol]=0;
		for(int iMono=0; iMono<nOrigMono; iMono++, data->nMono[iPol]++){
			if(iMono && str[iMono-1] == '0'){
				AddCharToTUV(str[iMono], tuv);
				data->nMono[iPol]--;
				continue;
			}
			for(int k=0; k<3; k++)
				data->tuv[iPol][data->nMono[iPol]][k] = tuv[k];
			TUV2DblXYZ(data->tuv[iPol][data->nMono[iPol]], &data->pos[iPol][data->nMono[iPol]]);
// 			printf("(%i %i %i)\n", tuv[0], tuv[1], tuv[2]); 
			AddCharToTUV(str[iMono], tuv);
		}
// 		printf("\n\n");
		data->nTotMono += data->nMono[iPol];
		if(str[data->nMono[iPol]] == 'f') data->polTypes[iPol] = POL_TYPE_LIN;
		else data->polTypes[iPol] = POL_TYPE_RING;
	}
	fclose(pFile);
	return data;
}


double ColorFrac(double frac){
	frac -= (int)frac;
	if(frac<0.25)
		return 4*frac;
	else if(frac<0.5)
		return 1;
	else if(frac<0.75)
		return 3-4*frac;
	else
		return 0;
}
	
// void GetColor(double frac, char* color){
// 	sprintf(color, "rgbf<%lf, %lf, %lf>", ColorFrac(frac), ColorFrac(frac+0.25), ColorFrac(frac+0.5));
// }

// void GetMonoColor(int charge, char* color){
// 	sprintf(color, "rgbf<%lf, %lf, %lf>", (double)(charge+1)/2,(double)(charge+1)/2,(double)(charge+1)/2);
// }

void GetSinglePolColor(Data* data, int iPol, int iMono, char* color){
	double rgbf[3];
	
	rgbf[0] = ColorFrac(iMono/(double)(data->nMono[iPol])     );
	rgbf[1] = ColorFrac(iMono/(double)(data->nMono[iPol])+0.25);
	rgbf[2] = ColorFrac(iMono/(double)(data->nMono[iPol])+0.5 );
	
	sprintf(color, "rgbf<%lf, %lf, %lf>", rgbf[0], rgbf[1], rgbf[2]);
}


void GetMonoColor(Data* data, int iPol, int iMono, char* color){
	double rgbf[3];
	double maxWhite=0.3;
	double monoR=iMono/(double)(data->nMono[iPol]);
	
	rgbf[0] = ColorFrac(iPol/(double)(data->nPol));
	rgbf[1] = ColorFrac(iPol/(double)(data->nPol)+0.25);
	rgbf[2] = ColorFrac(iPol/(double)(data->nPol)+0.5);
	
	for(int k=0; k<3; k++){
		double white = maxWhite*monoR;
		rgbf[k] = white+(1-white)*rgbf[k];
	}
	sprintf(color, "rgbf<%lf, %lf, %lf>", rgbf[0], rgbf[1], rgbf[2]);
}

Coor3 Min3(Coor3 a, Coor3 b){
	Coor3 res;
	
	res.x = a.x-b.x;
	res.y = a.y-b.y;
	res.z = a.z-b.z;
	return res;
}

Coor3 Add3(Coor3 a, Coor3 b){
	Coor3 res;
	
	res.x = a.x+b.x;
	res.y = a.y+b.y;
	res.z = a.z+b.z;
	return res;
}

Coor3 ScalMul3(double scal, Coor3 c){
	Coor3 res;
	
	res.x = c.x*scal;
	res.y = c.y*scal;
	res.z = c.z*scal;
	return res;
}

Coor3 AXPY(double a, Coor3 v1, Coor3 v2){
	Coor3 res;
	
	res.x = a*v1.x+v2.x;
	res.y = a*v1.y+v2.y;
	res.z = a*v1.z+v2.z;
	return res;
}

double Dot3(Coor3 a, Coor3 b){
	double res=0;
	
	res += a.x*b.x;
	res += a.y*b.y;
	res += a.z*b.z;
	return res;
}

Coor3 Normalize(Coor3 c){
	Coor3 cp;
	double norm = sqrt(Dot3(c,c));
	cp.x = c.x/norm;
	cp.y = c.y/norm;
	cp.z = c.z/norm;
	return cp;
}

void PrintCoor(Coor3 c){
	printf("(%lf %lf %lf)", c.x, c.y, c.z);
}

Coor3 GetOrthogonalVec(Coor3 unit){
	unsigned int rng[4];
	Coor3 c;
	Seed(rng, 1029471);
	do{
		c.x = 2*DRng(rng)-1;
		c.y = 2*DRng(rng)-1;
		c.z = 2*DRng(rng)-1;
		c = Normalize(c);
	} while (fabs(Dot3(c,unit)) > 0.5);
	
	c = AXPY(-Dot3(c,unit), unit, c);
	return Normalize(c);
}

void GetDoubleOrthogonal(Coor3 unit, Coor3* orth1, Coor3* orth2){
	unsigned int rng[4];
	
	*orth1 = GetOrthogonalVec(unit);
	
	Coor3 c;
	Seed(rng, 1029471);
	do{
		c.x = 2*DRng(rng)-1;
		c.y = 2*DRng(rng)-1;
		c.z = 2*DRng(rng)-1;
// 		PrintCoor(c); printf("\n");
		c = Normalize(c);
// 		printf("%lf %lf\n", Dot3(c,unit), Dot3(c, *orth1));
	} while (fabs(Dot3(c,unit)) > 0.5 || fabs(Dot3(c, *orth1)) > 0.5 );
	
	c = AXPY(-Dot3(c,unit), unit, c);
	*orth2 = AXPY(-Dot3(c,*orth1), *orth1, c);
	*orth2 = Normalize(*orth2);
}

/** angle in degrees. **/

Coor3 RotateWithAxis(Coor3 axis, double dAngle, Coor3 unit){
	Coor3 nUnit = {0,0,0};
	double angle = dAngle*2*PI/360;
	
	double c = cos(angle);
	double s = sin(angle);
	double C = 1-c;
	double x = axis.x;
	double y = axis.y;
	double z = axis.z;
	
	nUnit.x += unit.x*(x*x*C+c);
	nUnit.x += unit.y*(x*y*C-z*s);
	nUnit.x += unit.z*(x*z*C+y*s);
	
	nUnit.y += unit.x*(y*x*C+z*s);
	nUnit.y += unit.y*(y*y*C+c);
	nUnit.y += unit.z*(y*z*C-x*s);
	
	nUnit.z += unit.x*(z*x*C-y*s);
	nUnit.z += unit.y*(z*y*C+x*s);
	nUnit.z += unit.z*(z*z*C+c);
	
	return nUnit;
}

Camera GetCamera(Data* data, int polId){
	double cosFac;
	Camera cam;
	Coor3 unitVec = {0.72347, 0.502397, 0.224491};
	Coor3 yVec, z;
	double maxDist = 0.0;
	
	unitVec = Normalize(unitVec);
	
	cam.pointAt.x=0; cam.pointAt.y=0; cam.pointAt.z=0;
	
	///First get position to look at:
	
	long nAvg=0;
	for(int iPol=0; iPol<data->nPol; iPol++){
		if(polId != ALL_POLYMERS && iPol != polId) continue;
		for(int iMono=0; iMono<data->nMono[iPol]; iMono++){
			cam.pointAt = Add3(cam.pointAt, data->pos[iPol][iMono]);
			nAvg++;
		}
	}
	cam.pointAt.x /= nAvg;
	cam.pointAt.y /= nAvg;
	cam.pointAt.z /= nAvg;
	
	cam.angle = 15;
	cosFac = tan((cam.angle/2.0-0.5)/360.0*2*PI);
	
	for(int iPol=0; iPol<data->nPol; iPol++){
		if(polId != ALL_POLYMERS && iPol != polId) continue;
		for(int iMono=0; iMono<data->nMono[iPol]; iMono++){
			z = Min3(data->pos[iPol][iMono], cam.pointAt);
			double dx = Dot3(z, unitVec);
		
			if(dx<0)
				dx = -sqrt(-dx);
			else
				dx = sqrt(dx);
		
			yVec = Min3(z, ScalMul3(dx, unitVec));
			double y = sqrt(Dot3(yVec, yVec));
			double x = y/cosFac;
			if(x+dx>maxDist)
				maxDist = x+dx;
		}
	}
	cam.camPos = AXPY((maxDist), unitVec, cam.pointAt);
	cam.angle += 0.5;
	cam.lightAngle=40;
	
	Coor3 orth = GetOrthogonalVec(unitVec);
	Coor3 lightUnit = RotateWithAxis(orth, cam.lightAngle, unitVec);
	
	cam.lightPos = AXPY(maxDist, lightUnit, cam.pointAt);
	
	Coor3 orth1, orth2;
	GetDoubleOrthogonal(lightUnit, &orth1, &orth2);
	
	cam.lightAreaVec1 = ScalMul3(0.1*maxDist, orth1);
	cam.lightAreaVec2 = ScalMul3(0.1*maxDist, orth2);
	
// 	PrintCoor(lightUnit); printf(" "); PrintCoor(orth1); printf(" "); PrintCoor(orth2); printf("\n");
	
// 	printf("%lf %lf %lf\n", Dot3(lightUnit, orth1), Dot3(lightUnit, orth2), Dot3(orth1,orth2));
	
	return cam;
}

void FWritePos(FILE* pFile, Coor3 coor){
	fprintf(pFile, " <%lf, %lf, %lf> ", coor.x, coor.y, coor.z);
}

int WriteToPovFile(Data* data, char* file, int polId){
	FILE* pFile;
	char color[1000];
	Camera cam;
	LatticeInit(data);
	double size=0.15;
	
	cam = GetCamera(data, polId);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(0);
	}
	
	fprintf(pFile, "light_source { <%lf, %lf, %lf> ", cam.lightPos.x, cam.lightPos.y, cam.lightPos.z);
	fprintf(pFile, "color White area_light <%lf,%lf,%lf>,", cam.lightAreaVec1.x, cam.lightAreaVec1.y, cam.lightAreaVec1.z);
	fprintf(pFile, "<%lf,%lf,%lf>, 3, 3\n\tadaptive 5\n}\n", cam.lightAreaVec2.x, cam.lightAreaVec2.y, cam.lightAreaVec2.z);
	
	fprintf(pFile, "camera { location "); FWritePos(pFile, cam.camPos);
	fprintf(pFile, "up y\n right x\n");
	fprintf(pFile, "angle %lf \n", cam.angle);
	fprintf(pFile, "look_at "); FWritePos(pFile, cam.pointAt);
	
	fprintf(pFile, "}\n");
	for(int iPol=0; iPol<data->nPol; iPol++){
		if(polId != ALL_POLYMERS && polId != iPol) continue;
		fprintf(pFile, "\nmerge {\n");
		for(int iMono=0; iMono<data->nMono[iPol]-1; iMono++){
			if(polId != ALL_POLYMERS)
				GetSinglePolColor(data, iPol, iMono, color);
			else
				GetMonoColor(data, iPol, iMono, color);
			
			int tr = data->tuv[iPol][iMono][0]-data->lat.minTUV[0];
			int ur = data->tuv[iPol][iMono][1]-data->lat.minTUV[1];
			int vr = data->tuv[iPol][iMono][2]-data->lat.minTUV[2];
			if(!data->lat.tuvLattice[tr][ur][vr]){
				FWriteMono(pFile, data->pos[iPol][iMono], color, size);
				data->lat.tuvLattice[tr][ur][vr]=1;
			}
			FWriteLink(pFile, data->pos[iPol][iMono], data->pos[iPol][iMono+1], color, size);
		}
		int tr = data->tuv[iPol][data->nMono[iPol]-1][0]-data->lat.minTUV[0];
		int ur = data->tuv[iPol][data->nMono[iPol]-1][1]-data->lat.minTUV[1];
		int vr = data->tuv[iPol][data->nMono[iPol]-1][2]-data->lat.minTUV[2];
		if(!data->lat.tuvLattice[tr][ur][vr]){
			FWriteMono(pFile, data->pos[iPol][data->nMono[iPol]-1], color, size);
			data->lat.tuvLattice[tr][ur][vr]=1;
		}			
		fprintf(pFile, "}\n");
	}
	fclose(pFile);
	return 0;
}

void FWriteLink(FILE* pFile, Coor3 coor1, Coor3 coor2, char* color, double size){
	fprintf(pFile, "object { cylinder {");
	FWritePos(pFile, coor1);
	fprintf(pFile, ", ");
	FWritePos(pFile, coor2);
	fprintf(pFile, ", %lf texture {pigment {color %s}} finish { ambient 0.1 diffuse 0.6} } } \n", size, color);
}

void FWriteMono(FILE* pFile, Coor3 coor, char* color, double size){
	fprintf(pFile, "object { sphere { "); 
	FWritePos(pFile, coor);
	fprintf(pFile, ", %lf texture {pigment {color %s}} finish {ambient 0.1 diffuse 0.6} } }\n", 1.5*size, color);
}
// 	"<%lf, %lf, %lf>, <%lf, %lf, %lf>, coor1.x, pol->coors[i].y, pol->coors[i].z, pol->coors[i+1].x, pol->coors[i+1].y, pol->coors[i+1].z, color); 

	

int TrashLine(FILE *pFile){
	char c='\0';
	int ret=1;
	while(c != '\n' && ret>0){
		ret = fscanf(pFile, "%c", &c);
// 		printf("%iv vs %i\n", c, '\n');
	}
	return 0;
}

int GetNColumns(FILE* pFile){
	int i;
	long newLinePos;
	TrashLine(pFile);
	newLinePos = ftell(pFile);
	rewind(pFile);
	
	for(i=0; ftell(pFile)<newLinePos; i++)
		fscanf(pFile, "%*s");
	printf("File contains %i columns\n", i-1);
	rewind(pFile);
	return i;
}
	