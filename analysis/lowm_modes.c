#include "lowm_modes.h"

double TRelax(int polSize, int polType){
	double tRouse, tInter, tRept, tau;
	
	if(polType == POL_TYPE_RING){
		tRouse = 0.24  *pow(polSize, 2.2);
		tInter = 0.042 *pow(polSize, 2.5);
		tRept  = 1.5e-3*pow(polSize, 2.97);
	}
	else{
		tRouse = 1.1    *pow(polSize, 2.2);
		tInter = 0.3    *pow(polSize, 2.5);
		tRept  = 1.26e-2*pow(polSize, 3.1);
	}
	
	tau = MAX(15e3, tRouse);
	tau = MAX(tau, tInter);
	tau = MAX(tau, tRept);
	return tau;
}

double TRelaxStretched(int polSize, int polType, double nTau){
	double tRouse, tInter, tRept, tau;
	
	if(polType == POL_TYPE_RING){
		tRouse = pow(nTau,1./0.85)*exp(-1.40/0.85)*pow(polSize, 1.88/0.85);
		tInter = tRouse;
		tRept  = pow(nTau,1./0.61)*exp(-5.02/0.61)*pow(polSize, 1.88/0.61);
	}
	else{
		tRouse = 1.1    *pow(polSize, 2.2);
		tInter = 0.3    *pow(polSize, 2.5);
		tRept  = 1.26e-2*pow(polSize, 3.1);
	}
	
	tau = MAX(15e3, tRouse);
	tau = MAX(tau, tInter);
	tau = MAX(tau, tRept);
	return tau;
}

int main(int argc, char** argv){
	SimProperties sp;
	char* dir;
	FILE* pStdout;
	Timer tCompute, tIO;
	PolyTimeLapse ptl;
// 	int TERM;
	
	if(argc<2){ 
		printf("Need basedir\n");
		exit(0);
	}
	
	dir = argv[1];
	
	pStdout = fopen("/dev/stdout", "w");
	InitRelPos();
	sprintf(sampleDir, "%s", dir);
	if(GetNUpdates(&sp, sampleDir)==0) return 0;
	SetSimProps(&sp, sampleDir);
	if(!sp.equilibrated)
		ptl.nTherm = (TRelaxStretched(sp.polSize, sp.polType,5))/sp.dT;
	else
		ptl.nTherm = 0;
	if(sp.dT < 1000){
		printf("Thermalization: %li x %li, tau=%lf\n", ptl.nTherm, sp.dT, TRelaxStretched(sp.polSize, sp.polType,5));
	}
	else if(sp.dT <1000000){
		printf("Thermalization: %li x %liK, tau=%lfK\n", ptl.nTherm, sp.dT/1000, TRelaxStretched(sp.polSize, sp.polType,5)/1e3);
	}
	else{
		printf("Thermalization: %li x %liM, tau=%lfM\n", ptl.nTherm, sp.dT/1000000, TRelaxStretched(sp.polSize, sp.polType,5)/1e6);
	}
	if(ptl.nTherm>=sp.nTime){
		printf("Samples not thermalized, continuing...\n");
		return 0;
	}
	InitArrays(&sp, &ptl);
	
	for(int iDev=0; iDev<sp.nDev; iDev++){
		for(int i=0; i<sp.nPol; i++){
			fprintf(pStdout, "\rPolymer %i/%i [%.1f%%]", i+iDev*sp.nPol+1, sp.nDev*sp.nPol, 100.0*(i+iDev*sp.nPol+1)/(double)(sp.nDev*sp.nPol)); 
			fflush(pStdout);
			TimerStart(&tIO);
			LoadPTL(&sp, &ptl, i, iDev);
			fprintf(pStdout, "[IO: %.2f ms] ", 1e3*TimerElapsed(&tIO));
			fflush(pStdout);
			TimerStart(&tCompute);
			AddAverages(&sp, &ptl);
			fprintf(pStdout, "[Compute: %.2f ms]", 1e3*TimerElapsed(&tCompute));
			fflush(pStdout);
		}
	}
	printf("\n");
	WriteAllFiles(&sp, &ptl);
// 		DestrArrays(&sp);
// 	}
// 	SpacDifDestr(&sd);
	return 0;
}


