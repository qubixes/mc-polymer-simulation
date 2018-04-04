#include "lowm_modes.h"

int main(int argc, char** argv){
	SimProperties sp;
	char* dir;
	FILE* pStdout;
	Timer tCompute, tIO;
	PolyTimeLapse ptl;
	
	if(argc<2){ 
		printf("Need basedir\n");
		exit(0);
	}
	dir = argv[1];
	
	if(argc>=3){
		sp.neFile = argv[2];
	}
	else
		sp.neFile = NULL;
	
	pStdout = fopen("/dev/stdout", "w");
	InitRelPos();
	sprintf(sampleDir, "%s", dir);
	if(GetNUpdates(&sp, sampleDir)==0) return 0;
	SetSimProps(&sp, sampleDir);
	
	if(sp.dT < 1000){
		printf("Thermalization: %i x %li, tau=%lf\n", sp.nTherm, sp.dT, TRelax(&sp)/1e3);
	}
	else if(sp.dT <1000000){
		printf("Thermalization: %i x %liK, tau=%lfK\n", sp.nTherm, sp.dT/1000, TRelax(&sp)/1e3);
	}
	else{
		printf("Thermalization: %i x %liM, tau=%lfM\n", sp.nTherm, sp.dT/1000000, TRelax(&sp)/1e6);
	}
	PTLAllocate(&sp, &ptl);
	sprintf(sp.resDir, "%s", sp.sampleDir);
	
	for(int iPol=0; iPol<sp.nPol; iPol++){
		for(int iDev=0; iDev<sp.nDev; iDev++){
			fprintf(pStdout, "\rPolymer %i/%i [%.1f%%]", iPol*sp.nDev+iDev+1, sp.nDev*sp.nPol, 100.0*(iPol*sp.nDev+iDev+1)/(double)(sp.nDev*sp.nPol)); 
			fflush(pStdout);
			TimerStart(&tIO);
			LoadPTL(&sp, &ptl, iPol, iDev);
			if(!sp.equalLengths || (iPol==0 && iDev==0))
				PTLInit(&sp, &ptl);
			fprintf(pStdout, "[IO: %.2f ms] ", 1e3*TimerElapsed(&tIO));
			fflush(pStdout);
			TimerStart(&tCompute);
			AddAverages(&sp, &ptl);
			fprintf(pStdout, "[Compute: %.2f ms]", 1e3*TimerElapsed(&tCompute));
			fflush(pStdout);
		}
		if(!sp.equalLengths){
			sprintf(sp.resDir, "%s/pol_%i", sp.sampleDir, iPol);
			char exec[20000];
			sprintf(exec, "mkdir -p %s", sp.resDir);
			system(exec);
			WriteAllFiles(&sp, &ptl);
		}
	}
	if(sp.equalLengths)
		WriteAllFiles(&sp, &ptl);
	printf("\n");
	return 0;
}


