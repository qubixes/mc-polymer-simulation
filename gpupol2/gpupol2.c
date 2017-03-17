#include <stdio.h>
#include <sys/time.h>
#define IS_MAIN

#include "oclpol.h"
#include "gpupol.h"


SimProperties sp;
GPULibContext gpuContext;
SimState ss[N_MAX_GPU];
GPUDeviceState devices[N_MAX_GPU];

int main(int argc, char** argv){
	struct timeval sTime, eTime;
	char fileOut[1000];
	FILE* pFile;
	
	SimPropDefaultInit(&sp);
	if(argc> 1) ///Length of a polymer
		sp.polLength = atoi(argv[1]);
	if(argc>2) ///Total Simulation time
		sp.tMax = atof(argv[2]);
	if(argc>3) ///Seed
		sp.seed = atoi(argv[3]);
	if(argc>4){ ///Output directory
		sp.dirOut = argv[4];
	}
	else sp.dirOut = NULL;
	if(argc>5)  sp.density = atof(argv[5]);
	if(argc>6)  sp.fastEq = atoi(argv[6]);
	if(argc>7)  sp.writeInterval = (long)atof(argv[7]);
	if(argc>8)  sp.LT = atoi(argv[8]);
	if(argc>9)  sp.LU = atoi(argv[9]);
	if(argc>10) sp.LV = atoi(argv[10]);
	if(argc>11) sp.equilibrated = atoi(argv[11]);
	if(argc>12) sp.double_step = atoi(argv[12]);
	
// 	printf("Box: %ix%ix%i\n", sp.LT, sp.LU, sp.LV);
// 	GPULibInit(&sp, ss, devices, &gpuContext);
	
	LoadPolymers(&sp, ss);
	GPULibLoadBuffers(&sp, ss, devices, &gpuContext);
	///Find out if there are files in the directory, if so what is the longest time one.
	
	pFile = fopen("/dev/stdout", "w");
	char cleanStr[100];
	sprintf(cleanStr, "\b \b\b \b\b \b\b \r");
	double totStep = sp.tMax-sp.curT;
	for(;sp.curT < sp.tMax; sp.curT += sp.writeInterval){
		gettimeofday(&sTime, NULL);
		if(sp.fastEq) RedistribSL(ss, &sp);
		GPULibRun(&sp, ss, devices, &gpuContext, sp.writeInterval);
		gettimeofday(&eTime, NULL);
		double totT = eTime.tv_sec-sTime.tv_sec + (eTime.tv_usec-sTime.tv_usec)/1000000.0;
		fprintf(pFile, "%s[%.1lf%%] Time elapsed = %.2f s at %.2f MF/s", cleanStr, 100*(1-(sp.tMax-sp.curT-sp.writeInterval)/totStep), totT, sp.writeInterval*sp.L*sp.nDevices/totT/1000000.0);
		fflush(pFile);
		for(int i=0; i<sp.nDevices; i++){
			sprintf(fileOut,"%s/t=%li_dev=%i.res", sp.dirOut, sp.curT+sp.writeInterval, i);
			WriteLatticeFile(&sp, ss+i, fileOut);
		}
	}
	fprintf(pFile, "\nSimulation complete\n");
	fclose(pFile);
	WriteSimulationSettings(&sp, ss);
	GPULibRelease(&sp, devices, &gpuContext);
	return 0;
}

