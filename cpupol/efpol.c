#define IS_MAIN
#include "efpol.h"

int main(int argc, char** argv){
// 	cfg.polSize = 35;
// 	cfg.nTime = 200;
	cfg.tMax = 1000;
// 	cfg.interval = 1000;
// 	cfg.nEq = 0;
	cfg.seed = 37249131;
	cfg.density = 1.0;
// 	cfg.polSize=4;
// 	cfg.eqInterval = 100;
// 	cfg.efPolId=-1;
// 	cfg.EField=0.025;
	cfg.nRun=10;
	long nStep=0;
// 	char outFile[1000];
// 	char inFile[1000];
	
	
	if(argc <=3){
		printf("In testing mode, using defaults\n");
		cfg.dir = malloc(sizeof(char)*100);
		cfg.interval = 100;
		sprintf(cfg.dir, "./data/test");
	}
	else{
// 		cfg.polSize = atoi(argv[1]);
// 		cfg.nTime = (long)atof(argv[2]);
		cfg.seed = (uint)atol(argv[1]);
		cfg.dir= argv[2];
		cfg.density = atof(argv[3]);
// 		cfg.interval = (long)atof(argv[6]);
// 		cfg.EField = atof(argv[7]);
	}
	ConfigInit();
	Timer tAll, t;
	nStep += cfg.totT;
	
	int polSize=4;
	int doubleSteps=4;
	int BLStart=3;
	int LStart=(1<<BLStart);
	int nPol=MAX(GetNPol(cfg.density, polSize, LStart*LStart*LStart),1);
	for(int L=LStart, i=0; i<doubleSteps; L*=2, i++){
		CSInit(cs+i, BLStart+i, polSize, nPol);
		ResultInit(&cfg, cs+i, res+i);
		polSize *= 8;
	}
	SetTransTable(&cs[0].con);

	TimerStart(&tAll);

	CreatePolymers(&cfg, &cs[0]);
	
	if(CheckIntegrity(cs)) return 192;
// 	TimerStart(&t);
	
	CurState csBackup;
	CSInit(&csBackup, BLStart, polSize, nPol);
	
	DoMCStep(&cs[0], 10000);
	CopyState(&cs[0], &csBackup);

	for(int iRun=0; iRun<cfg.nRun; iRun++){
		CopyState(&csBackup, &cs[0]);
		if(CheckIntegrity(cs)){
			printf("Error after copy\n");
			exit(0);
		}
		for(int iStep=0; iStep<doubleSteps; iStep++){
			TimerStart(&t);
			DoMCStep(cs+iStep, cfg.tMax);
			CheckIntegrity(cs+iStep);
// 			TimerEnd(&t);
			printf("\b \b\rrun = %i/%i, L=%i, [%.1lf%%] at %.2lf M/s", iRun+1, cfg.nRun, cs[iStep].con.L, 100*(iRun*doubleSteps+iStep)/(double)(cfg.nRun*doubleSteps), 1e-6*cfg.tMax*cs[iStep].con.LAT_SIZE/TimerElapsed(&t));
			fflush(NULL);
			MeasVars(cs+iStep, res+iStep);
			if(iStep<doubleSteps-1){
				CheckIntegrity(cs+iStep);
				DoubleLattice(cs+iStep, cs+iStep+1);
				if(CheckIntegrity(cs+iStep)){
					printf("Error after doubling (source)\n");
					exit(0);
				}
				if(CheckIntegrity(cs+iStep+1)){
					printf("Error after doubling (destination)\n");
					exit(0);
				}
			}
		}
	}
	
	for(int i=0; i<doubleSteps; i++)
		WriteResults(&cfg, &cs[i], &res[i]);
	
// 	printf("\n\nTime elapsed: %.2lf s at %.2lf M/s\n", TimerElapsed(&tAll), 1e-6*nStep/TimerElapsed(&tAll));
	if(CheckIntegrity(cs)) return 192;
	printf("Final sanity check succeeded\n");
	return 0;
}

