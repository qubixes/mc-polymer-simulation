#define IS_MAIN
#include "efpol.h"

int main(int argc, char** argv){
	cfg.tMax = 1000;
	cfg.seed = 37249131;
	cfg.density = 1.0;
	cfg.nRun=1;
	cfg.dblStep=0;
	cfg.LStart=16;
	cfg.initPolSize=4;
	cfg.polModel = SL_QUAD;
	char outFile[2000];
	
	if(argc <=8){
		printf("In testing mode, using defaults\n");
		cfg.dir = malloc(sizeof(char)*100);
		cfg.interval = 100;
		sprintf(cfg.dir, "./data/test");
	}
	else{
		cfg.seed = (uint)atol(argv[1]);
		cfg.dir= argv[2];
		cfg.density = atof(argv[3]);
		cfg.tMax = atof(argv[4]);
		cfg.interval = atof(argv[5]);
		cfg.initPolSize = atoi(argv[6]); ///Actually the final polymer size, see below.
		cfg.dblStep = atoi(argv[7]);
		cfg.LStart = atoi(argv[8]); ///Final lattice size.
		if(!strcmp(argv[9], "sl_equal"))
			cfg.polModel = SL_EQUAL;
		else if(!strcmp(argv[9], "sl_double"))
			cfg.polModel = SL_DOUBLE;
		else if(!strcmp(argv[9], "sl_quad"))
			cfg.polModel = SL_QUAD;
		else{
			printf("Unknown polymer model, aborting!\n");
			exit(193);
		}
		for(int iDouble=cfg.dblStep; iDouble>0; iDouble--){
			cfg.LStart/=2;
			cfg.initPolSize/=8;
		}
		if(cfg.LStart<2 || cfg.initPolSize<2){
			printf("Error polymer or lattice too small\n");
			return 212;
		}
	}
	ConfigInit();
	Timer timer;
	
	int polSize=cfg.initPolSize;
	int BLStart;
	for(BLStart=0; BLStart<=9; BLStart++){
		if((1<<BLStart)>=cfg.LStart)
			break;
	}
	int LStart=(1<<BLStart);
	int nPol=MAX(GetNPol(cfg.density, polSize, LStart*LStart*LStart),1);
	for(int L=LStart, i=0; i<=cfg.dblStep; L*=2, i++){
		CSInit(cs+i, BLStart+i, polSize, nPol, cfg.dir);
		ResultInit(&cfg, cs+i, res+i);
		polSize *= 8;
	}
	SetTransTable(&cs[0].con, cfg.polModel);
	CreatePolymers(&cfg, &cs[0]);
	
	if(CheckIntegrity(cs)) return 192;
	
	CurState csBackup;
	CSInit(&csBackup, BLStart, cfg.initPolSize, nPol, NULL);
	
// 	DoMCStep(&cs[0], 10000);
	CopyState(&cs[0], &csBackup);
	
	long totStep=0, curStep=0;
	for(int iStep=0; iStep<=cfg.dblStep; iStep++){
		totStep += cfg.nRun*cfg.tMax*cs[iStep].polSize*cs[iStep].nPol;
	}
	
	
	for(int iRun=0; iRun<cfg.nRun; iRun++){
		CopyState(&csBackup, &cs[0]);
		if(CheckIntegrity(cs)){
			printf("Error after copy\n");
			exit(0);
		}
		for(int iStep=0; iStep<=cfg.dblStep; iStep++){
			for(long t=0; t<cfg.tMax; t+= cfg.interval){
				sprintf(outFile, "%s/t=%li_dev=%i.res", cs[iStep].dir, t, 0);
				WriteLatticeFile(cs+iStep, outFile);
				TimerStart(&timer);
				DoMCStep(cs+iStep, cfg.interval);
				curStep += cfg.interval*cs[iStep].polSize*cs[iStep].nPol;
				printf("\b \b\b \b\rrun = %i/%i, L=%i, [%.1lf%%] at %.2lf M/s", iRun+1, cfg.nRun, cs[iStep].con.L, 100*curStep/(double)totStep, 1e-6*cfg.interval*cs[iStep].polSize*cs[iStep].nPol/TimerElapsed(&timer));
				fflush(NULL);
			}
			sprintf(outFile, "%s/t=%li_dev=%i.res", cs[iStep].dir, ((cfg.tMax-1)/cfg.interval+1)*cfg.interval, 0);
			WriteLatticeFile(cs+iStep, outFile);
			CheckIntegrity(cs+iStep);
			if(iStep<cfg.dblStep){
				DoubleLattice(cs+iStep, cs+iStep+1);
			}
		}
	}
	
	sleep(1);
	for(int i=0; i<=cfg.dblStep; i++)
		WriteSimulationSettings(&cfg, cs+i, i);
	
	if(CheckIntegrity(cs)) return 192;
	printf("\nFinal sanity check succeeded\n");
	return 0;
}

