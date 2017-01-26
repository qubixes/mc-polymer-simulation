#include "lowm_modes.h"

void InitArrays(SimProperties* sp, PolyTimeLapse* ptl){
	SpacDifInit(&sd);
	PTLInit(sp, ptl, &sd);
	allFiles = malloc(sizeof(char*)*sp->nTime);
	filePos = malloc(sizeof(fpos_t*)*sp->nTime);
	strIn = malloc(sizeof(char)*(sp->polSize+1));
	t = malloc(sizeof(int)*(sp->polSize));
	u = malloc(sizeof(int)*(sp->polSize));
	v = malloc(sizeof(int)*(sp->polSize));
}

void InitRelPos(){
	int t,u,v,w, iRel=0;
	for(int i=0; i<16; i++){
		if(!IsValid(i)) continue;
		t = i%2; u=(i/2)%2; v = (i/4)%2; w = i/8;
		t -= w; u -= w; v -= w;
		tuvRelPos[iRel][0] = t;
		tuvRelPos[iRel][1] = u;
		tuvRelPos[iRel][2] = v;
// 		printf("[%i %i %i]\n", t,u,v);
		iRel++;
	}
}

void PTLInit(SimProperties* sp, PolyTimeLapse* ptl, SpacDif* sd){
	ptl->sinfac = malloc(sizeof(double*)*(sp->polSize/2+1));
	ptl->cosfac = malloc(sizeof(double*)*(sp->polSize/2+1));
	
	for(int p=0;p<=sp->polSize/2;p++){
		ptl->sinfac[p] = malloc(sizeof(double)*sp->polSize);
		ptl->cosfac[p] = malloc(sizeof(double)*sp->polSize);
		for(int n=0;n<sp->polSize;n++){
			ptl->sinfac[p][n]=sqrt(1.0/sp->polSize)*sin((2*PI*(n+0.25)*p)/(double)sp->polSize);
			ptl->cosfac[p][n]=sqrt(1.0/sp->polSize)*cos((2*PI*(n+0.25)*p)/(double)sp->polSize);
		}
	}
	int step=1; ptl->nModes=0; ptl->modeList = malloc(sizeof(int)*sp->polSize/2);
	for(int p=0; p<sp->polSize/2; p+=step){
		ptl->modeList[ptl->nModes++] = p;
		step = MAX(1, p/5);
	}
	
	ptl->pointDensity = malloc(sizeof(int)*sp->polSize);
	for(int i=0; i<sp->polSize; i++) ptl->pointDensity[i]=0;
	
	
	int maxTUV;
	
	maxTUV = MIN(sp->polSize, 2*LT);
	
	for(int t=-maxTUV; t<maxTUV; t++){
		for(int u=-maxTUV; u<maxTUV; u++){
			for(int v=-maxTUV; v<maxTUV; v++){
				int x = t+u-v;
				int y = t-u;
				int z = v;
				
				int rSq = x*x+y*y+z*z;
				double r = sqrt(rSq/2.0);
				int bin = (int)r;
				if(bin<sp->polSize)
					ptl->pointDensity[bin]++;
			}
		}
	}
	
// 	for(int i=0; i<sp->polSize; i++)
// 		printf("%i %i %lf\n", i, ptl->pointDensity[i], ptl->pointDensity[i]/(double)(i*i));
// 	exit(0);
	
	ptl->genomList = malloc(sizeof(DInt)*sp->polSize*sp->polSize);
	ptl->nGenom = 0;
	int di=20;
	int gMax;
	if(sp->polType == POL_TYPE_LIN) gMax=sp->polSize-1;
	else gMax = sp->polSize/2;
	for(int g=1, dg=1; g<=gMax; g+=dg){
		di = MAX(1,MAX(sp->polSize/30, dg/3));
		for(int i=0; i<sp->polSize; i+=di){
			if(i+g >= sp->polSize && sp->polType == POL_TYPE_LIN) continue;
			int j = (i+g)%sp->polSize;
// 			if(j>i){
				ptl->genomList[ptl->nGenom].x = i;
				ptl->genomList[ptl->nGenom++].y = j;
// 			}
// 			else{
// 				ptl->genomList[ptl->nGenom].x = j;
// 				ptl->genomList[ptl->nGenom++].y = i;
// 			}
		}
		dg = MAX(1,MIN(gMax-g, g/10));
	}
	printf("nGenom=%i\n", ptl->nGenom);
	ptl->genomProb = malloc(sizeof(long*)*(gMax+1));
	ptl->genomR = malloc(sizeof(double*)*(gMax+1));
	ptl->genomCount = malloc(sizeof(long)*(gMax+1));
	for(int i=0; i<(gMax+1); i++){
		ptl->genomProb[i] = malloc(sizeof(long)*sp->polSize);
		ptl->genomR[i] = malloc(sizeof(double)*sp->polSize);
		ptl->genomCount[i] = 0;
		for(int j=0; j<sp->polSize; j++){
			ptl->genomProb[i][j]=0;
			ptl->genomR[i][j]=0;
		}
	}
	
	ptl->sqrtList = malloc(sizeof(double)*sp->polSize*sp->polSize*2);
	for(int i=0; i<sp->polSize*sp->polSize*2; i++)
		ptl->sqrtList[i] = sqrt(i/2.0);
	
	ptl->nEqd = sp->nTime-ptl->nTherm;
	ptl->rGyrT = malloc(sizeof(double)*sp->nTime); 
	ptl->avgUnitCor = malloc(sizeof(double)*sp->polSize);
	ptl->avgGenom = malloc(sizeof(double)*sp->polSize);
	ptl->cmsDif = malloc(sizeof(double)*ptl->nEqd);
	ptl->smDif = malloc(sizeof(double)*ptl->nEqd);
	ptl->mmDif = malloc(sizeof(double)*ptl->nEqd);
	ptl->emDif = malloc(sizeof(double)*ptl->nEqd);
	ptl->avgSL = malloc(sizeof(double)*sp->nTime);
	ptl->avgRee = malloc(sizeof(double)*ptl->nEqd);
	ptl->avgShearMod = malloc(sizeof(double)*ptl->nEqd);
	ptl->polys = malloc(sizeof(PolyConfig)*ptl->nEqd);
	ptl->avgRGyr = 0;
	ptl->pc = malloc(sizeof(double*)*sp->polSize);
	
	
	ptl->L = 400;
	ptl->lattice = (LatPoint*)malloc(sizeof(LatPoint)*ptl->L*ptl->L*ptl->L)+(ptl->L*ptl->L*ptl->L)/2;
	for(int i=0; i<ptl->L*ptl->L*ptl->L; i++)
		ptl->lattice[i-(ptl->L*ptl->L*ptl->L)/2].nOcc=0;
	for(int i=0; i<sp->polSize; i++){
		ptl->avgUnitCor[i] = 0;
		ptl->avgGenom[i] = 0;
		ptl->pc[i] = malloc(sizeof(double)*sp->polSize);
		for(int j=0; j<sp->polSize; j++)
			ptl->pc[i][j]=0;
	}
	
	for(int i=0; i<ptl->nEqd; i++){
		ptl->cmsDif[i]=0;
		ptl->smDif[i]=0;
		ptl->emDif[i]=0;
		ptl->mmDif[i]=0;
		ptl->avgShearMod[i]=0;
		ptl->avgRee[i]=0;
	}

	for(long i=0; i<sp->nTime; i++){
		ptl->rGyrT[i]=0;
		ptl->avgSL[i]=0;
	}		

	
	ptl->avgModesStat = malloc(sizeof(double*)*(sp->polSize/2+1));
	for(int i=0; i<=sp->polSize/2; i++){
		ptl->avgModesStat[i] = malloc(sizeof(double)*(sp->polSize/2+1));
		for(int j=0; j<=sp->polSize/2; j++)
			ptl->avgModesStat[i][j] = 0;
	}
	
	ptl->avgModesDyn = malloc(sizeof(double*)*(sp->polSize/2+1));
	for(int i=0; i<=sp->polSize/2; i++){
		ptl->avgModesDyn[i] = malloc(sizeof(double)*ptl->nEqd);
		for(int j=0; j<ptl->nEqd; j++)
			ptl->avgModesDyn[i][j] = 0;
	}
	
	
	for(int i=0; i<sp->nTime; i++) 
		PCInit(sp, &(ptl->polys[i]), ptl->sinfac, ptl->cosfac, ptl->modeList, ptl->nModes);
	
	TTableInit(sp, ptl);
	
// 	printf("Need an approximate %.2lf GB of memory\n", (double)sd->nSDPoints*sp->nDev*ptl->tTable.nTDT*3*sizeof(double)*1e-9); 
	long nAlloc=0;
	ptl->avgSpacDif = malloc(sizeof(double***)*sd->nSDPoints);
	for(int i=0; i<sd->nSDPoints; i++){
		ptl->avgSpacDif[i] = malloc(sizeof(double**)*sp->nDev);
		for(int j=0; j<sp->nDev; j++){
			ptl->avgSpacDif[i][j] = malloc(sizeof(double*)*ptl->tTable.nTDT);
			for(int k=0; k<ptl->tTable.nTDT; k++){
				ptl->avgSpacDif[i][j][k]=malloc(sizeof(double)*3);
				nAlloc += sizeof(double)*3;
				for(int l=0; l<3; l++)
					ptl->avgSpacDif[i][j][k][l]=0;
			}
		}
// 		printf("%i/%i: %.1lf MB allocated vs %.1lf MB\n", i, sd->nSDPoints, 1e-6*nAlloc, sp->nDev*ptl->tTable.nTDT*3*sizeof(double)*1e-6);
	}
	
	
	ptl->sAvgSpacMode = malloc(sizeof(double***)*SPAC_MODES);
	ptl->cAvgSpacMode = malloc(sizeof(double***)*SPAC_MODES);
	for(int p=0; p<SPAC_MODES; p++){
		ptl->sAvgSpacMode[p] = malloc(sizeof(double**)*ptl->tTable.nTDT);
		ptl->cAvgSpacMode[p] = malloc(sizeof(double**)*ptl->tTable.nTDT);
		for(int i=0; i<ptl->tTable.nTDT; i++){
			ptl->sAvgSpacMode[p][i] = malloc(sizeof(double*)*3);
			ptl->cAvgSpacMode[p][i] = malloc(sizeof(double*)*3);
			for(int k=0; k<3; k++){
				ptl->sAvgSpacMode[p][i][k] = malloc(sizeof(double)*3);
				ptl->cAvgSpacMode[p][i][k] = malloc(sizeof(double)*3);
				for(int j=0; j<3; j++){
					ptl->sAvgSpacMode[p][i][k][j] = 0;
					ptl->cAvgSpacMode[p][i][k][j] = 0;
				}
			}
		}
	}
	
	if(!sp->updAvgPos)
		ReadAvgPos(sp,ptl);
	else{
		ptl->avgPosition = malloc(sizeof(double*)*sp->nTime);
		for(int t=0; t<sp->nTime; t++)
			ptl->avgPosition[t] = malloc(sizeof(double)*3);
	}
}

int GetNUpdates(SimProperties* sp, char* sampleDir){
	int nUpd=0;
	char exec[10000];
	
	sp->updRouseDyn=0;
	sp->updRouseStat=0;
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/rgyr.dat -o %s/simulation_settings.txt -nt %s/rgyr_time.dat", sampleDir, sampleDir, sampleDir, sampleDir);
	sp->updGyr = !system(exec); if(sp->updGyr){ nUpd++; printf("Updating rgyr\n");}
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/rouse_stat.dat", sampleDir, sampleDir);
	sp->updRouseStat = !system(exec); if(sp->updRouseStat){ nUpd++; printf("Updating static rouse modes\n");}
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/rouse_dyn.dat", sampleDir, sampleDir);
	sp->updRouseDyn = !system(exec); if(sp->updRouseDyn){ nUpd++; printf("Updating dynamic rouse modes\n");}
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/genom.dat -o %s/simulation_settings.txt -nt %s/pgenom.dat", sampleDir, sampleDir, sampleDir, sampleDir);
	sp->updGenom = !system(exec); if(sp->updGenom){ nUpd++; printf("Updating gen dist\n");}
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/ucor.dat", sampleDir, sampleDir);
	sp->updUnit = !system(exec); if(sp->updUnit){ nUpd++; printf("Updating unit correlation\n");}
	
// 	sprintf(exec, "test %s/simulation_settings.txt -nt %s/spac_rouse.dat", sampleDir, sampleDir);
// 	sp->updSPRouse = !system(exec); if(sp->updSPRouse){ nUpd++; printf("Updating spatial rouse modes\n");}
	sp->updSPRouse = 0;
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/spac_dif.dat", sampleDir, sampleDir);
	sp->updSpacDif = !system(exec); if(sp->updSpacDif){ nUpd++; printf("Updating spatial diffusion\n");}
// 	sp->updSpacDif=0;
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/slrat.dat", sampleDir, sampleDir);
	sp->updSL = !system(exec); if(sp->updSL){ nUpd++; printf("Updating SL ratio\n");}
	
	sprintf(exec, "test %s/cmsdif_raw.dat -nt %s/cmsdif.dat", sampleDir, sampleDir);
	sp->updDif = !system(exec); if(sp->updDif){ nUpd++; printf("Updating CMS/MM/EM/SM diffusion\n");}
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/cmsdif_raw.dat", sampleDir, sampleDir);
	sp->updAvgPos = !system(exec); if(sp->updAvgPos){sp->updDif=1; nUpd++; printf("Updating raw CMS/MM/EM/SM diffusion\n");}
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/shearmod.dat", sampleDir, sampleDir);
	sp->updShearMod = !system(exec); if(sp->updShearMod){ nUpd++; printf("Updating shear modulus\n");}	
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/ree.dat", sampleDir, sampleDir);
	sp->updRee = !system(exec); if(sp->updRee){ nUpd++; printf("Updating R_ee\n");}	
	
	sprintf(exec, "test %s/simulation_settings.txt -nt %s/pc.dat", sampleDir, sampleDir);
	sp->updPC = !system(exec); if(sp->updPC){ nUpd++; printf("Updating p_c\n");}	

	return nUpd;
}

void TTableInit(SimProperties* sp, PolyTimeLapse* ptl){
	int i=0;
	int allocStep=10000;
	int nMaxStep=200;
	int nAlloc = allocStep;
	TDTTable* tTable = &(ptl->tTable);
	tTable->dt = malloc(nAlloc*sizeof(int));
	tTable->t = malloc(nAlloc*sizeof(int));
	
	for(int dt=1, dtStep=1; dt<ptl->nEqd; dt += dtStep){
		int tStep = MAX(MAX(1,dt/10), ptl->nEqd/(double)nMaxStep);
		for(int t=ptl->nTherm; t+dt<sp->nTime; t+=tStep ){
			tTable->dt[i] = dt;
			tTable->t[i] = t;
			i++;
			if(i>=nAlloc){
				nAlloc += allocStep;
				tTable->dt = realloc(tTable->dt, nAlloc*sizeof(int));
				tTable->t = realloc(tTable->t, nAlloc*sizeof(int));
			}
		}
		dtStep = MAX(1, dt/10);
	}
	tTable->nTDT = i;
	printf("Using %i measuring points\n", tTable->nTDT);
}

void TTableDestr(SimProperties* sp, PolyTimeLapse* ptl){
	TDTTable* tTable = &(ptl->tTable);
	free(tTable->t); free(tTable->dt);
}

void PCInit(SimProperties* sp,  PolyConfig* pc, double** sinfac, double** cosfac, int* modeList, int nModes){
	pc->x = malloc(sizeof(int)*sp->polSize);
	pc->y = malloc(sizeof(int)*sp->polSize);
	pc->z = malloc(sizeof(int)*sp->polSize);
	
	pc->t = malloc(sizeof(int)*sp->polSize);
	pc->u = malloc(sizeof(int)*sp->polSize);
	pc->v = malloc(sizeof(int)*sp->polSize);
	
	pc->sxmode = malloc(sizeof(double)*(sp->polSize/2+1));
	pc->symode = malloc(sizeof(double)*(sp->polSize/2+1));
	pc->szmode = malloc(sizeof(double)*(sp->polSize/2+1));
	pc->cxmode = malloc(sizeof(double)*(sp->polSize/2+1));
	pc->cymode = malloc(sizeof(double)*(sp->polSize/2+1));
	pc->czmode = malloc(sizeof(double)*(sp->polSize/2+1));
	
	pc->unitCor = malloc(sizeof(double)*sp->polSize);
	pc->genom = malloc(sizeof(double)*sp->polSize);
	pc->rGyr=0;
	
	pc->sinfac = sinfac;
	pc->cosfac = cosfac;
	pc->modeList = modeList;
	pc->nModes = nModes;
}

void DestrArrays(SimProperties* sp, PolyTimeLapse* ptl){
	PTLDestr(sp, ptl);
	free(allFiles); free(filePos); free(strIn);
	free(t); free(u); free(v);
}

void PTLDestr(SimProperties* sp, PolyTimeLapse* ptl){
	free(ptl->avgUnitCor); free(ptl->avgGenom);
	free(ptl->avgShearMod);
	for(int i=0; i<sp->nTime; i++)
		PCDestr(sp, &(ptl->polys[i]));
	free(ptl->polys);
	
	for(int i=0; i<sp->polSize/2; i++){
		free(ptl->sinfac[i]); free(ptl->cosfac[i]);
		free(ptl->avgModesDyn[i]);
		free(ptl->avgModesStat[i]);
	}
	free(ptl->avgModesDyn);
	free(ptl->avgModesStat);
	free(ptl->sinfac); free(ptl->cosfac);
	free(ptl->cmsDif);
	free(ptl->avgRee);
	free(ptl->modeList);
	
	for(int i=0; i<sp->polSize; i++){
		free(ptl->pc[i]);
	}
	free(ptl->pc);
	
	for(int p=0; p<SPAC_MODES; p++){
		for(int i=0; i<ptl->tTable.nTDT; i++){
			for(int k=0; k<3; k++){
				free(ptl->sAvgSpacMode[p][i][k]); free(ptl->cAvgSpacMode[p][i][k]);
			}
			free(ptl->sAvgSpacMode[p][i]); free(ptl->cAvgSpacMode[p][i]);
		}
		free(ptl->sAvgSpacMode[p]); free(ptl->cAvgSpacMode[p]);
	}
	for(int i=0; i<sd.nSDPoints; i++){
		for(int j=0; j<sp->nDev; j++){
			for(int k=0; k<ptl->tTable.nTDT; k++){
				free(ptl->avgSpacDif[i][j][k]);
			}
			free(ptl->avgSpacDif[i][j]);
		}
		free(ptl->avgSpacDif[i]);
	}
	free(ptl->avgSpacDif);
	
	free(ptl->sAvgSpacMode); free(ptl->cAvgSpacMode);
	TTableDestr(sp, ptl);
// 	SpacDifDestr(&(sd));
}

void PCDestr(SimProperties* sp, PolyConfig* pc){
	free(pc->x); free(pc->y); free(pc->z);
	free(pc->sxmode); free(pc->symode); free(pc->szmode);
	free(pc->cxmode); free(pc->cymode); free(pc->czmode);
	free(pc->unitCor); free(pc->genom);
}

long GetDT(char* sampleDir, char** firstFile, long* firstT){
	char exec[10000];
	char* file = malloc(10000*sizeof(char));
	char* retFile = malloc(10000*sizeof(char));
	FILE* pPipe;
	int i;
	long t;
	long minT=(long)1e13, secMinT=(long)1e14;
	
	*firstFile = retFile;
	sprintf(exec, "ls %s | grep 't=' |grep 'dev=0'", sampleDir);
	pPipe = popen(exec, "r");
	printf("Command: %s\n", exec);
	if(!pPipe){
		printf("error opening pipe\n");
		exit(0);
	}
	while(fscanf(pPipe, "%s", file)>0){
		i=0; 
		while(file[i] != '_' && file[i] != '\0') i++;
		if(file[i] != '_') continue;
		file[i]='\0';
		t = atoi(file+2);
		file[i]='_';
		if(t<minT){
			secMinT = minT;
			minT = t;
			strcpy(retFile, file);
			*firstT = t;
		}
		else if(t<secMinT){
			secMinT = t;
		}
	}
	pclose(pPipe);
	free(file);
	return secMinT-minT;
}

void SetSimProps(SimProperties* sp, char* sampleDir){
	FILE *pFile;
	char filename[10000];
	char* firstFile;
	char polType[100];
	char exec[1000];
	int nFiles;
	
	sp->sampleDir = sampleDir;
	printf("Sample directory: %s\n", sampleDir);
	sp->dT = GetDT(sampleDir, &firstFile, &(sp->tStart));
	sprintf(filename, "%s/%s", sampleDir, firstFile);
	printf("file = %s\n", filename);
	printf("dT = %li\n", sp->dT);
	pFile=fopen(filename,"r");
	fscanf(pFile, "%*s %i", &LT);
	fscanf(pFile, "%*s %i", &LU);
	fscanf(pFile, "%*s %i", &LV);
	sp->LT = LT; sp->LU=LU; sp->LV=LV;
	sp->LSIZE=sp->LT*sp->LU*sp->LV;
	if(sp->LT > MAX_LT || sp->LU > MAX_LU || sp->LV > MAX_LV){
		printf("Box is too large: (%i, %i, %i) vs (%i,%i,%i)\n", sp->LT, sp->LU, sp->LV, MAX_LT, MAX_LU, MAX_LV);
		exit(0);
	}
	
	fscanf(pFile, "%*s %i", &sp->nPol);
	fscanf(pFile, "%*s %*i");
	fscanf(pFile, "%*s %i", &sp->polSize);
	printf("nPol = %i\n", sp->nPol);
	printf("polSize=%i\n", sp->polSize);
	fclose(pFile);
	
	sprintf(exec, "ls %s | grep 't=' | grep 'dev=0' | wc -w", sampleDir);
	pFile = popen(exec, "r");
	fscanf(pFile, "%i", &sp->nTime);
	printf("Number of files detected: %i\n", sp->nTime);
	pclose(pFile);
	
	sprintf(exec, "ls %s | grep 't=' | wc -w", sampleDir);
	pFile = popen(exec, "r");
	fscanf(pFile, "%i", &nFiles);
	if(nFiles%sp->nTime != 0){
		printf("Error: uneven number of files for the different devices\n");
		exit(123);
	}
	sp->nDev = nFiles/sp->nTime;
	printf("Number of devices detected: %i\n", sp->nDev);
	pclose(pFile);
	
	sprintf(filename, "%s/simulation_settings.txt", sampleDir);
	pFile = fopen(filename, "r");
	char buf[400];
	for(int i=0; i<4; i++) fgets(buf, 400, pFile);
	fscanf(pFile, "%*s %*s %s", polType);
	for(int i=0; i<7; i++){
		fgets(buf, 400, pFile);
	}
	sp->equilibrated=0;
	fscanf(pFile, "%*s %*s %i", &sp->equilibrated);
	if(!strcmp(polType, "ring")){
		sp->polType = POL_TYPE_RING;
	}
	else if (!strcmp(polType, "lin")){
		sp->polType = POL_TYPE_LIN;
	}
	else{
		printf("Error: unknown polymer type %s\n", polType);
		exit(192);
	}
}

void InitFilePos(SimProperties* sp, int devId){
	char file[10000];
	FILE* pFile;
	
	for(int i=0; i<sp->nTime; i++){
		sprintf(file, "%s/t=%li_dev=%i.res", sp->sampleDir, i*sp->dT, devId);
		allFiles[i] = malloc(sizeof(char)*(strlen(file)+1));
		strcpy(allFiles[i], file);
		pFile = fopen(file, "r");
		if(!pFile) printf("woops: not opening file %s (%i)\n", allFiles[i], sp->nTime);
		
		for(int i=0; i<5; i++) fscanf(pFile, "%*s %*s");
		
		fgetpos(pFile, &(filePos[i]));
		fclose(pFile);
	}
}

void SpacDifInit(SpacDif* sd){
	int balVol;
// 	IDouble* distArray = malloc(sizeof(IDouble)*LT*LU*LV);
	int coor1[3];
	int maxMult=8;
	int LSIZE=LT*LU*LV;
	int *distArray = malloc(sizeof(int)*LSIZE);
	int *occLat = malloc(sizeof(int)*LSIZE);
	for(int i=0; i<LSIZE; i++){
		distArray[i]=0; occLat[i]=0;
	}
	
	sd->nSDPoints = 0;
	for(int multi=1; multi<=maxMult; multi++){
		if(LT%multi) continue;
		sd->nSDPoints += multi*multi*multi;
	}
	
	for(int i=0; i<LSIZE; i++) sd->ballList[i]=NULL;
	
	sd->sDPoints = malloc(sizeof(SpacDifPoint)*sd->nSDPoints);
	
	int sdId=0;
	long memUsed=0;
	for(int multi=1; multi<=maxMult; multi++){
		if(LT%multi) continue;
		balVol = pow(LT/multi, 3);
		
		int tuvStart = (LT/multi)/2;
		
		for(int t=tuvStart; t<LT; t+=LT/multi){
			for(int u=tuvStart; u<LT; u+=LU/multi){
				for(int v=tuvStart; v<LT; v += LV/multi){
					coor1[0]=t; coor1[1]=u; coor1[2]=v;
					int nInside = GetNClosestNeigh(occLat, distArray,coor1, balVol);
					sd->sDPoints[sdId].nCoorInside = nInside;
					for(int i=0; i<nInside; i++){
						LList* tLink = malloc(sizeof(LList));
						memUsed += sizeof(LList);
						tLink->next = sd->ballList[distArray[i]];
						tLink->sdId = sdId;
						sd->ballList[distArray[i]] = tLink;
					}
					sdId++;
				}
			}
		}
// 		printf("multi=%i, memUsed=%li, sdIDUsed=%i\n", multi, memUsed, sdId);
	}
// 	printf("%i spatial diffusion measure points\n", sd->nSDPoints);
// 	exit(0);
	printf("Warning: this function [SpacDifInit] only works when all lattices have the same size.\n");
	free(occLat); free(distArray);
}

void SpacDifDestr(SpacDif* sd){
	free(sd->sDPoints);
	for(int i=0; i<LT*LU*LV; i++){
		LList* tLink, *tNext;
		tLink = sd->ballList[i];
		while(tLink){
			tNext = tLink->next;
			free(tLink);
			tLink = tNext;
		}
	}
}

int GetNClosestNeigh(int* occLat, int* retList, int tuv[3], int volume){
	
	int nRet=1;
	int curConsider=0;
	int lastDChange=1;
	retList[0] = TUV2Coor(tuv);
	occLat[TUV2Coor(tuv)]=1;
	int LSIZE=LT*LU*LV;
	
	while(nRet<LSIZE && curConsider<nRet){
		if(nRet >= volume && curConsider == lastDChange) break;
		if(curConsider==lastDChange) lastDChange=nRet;
		if(lastDChange >= volume) break;
		int coor = retList[curConsider];
		for(int iRel=0; iRel<12; iRel++){
			int newCoor = AddCoor2TUV(coor, tuvRelPos[iRel]);
			if(newCoor > LSIZE) {printf("wtf?\n"); exit(0);}
			if(!occLat[newCoor]){
				occLat[newCoor]=1; 
				retList[nRet++]=newCoor;
			}
		}
		curConsider++;
	}
	
	for(int i=0; i<nRet; i++) occLat[retList[i]]=0;
	return nRet;
}