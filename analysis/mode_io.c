#include "lowm_modes.h"

void WriteAllFiles(SimProperties* sp, PolyTimeLapse *ptl){
	if(sp->updRouseStat) WriteModesStat(sp, ptl);
	if(sp->updRouseDyn)  WriteModesDyn(sp, ptl);
	if(sp->updGyr)       WriteGyration(sp, ptl);
	if(sp->updGenom)     WriteGenom(sp, ptl);
	if(sp->updUnit)      WriteUnitCor(sp, ptl);
	if(sp->updSL)        WriteSL(sp, ptl);
	if(sp->updSPRouse)   WriteSpaceRouse(sp, ptl);
	if(sp->updSpacDif)   WriteSpacDif(sp,ptl);
	if(sp->updDif) {
// 		printf("nEqd=%li\n", ptl->nEqd);
		WriteDiff(sp, ptl->cmsDif, "cms", ptl->nEqd);
		WriteDiff(sp, ptl->smDif, "sm", ptl->nEqd);
		WriteDiff(sp, ptl->emDif, "em", ptl->nEqd);
		WriteDiff(sp, ptl->mmDif, "mm", ptl->nEqd);
	}
	if(sp->updPC) WriteContactProbability(sp, ptl);
	if(sp->updShearMod) WriteShearMod(sp,ptl);
	if(sp->updRee) WriteRee(sp,ptl);
	if(sp->updAvgPos) WriteAvgPos(sp,ptl);
}

void WriteContactProbability(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	sprintf(file, "%s/pc.dat", sp->sampleDir);
	FILE* pFile = fopen(file, "w");
	for(int iMono=0; iMono<sp->polSize; iMono++){
		for(int jMono=0; jMono<sp->polSize; jMono++){
			if(iMono>jMono)
				fprintf(pFile, "%le ", ptl->pc[iMono][jMono]/(double)(sp->nPol*sp->nDev*ptl->nEqd));
			else
				fprintf(pFile, "%le ", ptl->pc[jMono][iMono]/(double)(sp->nPol*sp->nDev*ptl->nEqd));
		}
		fprintf(pFile, "\n");
	}
}

void WriteRee(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
// 	printf("Writing shearmod\n");
	sprintf(file, "%s/ree.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int dt=0; dt<ptl->nEqd; dt++){
		if(ptl->avgRee[dt])
			fprintf(pFile, "%li %le\n", dt*sp->dT, ptl->avgRee[dt]/(sp->nPol)/sp->nDev);
	}
	fclose(pFile);
}


void WriteShearMod(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
// 	printf("Writing shearmod\n");
	sprintf(file, "%s/shearmod.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int dt=1; dt<ptl->nEqd; dt++){
		if(ptl->avgShearMod[dt])
			fprintf(pFile, "%li %le\n", dt*sp->dT, ptl->avgShearMod[dt]/(sp->nPol)/sp->nDev);
	}
	
	fclose(pFile);
}

void WriteGyration(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
	
	sprintf(file, "%s/rgyr.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "%le\n", ptl->avgRGyr/(sp->nPol*ptl->nEqd*sp->nDev));
	fclose(pFile);
	
	sprintf(file, "%s/rgyr_time.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int t=0; t<sp->nTime; t++){
		fprintf(pFile, "%li %le\n", sp->dT*t, ptl->rGyrT[t]/(sp->nPol*sp->nDev));
	}
	fclose(pFile);
}

void WriteGenom(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	int g;
	FILE* pFile;
	
	sprintf(file, "%s/genom.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(g=1; g<sp->polSize; g++){
		if(ptl->avgGenom[g])
			fprintf(pFile, "%i %le\n", g, ptl->avgGenom[g]/(double)ptl->genomCount[g]);
	}
	
	fclose(pFile);
	
	sprintf(file, "%s/pgenom.dat", sp->sampleDir);
	
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file: %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "# ");
	for(int g=0; g<sp->polSize; g++){
		if(ptl->avgGenom[g])
			fprintf(pFile, "%i ", g);
	}
	fprintf(pFile, "\n");
	for(int i=0; i<sp->polSize; i++){
		for(int g=0; g<sp->polSize; g++){
			if(ptl->avgGenom[g]){
				if(ptl->genomProb[g][i]==0)
					fprintf(pFile, "nan nan ");
				else{
					fprintf(pFile, "%le ", ptl->genomProb[g][i]/(double)(ptl->genomCount[g]*ptl->pointDensity[i]));
					fprintf(pFile, "%le ", ptl->genomR[g][i]/(double)ptl->genomProb[g][i]);
				}
			}
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
	
}

void WriteSpaceRouse(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	int count=0;
	int prevDT=1, dt;
	double sqModes[SPAC_MODES];
	sprintf(file, "%s/spac_rouse.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	
	for(int p=1; p<SPAC_MODES; p++) sqModes[p]=0;
	
	for(int i=0; i<=ptl->tTable.nTDT; i++){
		if(i<ptl->tTable.nTDT) dt = ptl->tTable.tdt[i].dt;
		else dt = -1;
// 		printf("[%i, %i]\n", dt, ptl->tTable.t[i]);
		if(dt != prevDT){
			
			fprintf(pFile, "%li ", prevDT*sp->dT);
			for(int p=1; p<SPAC_MODES; p++){
				fprintf(pFile, "%lf ", sqModes[p]/(double)(count*sp->nPol*sp->nDev));
				sqModes[p]=0;
			}
			fprintf(pFile, "\n");
			count=0;
			if(i==ptl->tTable.nTDT) break;
			prevDT = dt;
		}

		for(int p=1; p<SPAC_MODES; p++){
			for(int k=0; k<3; k++){
				for(int j=0; j<3; j++){
					sqModes[p] += ptl->sAvgSpacMode[p][i][k][j]*ptl->sAvgSpacMode[p][i][k][j];
					sqModes[p] += ptl->cAvgSpacMode[p][i][k][j]*ptl->cAvgSpacMode[p][i][k][j];
				}
				count +=2;
			}
		}
	}
	fclose(pFile);
}

void WriteUnitCor(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	int g;
	FILE* pFile;
	double* unitCor = ptl->avgUnitCor;
	
	sprintf(file, "%s/ucor.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	int gmax;
	if(sp->polType == POL_TYPE_LIN)
		gmax = sp->polSize-1;
	else
		gmax = sp->polSize/2;
	
	for(g=1; g<=gmax; g++){
		if(unitCor[g])
			fprintf(pFile, "%i %le\n", g, unitCor[g]/(sp->nPol*ptl->nEqd)/sp->nDev);
	}
	
	fclose(pFile);
}

void WriteModesStat(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	
	sprintf(file, "%s/rouse_stat.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile)
		return;
	
	fprintf(pFile, "# ");
	for(int ip=0; ip<ptl->nModes; ip++) fprintf(pFile, "%i ", ptl->modeList[ip]);
	fprintf(pFile, "\n");
	for(int ip=0; ip<ptl->nModes; ip++){
		int p = ptl->modeList[ip];
		for(int iq=0; iq<ptl->nModes; iq++){
			int q = ptl->modeList[iq];
			fprintf(pFile, "%le ", ptl->avgModesStat[p][q]/(double)(sp->nPol*sp->nDev));
// 			printf("(%i,%i,%lf) ", p,q,ptl->avgModesStat[p][q]);
		}
		fprintf(pFile, "\n");
// 		printf("\n");
	}
	fclose(pFile);
}


void WriteModesDyn(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	
	sprintf(file, "%s/rouse_dyn.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	if(!pFile) return;
	
	fprintf(pFile, "# ");
	for(int ip=0; ip<ptl->nModes; ip++) fprintf(pFile, "%i ", ptl->modeList[ip]);
	fprintf(pFile, "\n");
	
	for(int dt=0; dt<ptl->nEqd; dt++){
		if(!ptl->avgModesDyn[0][dt])
			continue;
		fprintf(pFile, "%li ", sp->dT*dt);
		for(int ip=0; ip<ptl->nModes; ip++){
			int p = ptl->modeList[ip];
			fprintf(pFile, "%le ", ptl->avgModesDyn[p][dt]/(double)(sp->nPol*sp->nDev));
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}



void WriteDiff(SimProperties* sp, double* dif, char* base, long nEqd){
	char file[1000];
	FILE* pFile;
	
	if(sp->updAvgPos)
		sprintf(file, "%s/%sdif_raw.dat", sp->sampleDir, base);
	else
		sprintf(file, "%s/%sdif.dat", sp->sampleDir, base);
	pFile = fopen(file, "w"); if(!pFile) return;

	for(int dt=1; dt<nEqd; dt++){
		if(!dif[dt]) continue;
		fprintf(pFile, "%li %le\n", dt*sp->dT,  dif[dt]/(double)(sp->nPol*sp->nDev));
	}
	fclose(pFile);
}


void WriteCMSDiff(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
	
	double* cmsDif = ptl->cmsDif;
	
	if(sp->updAvgPos)
		sprintf(file, "%s/cmsdif_raw.dat", sp->sampleDir);
	else
		sprintf(file, "%s/cmsdif.dat", sp->sampleDir);
	pFile = fopen(file, "w"); if(!pFile) return;
	
	for(int dt=0; dt<ptl->nEqd; dt++){
		if(!cmsDif[dt]) continue;
		fprintf(pFile, "%li %le\n", dt*sp->dT, cmsDif[dt]/(double)(sp->nPol*sp->nDev));
	}
	fclose(pFile);
}
	
void WriteSL(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
	double* slRat = ptl->avgSL;
	
	sprintf(file, "%s/slrat.dat", sp->sampleDir);
	pFile = fopen(file, "w"); if(!pFile) return;
	
	for(int t=0; t<sp->nTime; t++){
		fprintf(pFile, "%li %le\n", t*sp->dT,  slRat[t]/(double)(sp->nPol*sp->nDev));
	}
	fclose(pFile);
}

void LoadPTL(SimProperties* sp, PolyTimeLapse* ptl, int polId, int devId){
	char file[1000];
	int tLast, uLast, vLast, step;
	
	ptl->polId = polId; ptl->devId = devId;
	sprintf(file, "%s/ptl/pol=%i_dev=%i.res", sp->sampleDir, polId, devId);
	FILE* pFile = fopen(file, "r");
	if(!pFile){ printf("Error loading file %s\n", file); exit(0);}
	
	for(int i=0; i<3; i++) fscanf(pFile, "%*s %*s");
	
// 	for(int i=0; i<THERM; i++) fscanf(pFile, "%*i %*i %*i %*s");
	for(int i=0; i<sp->nTime; i++){
		fscanf(pFile, "%d %d %d", t, u, v);
		fscanf(pFile, "%s", strIn);

		if(i != 0){
			while(t[0]-tLast >  LT/2) t[0] -= LT;
			while(t[0]-tLast < -LT/2) t[0] += LT;
			while(u[0]-uLast >  LU/2) u[0] -= LU;
			while(u[0]-uLast < -LU/2) u[0] += LU;
			while(v[0]-vLast >  LU/2) v[0] -= LV;
			while(v[0]-vLast < -LU/2) v[0] += LV;
		}
		tLast = t[0];
		uLast = u[0];
		vLast = v[0];
		
		if(strlen(strIn) != sp->polSize){
			printf("\n\nMeh: %li vs %i\nfile=%s\n\n", strlen(strIn), sp->polSize, file);
			exit(0);
		}
		
		for(int iMono=0; iMono<sp->polSize-1; iMono++){
			step = CharToHex(strIn[iMono]);
			t[iMono+1] = t[iMono]+((step>>0)&0x1)-((step>>3)&0x1);
			u[iMono+1] = u[iMono]+((step>>1)&0x1)-((step>>3)&0x1);
			v[iMono+1] = v[iMono]+((step>>2)&0x1)-((step>>3)&0x1);
		}
		
		for(int iMono=0; iMono<sp->polSize; iMono++){
// 			ptl->polys[i-THERM].x[iMono] = v[iMono]-t[iMono];
// 			ptl->polys[i-THERM].y[iMono] = t[iMono]+v[iMono]-u[iMono];
// 			ptl->polys[i-THERM].z[iMono] = u[iMono];
			
			ptl->polys[i].x[iMono] = t[iMono]+u[iMono]-v[iMono];
			ptl->polys[i].y[iMono] = t[iMono]-u[iMono];
			ptl->polys[i].z[iMono] = v[iMono];
			
			ptl->polys[i].t[iMono] = t[iMono];
			ptl->polys[i].u[iMono] = u[iMono];
			ptl->polys[i].v[iMono] = v[iMono];
		}
	}
	fclose(pFile);
}

double** ReadAvgPos(SimProperties* sp){
	char file[1000];
	
	sprintf(file, "%s/avg_pos.dat", sp->sampleDir);
	FILE* pFile = fopen(file, "r");
	double** avgPosition = malloc(sizeof(double*)*sp->nTime);
	double x,y,z;
	int t=0;
	while( fscanf(pFile, "%le %le %le",  &x, &y, &z) == 3 ){
		avgPosition[t] = malloc(sizeof(double)*3);
		avgPosition[t][0] = x;
		avgPosition[t][1] = y;
		avgPosition[t][2] = z;
		t++;
	}
	if(t!= sp->nTime){
		printf("Error reading file %s: wrong number of lines: have %i, need %i\n", file, t, sp->nTime);
		exit(195);
	}
	fclose(pFile);
	return avgPosition;
}

void WriteAvgPos(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	sprintf(file, "%s/avg_pos.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	
	for(int t=0; t<sp->nTime; t++){
		for(int k=0; k<3; k++)
			fprintf(pFile, "%le ", ptl->avgPosition[t][k]/(double)(sp->nPol*sp->nDev));
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

void WriteSpacDif(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	sprintf(file, "%s/spac_dif.dat", sp->sampleDir);
	pFile = fopen(file, "w");
	double difFac, difBox;
	
	fprintf(pFile, "# ");
	for(int i=0; i<sd.nSDPoints; i++){
		if(i==0 || sd.sDPoints[i].nCoorInside != sd.sDPoints[i-1].nCoorInside){
			fprintf(pFile, "%i ", sd.sDPoints[i].nCoorInside);
		}
	}
	fprintf(pFile, "\n");
	
	
	for(int tBlockS=0; tBlockS<ptl->tTable.nTDT; ){
		int tBlockE=tBlockS+1;
		while(tBlockE<ptl->tTable.nTDT && ptl->tTable.tdt[tBlockS].dt == ptl->tTable.tdt[tBlockE].dt) tBlockE++;
		fprintf(pFile, "%li ", ptl->tTable.tdt[tBlockS].dt*sp->dT);
		for(int sBlockS=0; sBlockS<sd.nSDPoints; ){
			int sBlockE=sBlockS+1;
			while(sBlockE<sd.nSDPoints && sd.sDPoints[sBlockS].nCoorInside == sd.sDPoints[sBlockE].nCoorInside) 
				sBlockE++;
			
// 			printf("t/dt block [%i-%i], sd block [%i-%i]\n", tBlockS, tBlockE, sBlockS, sBlockE);
			double sumDif=0; long num=0;
			for(int itdt=tBlockS; itdt<tBlockE; itdt++){
				for(int sdId=sBlockS; sdId<sBlockE; sdId++){
					for(int devId=0; devId<sp->nDev; devId++){
						for(int i=0; i<3; i++){
							difFac=ptl->avgSpacDif[sdId][devId][itdt][i]/(double)(sd.sDPoints[sBlockS].nCoorInside);
							if(sdId == 0)
								difBox=0;
							else
								difBox=ptl->avgSpacDif[0][devId][itdt][i]/(double)(sd.sDPoints[0].nCoorInside);
							sumDif += pow(difFac-difBox, 2);
// 							if(ptl->tTable.dt[itdt]==1)printf("sumdif=%lf\n", sumDif);
						}
						num++;
					}
				}
			}
			fprintf(pFile, "%le ", sumDif/(double)(num));
			sBlockS=sBlockE;
		}
		tBlockS = tBlockE;
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}