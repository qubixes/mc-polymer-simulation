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
	if(sp->updMagDip)    WriteMagDip(sp,ptl);
	if(sp->updDif) {
		WriteDiff(sp, ptl, ptl->cmsDif, "cms", ptl->nEqd);
		WriteDiff(sp, ptl, ptl->smDif, "sm", ptl->nEqd);
		WriteDiff(sp, ptl, ptl->emDif, "em", ptl->nEqd);
		WriteDiff(sp, ptl, ptl->mmDif, "mm", ptl->nEqd);
	}
	if(sp->updPC){
		WriteContactProbability(sp, ptl);
		WriteAvgContactProbability(sp, ptl);
	}
	if(sp->updShearMod) WriteShearMod(sp,ptl);
	if(sp->updRee) WriteRee(sp,ptl);
	if(sp->updAvgPos) WriteAvgPos(sp,ptl);
}



void WriteContactProbability(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	sprintf(file, "%s/pc.dat", sp->resDir);
	int nMono = ptl->polys[0].nMono;
	
	FILE* pFile = fopen(file, "w");
	fprintf(pFile, "#Bin= %li\n", ptl->pcBins);
	for(int iMono=0; iMono<(nMono-1)/ptl->pcBins+1; iMono++){
		for(int jMono=0; jMono<(nMono-1)/ptl->pcBins+1; jMono++){
			if(iMono>jMono)
				fprintf(pFile, "%le ", ptl->pc[iMono][jMono]/(double)(ptl->nPolAdded*ptl->nEqd*ptl->pcBins*ptl->pcBins));
			else if(iMono == jMono && ptl->pcBins > 1){
				fprintf(pFile, "%le ", ptl->pc[iMono][jMono]/(double)(ptl->nPolAdded*ptl->nEqd*((ptl->pcBins*ptl->pcBins-ptl->pcBins)/2)));
			}
			else{
				fprintf(pFile, "%le ", ptl->pc[jMono][iMono]/(double)(ptl->nPolAdded*ptl->nEqd*ptl->pcBins*ptl->pcBins));
			}
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

void WriteAvgContactProbability(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	sprintf(file, "%s/pc_avg.dat", sp->resDir);
	FILE* pFile = fopen(file, "w");
	int nMono = ptl->polys[0].nMono;
	for(long g=0; g<nMono; g++){
		fprintf(pFile, "%li %le\n", g, ptl->pcAvg[g]/(double)(ptl->nPolAdded*ptl->nEqd*(nMono-g)));
	}
	fclose(pFile);
	
	sprintf(file, "%s/pc_avg_ss.dat", sp->resDir);
	pFile = fopen(file, "w");
	for(int g=0; g<nMono; g++){
		fprintf(pFile, "%i %le\n", g, ptl->pcssAvg[g]/(double)(ptl->nPolAdded*ptl->nEqd*(nMono-g)));
	}
	fclose(pFile);
}

void WriteMagDip(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	sprintf(file, "%s/magdip.dat", sp->resDir);
	
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "#Avg_norm= %lf\n", ptl->avgMagDip/ptl->nPolAdded);
	for(int dt=0; dt<ptl->nEqd; dt++){
		if(ptl->magDipCor[dt])
			fprintf(pFile, "%li %le\n", dt*sp->dT, ptl->magDipCor[dt]/ptl->nPolAdded);
	}
	fclose(pFile);
	
	sprintf(file, "%s/magdip_time.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file for writing: %s\n", file);
		exit(192);
	}
	for(int t=0; t<sp->nTime; t++){
		fprintf(pFile, "%li %le\n", t*(long)sp->dT, ptl->magDipTime[t]/ptl->nPolAdded);
	}
	fclose(pFile);
	
	Histogram* allHist = HistogramSum(ptl->magDipHist+ptl->nTherm, sp->nTime-ptl->nTherm);
	
	int maxBin=0;
	for(int iBin=allHist->nBin-1; iBin>0; iBin--){
		if(allHist->count[iBin]){
			maxBin = iBin+1;
			break;
		}
	}
	
	sprintf(file, "%s/magdip_hist.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int iBin=0; iBin<maxBin; iBin++){
		if(allHist->count[iBin])
			fprintf(pFile, "%lf %lf %lf\n", (iBin+0.5)*allHist->dBin, allHist->avgVal[iBin]/allHist->count[iBin], allHist->count[iBin]/(double)allHist->totCount);
		else
			fprintf(pFile, "%lf %lf %lf\n", (iBin+0.5)*allHist->dBin, 0.0,0.0);
	}
	fclose(pFile);
	
	maxBin=0;
	for(int iTime=0; iTime<sp->nTime; iTime++){
		for(int iBin=ptl->magDipHist[iTime].nBin-1; iBin>=maxBin; iBin--){
			if(ptl->magDipHist[iTime].count[iBin]){
				maxBin = iBin+1;
				break;
			}
		}
	}
	
	sprintf(file, "%s/magdip_hist_time.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int iBin=0; iBin<maxBin; iBin++){
		fprintf(pFile, "%lf ", (iBin+0.5));
		for(int iTime=0; iTime<sp->nTime; iTime++){
			fprintf(pFile, "%lf ", ptl->magDipHist[iTime].count[iBin]/(double)ptl->magDipHist[iTime].totCount);
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}



void WriteRee(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
// 	printf("Writing shearmod\n");
	sprintf(file, "%s/ree.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int dt=0; dt<ptl->nEqd; dt++){
		if(ptl->avgRee[dt])
			fprintf(pFile, "%li %le\n", dt*sp->dT, ptl->avgRee[dt]/ptl->nPolAdded);
	}
	fclose(pFile);
}


void WriteShearMod(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
// 	printf("Writing shearmod\n");
	sprintf(file, "%s/shearmod.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int dt=1; dt<ptl->nEqd; dt++){
		if(ptl->avgShearMod[dt])
			fprintf(pFile, "%li %le\n", dt*sp->dT, ptl->avgShearMod[dt]/ptl->nPolAdded);
	}
	
	fclose(pFile);
}

void WriteGyration(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
	
	sprintf(file, "%s/rgyr.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "%le\n", ptl->avgRGyr/(ptl->nEqd*ptl->nPolAdded));
	fclose(pFile);
	
	sprintf(file, "%s/rgyr_time.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int t=0; t<sp->nTime; t++){
		fprintf(pFile, "%li %le\n", sp->dT*t, ptl->rGyrT[t]/ptl->nPolAdded);
	}
	fclose(pFile);
}

void WriteGenom(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	int g;
	int nMono = ptl->polys[0].nMono;
	FILE* pFile;
	
	sprintf(file, "%s/genom.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	for(int ig=0; ig<ptl->nIg; ig++){
		g = ptl->genomIdList[ig];
		if(ptl->avgGenom[ig])
			fprintf(pFile, "%i %le\n", g, ptl->avgGenom[ig]/(double)ptl->genomCount[ig]);
	}
	
	fclose(pFile);
	
	sprintf(file, "%s/pgenom.dat", sp->resDir);
	
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file: %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "# ");
	for(int g=0; g<nMono; g++){
		if(ptl->avgGenom[g])
			fprintf(pFile, "%i ", g);
	}
	fprintf(pFile, "\n");
	for(int i=0; i<ptl->nGenomBin; i++){
		for(int ig=0; ig<ptl->nIg; ig++){
			int g = ptl->genomIdList[ig];
			if(ptl->avgGenom[g]){
				if(ptl->genomProb[ig][i]==0)
					fprintf(pFile, "nan nan ");
				else{
					fprintf(pFile, "%le ", ptl->genomProb[ig][i]/(double)(ptl->genomCount[ig]*ptl->pointDensity[i]));
					fprintf(pFile, "%le ", ptl->genomR[ig][i]/(double)ptl->genomProb[ig][i]);
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
	sprintf(file, "%s/spac_rouse.dat", sp->resDir);
	pFile = fopen(file, "w");
	
	for(int p=1; p<SPAC_MODES; p++) sqModes[p]=0;
	
	for(int i=0; i<=ptl->tTable.nTDT; i++){
		if(i<ptl->tTable.nTDT) dt = ptl->tTable.tdt[i].dt;
		else dt = -1;
// 		printf("[%i, %i]\n", dt, ptl->tTable.t[i]);
		if(dt != prevDT){
			
			fprintf(pFile, "%li ", prevDT*sp->dT);
			for(int p=1; p<SPAC_MODES; p++){
				fprintf(pFile, "%lf ", sqModes[p]/(double)(count*ptl->nPolAdded));
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
	int polSize = ptl->polys[0].polSize;
	
	sprintf(file, "%s/ucor.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	int gmax;
	if(ptl->polys[0].polType == POL_TYPE_LIN)
		gmax = polSize;
	else
		gmax = polSize/2;
	
	for(g=1; g<=gmax; g++){
		if(unitCor[g])
			fprintf(pFile, "%i %le\n", g, unitCor[g]/(ptl->nEqd*ptl->nPolAdded));
	}
	
	fclose(pFile);
}

void WriteModesStat(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	
	sprintf(file, "%s/rouse_stat.dat", sp->resDir);
	pFile = fopen(file, "w");
	if(!pFile)
		return;
	
	fprintf(pFile, "# ");
	for(int ip=0; ip<ptl->nModes; ip++) fprintf(pFile, "%i ", ptl->modeList[ip]);
	fprintf(pFile, "\n");
	for(int ip=0; ip<ptl->nModes; ip++){
// 		int p = ptl->modeList[ip];
		for(int iq=0; iq<ptl->nModes; iq++){
// 			int q = ptl->modeList[iq];
			fprintf(pFile, "%le ", ptl->avgModesStat[ip][iq]/(double)(ptl->nPolAdded));
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
	
	sprintf(file, "%s/rouse_dyn.dat", sp->resDir);
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
// 			int p = ptl->modeList[ip];
			fprintf(pFile, "%le ", ptl->avgModesDyn[ip][dt]/(double)(ptl->nPolAdded));
		}
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}



void WriteDiff(SimProperties* sp, PolyTimeLapse* ptl, double* dif, char* base, long nEqd){
	char file[1000];
	FILE* pFile;
	
	if(sp->updAvgPos)
		sprintf(file, "%s/%sdif_raw.dat", sp->resDir, base);
	else
		sprintf(file, "%s/%sdif.dat", sp->resDir, base);
	pFile = fopen(file, "w"); if(!pFile) return;

	for(int dt=1; dt<nEqd; dt++){
		if(!dif[dt]) continue;
		fprintf(pFile, "%li %le\n", dt*sp->dT,  dif[dt]/(double)(ptl->nPolAdded));
	}
	fclose(pFile);
}


void WriteCMSDiff(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
	
	double* cmsDif = ptl->cmsDif;
	
	if(sp->updAvgPos)
		sprintf(file, "%s/cmsdif_raw.dat", sp->resDir);
	else
		sprintf(file, "%s/cmsdif.dat", sp->resDir);
	pFile = fopen(file, "w"); if(!pFile) return;
	
	for(int dt=0; dt<ptl->nEqd; dt++){
		if(!cmsDif[dt]) continue;
		fprintf(pFile, "%li %le\n", dt*sp->dT, cmsDif[dt]/(double)(ptl->nPolAdded));
	}
	fclose(pFile);
}
	
void WriteSL(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	FILE* pFile;
	double* slRat = ptl->avgSL;
	
	sprintf(file, "%s/slrat.dat", sp->resDir);
	pFile = fopen(file, "w"); if(!pFile) return;
	
	for(int t=0; t<sp->nTime; t++){
		fprintf(pFile, "%li %le\n", t*sp->dT,  slRat[t]/(double)(ptl->nPolAdded));
	}
	fclose(pFile);
}

void LoadPTL(SimProperties* sp, PolyTimeLapse* ptl, int polId, int devId){
	char file[1000];
	int tLast, uLast, vLast, step;
	long nMono=-1;
	
	ptl->polId = polId; ptl->devId = devId;
	sprintf(file, "%s/ptl/pol=%i_dev=%i.res", sp->sampleDir, polId, devId);
	FILE* pFile = fopen(file, "r");
	if(!pFile){ printf("Error loading file %s\n", file); exit(0); }
	
	for(int i=0; i<2; i++) fscanf(pFile, "%*s %*s");
	fscanf(pFile, "%*s %li", &nMono);
	
	char* strIn = malloc(sizeof(char)*(nMono+1));
	int* t = malloc(sizeof(int)*nMono);
	int* u = malloc(sizeof(int)*nMono);
	int* v = malloc(sizeof(int)*nMono);
	
	double tAvgLast, uAvgLast, vAvgLast;
	for(int iTime=0; iTime<sp->nTime; iTime++){
		fscanf(pFile, "%d %d %d", t, u, v);
		fscanf(pFile, "%s", strIn);

		if(iTime != 0){
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
		
		if(strlen(strIn) != nMono){
			printf("\n\nError reading PTL file (rerun?): %li vs %li\nfile=%s\n\n", strlen(strIn), nMono, file);
			exit(0);
		}
		
		for(int iMono=0; iMono<nMono-1; iMono++){
			step = CharToHex(strIn[iMono]);
			t[iMono+1] = t[iMono]+((step>>0)&0x1)-((step>>3)&0x1);
			u[iMono+1] = u[iMono]+((step>>1)&0x1)-((step>>3)&0x1);
			v[iMono+1] = v[iMono]+((step>>2)&0x1)-((step>>3)&0x1);
		}
		
		double tAvg=0, uAvg=0, vAvg=0;
		for(int iMono=0; iMono<nMono; iMono++){
			tAvg += t[iMono]/(double)nMono;
			uAvg += u[iMono]/(double)nMono;
			vAvg += v[iMono]/(double)nMono;
		}
		
		int dt=0, du=0, dv=0;
		if(iTime != 0){
			while(tAvg-tAvgLast+dt >  LT/2) dt -= LT;
			while(tAvg-tAvgLast+dt < -LT/2) dt += LT;
			while(uAvg-uAvgLast+du >  LU/2) du -= LU;
			while(uAvg-uAvgLast+du < -LU/2) du += LU;
			while(vAvg-vAvgLast+dv >  LU/2) dv -= LV;
			while(vAvg-vAvgLast+dv < -LU/2) dv += LV;
		}
		
		tAvgLast = tAvg+dt;
		uAvgLast = uAvg+du;
		vAvgLast = vAvg+dv;
		
		
		for(int iMono=0; iMono<nMono; iMono++){
			t[iMono] += dt;
			u[iMono] += du;
			v[iMono] += dv;
		}
		
// 		printf("nMono = %li\n", nMono);
// 		printf("maxMono = %li\n", sp->maxNMono);
		for(int iMono=0; iMono<nMono; iMono++){
			ptl->polys[iTime].x[iMono] = t[iMono]+u[iMono]-v[iMono];
			ptl->polys[iTime].y[iMono] = t[iMono]-u[iMono];
			ptl->polys[iTime].z[iMono] = v[iMono];
			
			ptl->polys[iTime].t[iMono] = t[iMono];
			ptl->polys[iTime].u[iMono] = u[iMono];
			ptl->polys[iTime].v[iMono] = v[iMono];
		}
		ptl->polys[iTime].nMono = nMono;
		if(strIn[nMono-1] == 'f'){
			ptl->polys[iTime].polType = POL_TYPE_LIN;
			ptl->polys[iTime].polSize = nMono-1;
		}
		else{
			ptl->polys[iTime].polType = POL_TYPE_RING;
			ptl->polys[iTime].polSize = nMono;
		}
	}
	fclose(pFile);
}

double** ReadAvgPos(SimProperties* sp, PolyTimeLapse* ptl){
	char file[1000];
	
	sprintf(file, "%s/avg_pos.dat", sp->sampleDir);
	FILE* pFile = fopen(file, "r");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	double** avgPosition = malloc(sizeof(double*)*sp->nTime);
	double x,y,z;
	int t=0;
	while( fscanf(pFile, "%le %le %le",  &x, &y, &z) == 3 ){
		if(t >= sp->nTime){
			printf("Error reading file %s: too many lines (>=%i)\n", file, t);
			exit(192);
		}
		
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
			fprintf(pFile, "%le ", ptl->avgPosition[t][k]/(double)(ptl->nPolAdded));
		fprintf(pFile, "\n");
	}
	fclose(pFile);
}

void WriteSpacDif(SimProperties* sp, PolyTimeLapse* ptl){
	FILE* pFile;
	char file[1000];
	sprintf(file, "%s/spac_dif.dat", sp->resDir);
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
