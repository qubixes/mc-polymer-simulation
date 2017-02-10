#include "lowm_modes.h"

void ComputeObsPoly(SimProperties* sp, PolyConfig* pcfg){
	ComputeCMS(sp, pcfg);
	ComputeTUVCMS(sp, pcfg);
	if(sp->updRouseStat || sp->updRouseDyn) ComputeModes(sp, pcfg);

	if(sp->updUnit)      ComputeUnitCor(sp, pcfg);
	if(sp->updGyr)       ComputeGyration(sp, pcfg);
// 	if(sp->updGenom)     ComputeGenom(sp, pcfg);
	if(sp->updSL)        ComputeSL(sp, pcfg);
	if(sp->updShearMod)  ComputeShearMod(sp,pcfg);
	if(sp->updRee)       ComputeRee(sp, pcfg);
}


void AddAverages(SimProperties* sp, PolyTimeLapse* ptl){
	PolyConfig* pcfg;
	Timer time;
	
	TimerStart(&time);
	for(pcfg = ptl->polys; pcfg < ptl->polys+sp->nTime; pcfg++)
		ComputeObsPoly(sp, pcfg);
	printf("[Compute: %.2lf ms] ", 1e3*TimerElapsed(&time)); fflush(NULL);
	TimerStart(&time);
	if(sp->updRouseStat) AddModesStat(sp,ptl);
	if(sp->updRouseDyn) AddModesDyn(sp,ptl);
	if(sp->updPC)      AddContactProbability(sp,ptl);
// 	if(sp->updGenom)   AddGenom(sp,ptl);
	if(sp->updGenom)   AddGenomNew2(sp,ptl);
	if(sp->updUnit)    AddUnit(sp,ptl);
	if(sp->updSPRouse) AddSpaceRouse(sp, ptl);
	if(sp->updSpacDif) AddSpacDif(sp,ptl,&sd);
	if(sp->updSL)      AddSL(sp,ptl);
	if(sp->updGyr)     AddRGyr(sp,ptl);
	if(sp->updDif){
		AddEM(sp,ptl); AddCMS(sp,ptl); AddSM(sp,ptl); 
		AddMM(sp,ptl); AddEM(sp,ptl);
	}
	if(sp->updAvgPos) AddAvgPos(sp,ptl);
	if(sp->updShearMod) AddShearMod(sp,ptl);
	if(sp->updRee) AddRee(sp, ptl);
	printf("[Add: %.2lf ms]", 1e3*TimerElapsed(&time)); fflush(NULL);
}


void ComputeRee(SimProperties* sp, PolyConfig* pcfg){
	int monoId;
	if(sp->polType == POL_TYPE_LIN)
		monoId = sp->polSize-1;
	else
		monoId = sp->polSize/2;
	pcfg->ree[0] = (pcfg->x[0]-pcfg->x[monoId])/sqrt(2);
	pcfg->ree[1] = (pcfg->y[0]-pcfg->y[monoId])/sqrt(2);
	pcfg->ree[2] = (pcfg->z[0]-pcfg->z[monoId])/sqrt(2);
}

void ComputeShearMod(SimProperties* sp, PolyConfig* pcfg){
	
	for(int i=0; i<3; i++) pcfg->stress[i]=0;
	
	for(int i=0; i<sp->polSize-1; i++){
		double dx = pcfg->x[i+1]-pcfg->x[i];
		double dy = pcfg->y[i+1]-pcfg->y[i];
		double dz = pcfg->z[i+1]-pcfg->z[i];
		pcfg->stress[0] += dx*dy;
		pcfg->stress[1] += dy*dz;
		pcfg->stress[2] += dx*dz;
// 		printf("[dx dy dz] %lf %lf %lf\n", dx*dy, dy*dz, dx*dz);
	}
	for(int i=0; i<3; i++)
		pcfg->stress[i] /=2;
// 	printf("%lf %lf %lf\n", pcfg->stress[0], pcfg->stress[1], pcfg->stress[2]);
}

void ComputeGenom(SimProperties* sp, PolyConfig* pcfg){
	int di=20;
	int dx, dy, dz;
	int gMax;
	if(sp->polType == POL_TYPE_LIN) gMax=sp->polSize-1;
	else gMax = sp->polSize/2;
	for(int g=1, dg=1; g<=gMax; g+=dg){
		int num=0;
		pcfg->genom[g]=0;
		di = MAX(1,MAX(sp->polSize/100, dg/10));
		for(int i=0; i<sp->polSize; i+=di){
			if(i+g >= sp->polSize && sp->polType == POL_TYPE_LIN) continue;
			int j = (i+g)%sp->polSize;
			
			dx = pcfg->x[i]-pcfg->x[j];
			dy = pcfg->y[i]-pcfg->y[j];
			dz = pcfg->z[i]-pcfg->z[j];
			
			pcfg->genom[g] += (dx*dx+dy*dy+dz*dz)/2.0;
			num++;
		}
		
		pcfg->genom[g] /= num;
		if(!num){
			printf("woops: g = %i, dg=%i\n", g, dg);
			exit(0);
		}
		dg = MAX(1,MIN(gMax/2-g, dg/10));
	}
}

void ComputeModes(SimProperties* sp, PolyConfig* pcfg){
	for(int ip=0; ip<pcfg->nModes; ip++){
		int p = pcfg->modeList[ip];
		pcfg->sxmode[p]=0; pcfg->symode[p]=0; pcfg->szmode[p]=0; 
		pcfg->cxmode[p]=0; pcfg->cymode[p]=0; pcfg->czmode[p]=0;
		for(int k=0;k<sp->polSize;k++){
			pcfg->sxmode[p]+=pcfg->sinfac[p][k]*pcfg->x[k];
			pcfg->symode[p]+=pcfg->sinfac[p][k]*pcfg->y[k];
			pcfg->szmode[p]+=pcfg->sinfac[p][k]*pcfg->z[k];
			pcfg->cxmode[p]+=pcfg->cosfac[p][k]*pcfg->x[k];
			pcfg->cymode[p]+=pcfg->cosfac[p][k]*pcfg->y[k];
			pcfg->czmode[p]+=pcfg->cosfac[p][k]*pcfg->z[k];
		}
	}
}

void ComputeSL(SimProperties* sp, PolyConfig* pcfg){
	int nSl=0;
	for(int i=0; i<sp->polSize-1; i++){
		if(pcfg->x[i] == pcfg->x[i+1] && pcfg->y[i] == pcfg->y[i+1] && pcfg->z[i] == pcfg->z[i+1]) 
			nSl++;
	}
	
	if(sp->polType == POL_TYPE_RING){
		int i=sp->polSize-1, j=0;
		if(pcfg->x[i] == pcfg->x[j] && pcfg->y[i] == pcfg->y[j] && pcfg->z[i] == pcfg->z[j])
			nSl++;
		pcfg->slRat = nSl/(double)sp->polSize;
	}
	else{
		pcfg->slRat = nSl/(double)(sp->polSize-1);
	}
}

void ComputeUnitCor(SimProperties* sp, PolyConfig* pcfg){
	int di=5;
	int dxi, dxj, dyi,dyj,dzi,dzj;
	int gMax;
	
	if(sp->polType == POL_TYPE_LIN)
		gMax =sp->polSize-1;
	else
		gMax =sp->polSize/2;
	
	for(int g=1, dg=1; g<=gMax; g+=dg){
		int num=0;
		pcfg->unitCor[g] = 0;
		di = MAX(1,MAX(sp->polSize/50, dg/10));
		for(int i=0; i<sp->polSize; i+=di){
			if(i+g+1 >= sp->polSize && sp->polType == POL_TYPE_LIN) continue;
			int j = (i+g)%sp->polSize;
			
			dxi = pcfg->x[i]-pcfg->x[(i+1)%sp->polSize];
			dyi = pcfg->y[i]-pcfg->y[(i+1)%sp->polSize];
			dzi = pcfg->z[i]-pcfg->z[(i+1)%sp->polSize];
			
			dxj = pcfg->x[j]-pcfg->x[(j+1)%sp->polSize];
			dyj = pcfg->y[j]-pcfg->y[(j+1)%sp->polSize];
			dzj = pcfg->z[j]-pcfg->z[(j+1)%sp->polSize];
			
			pcfg->unitCor[g] += (dxi*dxj+dyi*dyj+dzi*dzj)/2.0;
			num++;
		}
		
		pcfg->unitCor[g] /= num;
		
		dg = MAX(1,MIN(gMax-g, dg/10));
	}
}

void ComputeGyration(SimProperties* sp, PolyConfig* pcfg){
	double dx, dy, dz;
	double rg=0;
	for(int i=0; i<sp->polSize; i++){
		dx = pcfg->x[i]-pcfg->cms.x;
		dy = pcfg->y[i]-pcfg->cms.y;
		dz = pcfg->z[i]-pcfg->cms.z;
		rg += (dx*dx+dy*dy+dz*dz)/2;
	}
	rg /= sp->polSize;
	pcfg->rGyr = rg;
}

void ComputeCMS(SimProperties* sp, PolyConfig* pcfg){
	pcfg->cms.x=0.0; 
	pcfg->cms.y=0.0; 
	pcfg->cms.z=0.0;
	for(int i=0; i<sp->polSize; i++){
		pcfg->cms.x += pcfg->x[i];
		pcfg->cms.y += pcfg->y[i];
		pcfg->cms.z += pcfg->z[i];
	}
	pcfg->cms.x /= sp->polSize;
	pcfg->cms.y /= sp->polSize;
	pcfg->cms.z /= sp->polSize;
// 	printf("%lf %lf %lf\n", pcfg->cms.x, pcfg->cms.y, pcfg->cms.z);
}

void ComputeTUVCMS(SimProperties* sp, PolyConfig* pcfg){
	for(int k=0; k<3; k++) pcfg->tuvCMS[k] = 0.0;
	
	for(int i=0; i<sp->polSize; i++){
		pcfg->tuvCMS[0] += pcfg->t[i];
		pcfg->tuvCMS[1] += pcfg->u[i];
		pcfg->tuvCMS[2] += pcfg->v[i];
	}
	for(int k=0; k<3; k++) pcfg->tuvCMS[k] /= sp->polSize;
}

void AddRee(SimProperties* sp, PolyTimeLapse* ptl){
	for(int dt=0; dt<ptl->nEqd; dt++){
		double addRee=0; int num=0;
		for(int t=ptl->nTherm; t<sp->nTime-dt; t++){
			int t2 = t+dt;
			for(int i=0; i<3; i++){
				addRee += ptl->polys[t].ree[i]*ptl->polys[t2].ree[i];
				num++;
			}
		}
		ptl->avgRee[dt] += addRee/num;
	}
}

void AddAvgPos(SimProperties* sp, PolyTimeLapse* ptl){
	double invsqrt2 = 1/sqrt(2);
	for(int t=0; t<sp->nTime; t++){
		PolyConfig* pcfg = ptl->polys+t;
		ptl->avgPosition[t][0] += pcfg->cms.x*invsqrt2;
		ptl->avgPosition[t][1] += pcfg->cms.y*invsqrt2;
		ptl->avgPosition[t][2] += pcfg->cms.z*invsqrt2;
	}
}

void AddShearMod(SimProperties* sp, PolyTimeLapse* ptl){
	int dtStep=1;
	for(int dt=0; dt<ptl->nEqd; dt+=dtStep){
		double addShear=0; int num=0;
		for(int t=ptl->nTherm; t<sp->nTime-dt; t+= MAX(1,1)){
			PolyConfig* pcfg1 = ptl->polys+t;
			PolyConfig* pcfg2 = ptl->polys+t+dt;

			for(int i=0; i<3; i++){
				addShear += pcfg1->stress[i]*pcfg2->stress[i];
				num++;
			}
		}
		ptl->avgShearMod[dt] += addShear/num;
		dtStep = MAX(1, 1);
	}
}

void AddRGyr(SimProperties* sp, PolyTimeLapse* ptl){
	for(int t=0; t<sp->nTime; t++){
		ptl->rGyrT[t] += ptl->polys[t].rGyr;
		if(t>ptl->nTherm)
			ptl->avgRGyr += ptl->polys[t].rGyr;
	}
}

void AddSL(SimProperties* sp, PolyTimeLapse* ptl){
	for(int t=0; t<sp->nTime; t++)
		ptl->avgSL[t] += ptl->polys[t].slRat;
}

void AddGenom(SimProperties* sp, PolyTimeLapse* ptl){
	for(int t=ptl->nTherm; t<sp->nTime; t++){
		for(int i=0; i<sp->polSize; i++){
			ptl->avgGenom[i] += ptl->polys[t].genom[i];
		}
	}
}
/*
void AddGenomNew(SimProperties* sp, PolyTimeLapse* ptl){
	int dx, dy, dz;
	int dr;
	for(int iGenom=0; iGenom<ptl->nGenom; iGenom++){
		int i=ptl->genomList[iGenom].x;
		int j=ptl->genomList[iGenom].y;
		int g = (j-i+sp->polSize)%sp->polSize;
		for(int t=0; t<ptl->nMeas; t++){
			PolyConfig* pcfg = ptl->polys+t;
			
			dx = pcfg->x[i]-pcfg->x[j];
			dy = pcfg->y[i]-pcfg->y[j];
			dz = pcfg->z[i]-pcfg->z[j];
			dr = (dx*dx+dy*dy+dz*dz);
			
			ptl->avgGenom[g] += dr/2.0;
// 			printf("%i %i\n", g, dr);
			if(dr>sp->polSize*50) dr = sp->polSize*50-1;
			ptl->genomProb[g][dr]++;
			ptl->genomSample[g]++;
			
		}
	}
}
*/
void AddGenomNew2(SimProperties* sp, PolyTimeLapse* ptl){
	int dx, dy, dz;
	int dr;
	Timer time;
// 	printf("Doing %i operations\n", ptl->nGenom*ptl->nMeas);
	for(int iGenom=0; iGenom<ptl->nGenom; iGenom++){
		int i=ptl->genomList[iGenom].x;
		int j=ptl->genomList[iGenom].y;
		int g = (j-i+sp->polSize)%sp->polSize;
		TimerStart(&time);
		for(int t=ptl->nTherm; t<sp->nTime; t++){
			PolyConfig* pcfg = ptl->polys+t;
			
			dx = pcfg->x[i]-pcfg->x[j];
			dy = pcfg->y[i]-pcfg->y[j];
			dz = pcfg->z[i]-pcfg->z[j];
			dr = (dx*dx+dy*dy+dz*dz);
			
			ptl->avgGenom[g] += dr/2.0;
			double sqrtdr = ptl->sqrtList[dr];
			int bin = (int)sqrtdr;
			
			ptl->genomProb[g][bin]++;
			ptl->genomR[g][bin]+=sqrtdr;
			ptl->genomCount[g]++;
		}
// 		printf("g=%i, took %.2f ms\n", g, TimerElapsed(&time)*1e3);
	}
// 	printf("done\n");
}


// 	int di=20;
// 	int gMax;
// 			
// 		if(sp->polType == POL_TYPE_LIN) gMax=sp->polSize-1;
// 		else gMax = sp->polSize/2;
// 		for(int g=1, dg=1; g<=gMax; g+=dg){
// 			int num=0;
// 			pcfg->genom[g]=0;
// 			di = MAX(1,MAX(sp->polSize/100, dg/10));
// 			for(int i=0; i<sp->polSize; i+=di){
// 				if(i+g >= sp->polSize && sp->polType == POL_TYPE_LIN) continue;
// 				int j = (i+g)%sp->polSize;
// 				
// 				dx = pcfg->x[i]-pcfg->x[j];
// 				dy = pcfg->y[i]-pcfg->y[j];
// 				dz = pcfg->z[i]-pcfg->z[j];
// 				
// 				pcfg->genom[g] += (dx*dx+dy*dy+dz*dz)/2.0;
// 				num++;
// 			}
// 			
// 			pcfg->genom[g] /= num;
// 			if(!num){
// 				printf("woops: g = %i, dg=%i\n", g, dg);
// 				exit(0);
// 			}
// 			dg = MAX(1,MIN(gMax/2-g, dg/10));
// 		}
// 	}
// }

void AddUnit(SimProperties* sp, PolyTimeLapse* ptl){
	for(int t=ptl->nTherm; t<sp->nTime; t++){
		for(int i=0; i<sp->polSize; i++){
			ptl->avgUnitCor[i] += ptl->polys[t].unitCor[i];
		}
	}
}

///First Statics:
void AddModesStat(SimProperties* sp, PolyTimeLapse* ptl){
	for(int ip=0; ip<ptl->nModes; ip++){
		int p = ptl->modeList[ip];
		for(int iq=ip; iq<ptl->nModes; iq++){
			int q = ptl->modeList[iq];
			double crossAdd=0;
			for(int t=ptl->nTherm; t<sp->nTime; t++){
				crossAdd += ptl->polys[t].sxmode[p]*ptl->polys[t].sxmode[q];
				crossAdd += ptl->polys[t].symode[p]*ptl->polys[t].symode[q];
				crossAdd += ptl->polys[t].szmode[p]*ptl->polys[t].szmode[q];
				crossAdd += ptl->polys[t].cxmode[p]*ptl->polys[t].cxmode[q];
				crossAdd += ptl->polys[t].cymode[p]*ptl->polys[t].cymode[q];
				crossAdd += ptl->polys[t].czmode[p]*ptl->polys[t].czmode[q];
			}
			
			
			ptl->avgModesStat[p][q] += crossAdd/(double)(2*ptl->nEqd);
			if(q != p)
				ptl->avgModesStat[q][p] += crossAdd/(double)(2*ptl->nEqd);
		}
	}
}

///Then Dynamics:

void AddModesDyn(SimProperties* sp, PolyTimeLapse* ptl){
	int dtStep=1;
	PolyConfig* pcfg1, *pcfg2;
	for(int dt=0; dt<ptl->nEqd; dt+=dtStep){
		for(int ip=0; ip<ptl->nModes; ip++){
			int p = ptl->modeList[ip];
			double crossAdd = 0; int num=0;
			int dtMeas = MAX(1, dt/20);
			for(int t=ptl->nTherm; t<sp->nTime-dt; t+= dtMeas){
				int t2 = t+dt;
				pcfg1 = ptl->polys+t;
				pcfg2 = ptl->polys+t2;
					
				crossAdd += pcfg1->sxmode[p]*pcfg2->sxmode[p];
				crossAdd += pcfg1->symode[p]*pcfg2->symode[p];
				crossAdd += pcfg1->szmode[p]*pcfg2->szmode[p];
				crossAdd += pcfg1->cxmode[p]*pcfg2->cxmode[p];
				crossAdd += pcfg1->cymode[p]*pcfg2->cymode[p];
				crossAdd += pcfg1->czmode[p]*pcfg2->czmode[p];
				num++;
			}
			ptl->avgModesDyn[p][dt] += crossAdd/(double)(2*num);
		}
		dtStep = MAX(1, dt/10);
	}
}

void AddContactProbability(SimProperties* sp, PolyTimeLapse* ptl){
// 	int dx, dy, dz;
// 	double dr, pc;
	
	for(int t=ptl->nTherm; t<sp->nTime; t++){
		PolyConfig* pcfg = ptl->polys+t;
		int L = ptl->L;
		LatPoint* lattice = ptl->lattice;
		
		int ts = pcfg->t[0];
		int us = pcfg->u[0];
		int vs = pcfg->v[0];
		
		int t,u,v;
		for(int iMono=0; iMono<sp->polSize; iMono++){
			t = pcfg->t[iMono]-ts;
			u = pcfg->u[iMono]-us;
			v = pcfg->v[iMono]-vs;
			
			int pos = t+u*L+v*L*L;
			
			if(lattice[pos].nOcc)
				ptl->monoList[iMono]=lattice[pos].firstMono;
			else
				ptl->monoList[iMono]=-1;
			lattice[pos].firstMono=iMono;
			for(int curMono=ptl->monoList[iMono]; curMono >=0; curMono=ptl->monoList[curMono])
				ptl->pc[iMono][curMono]++;
			
			for(int i=0; i<12; i++){
				int dt = tuvRelPos[i][0];
				int du = tuvRelPos[i][1];
				int dv = tuvRelPos[i][2];
				
				int newPos = pos+dt+du*L+dv*L*L;
				if(newPos > ptl->L*ptl->L*ptl->L/2 || newPos<-ptl->L*ptl->L*ptl->L/2){
					printf("Oh dear!\n");
					printf("%i %i %i\n", t,u,v);
					printf("%i %i %i\n", dt, du, dv);
					printf("%i %i\n\n", pos, newPos);
					for(int i=0; i<sp->polSize; i++){
						printf("%i %i %i\n", pcfg->t[i]-ts, pcfg->u[i]-us, pcfg->v[i]-vs);
					}
					exit(192);
				}
				if(lattice[newPos].nOcc){
					for(int curMono=lattice[newPos].firstMono; curMono >=0; curMono=ptl->monoList[curMono])
						ptl->pc[iMono][curMono]++;
				}
			}
		}
		
		for(int iMono=0; iMono<sp->polSize; iMono++){
			t = pcfg->t[iMono]-ts;
			u = pcfg->u[iMono]-us;
			v = pcfg->v[iMono]-vs;
			
			int pos = t+u*L+v*L*L;
			lattice[pos].nOcc=0;
		}
	}
}


void AddMono(SimProperties* sp, PolyTimeLapse* ptl, int mono, double* res){
	int dtStep=1, num=0;
	double mAdd=0;
	for(int dt=0; dt<ptl->nEqd; dt+=dtStep){
		num=0; mAdd=0;
		for(int t=ptl->nTherm; t<sp->nTime-dt; t+= MAX(1,dt/10)){
			int t2 = t+dt;
			double dx=0, dy=0, dz=0;
			if(!sp->updAvgPos){
				dx = ptl->avgPosition[t][0]-ptl->avgPosition[t2][0];
				dy = ptl->avgPosition[t][1]-ptl->avgPosition[t2][1];
				dz = ptl->avgPosition[t][2]-ptl->avgPosition[t2][2];
			}

			mAdd += pow((ptl->polys[t].x[mono]-ptl->polys[t2].x[mono]-dx),2);
			mAdd += pow((ptl->polys[t].y[mono]-ptl->polys[t2].y[mono]-dy),2);
			mAdd += pow((ptl->polys[t].z[mono]-ptl->polys[t2].z[mono]-dz),2);
			num++;
		}
		res[dt] += mAdd/num;
		dtStep = MAX(1, dt/10);
	}
}

void AddCMS(SimProperties* sp, PolyTimeLapse* ptl){
	int dtStep=1, num=0;
	double cmsAdd=0;
	for(int dt=0; dt<ptl->nEqd; dt+=dtStep){
		num=0; cmsAdd=0;
		for(int t=ptl->nTherm; t<sp->nTime-dt; t+= MAX(1,dt/10)){
			int t2 = t+dt;
			double dx=0, dy=0, dz=0;
			if(!sp->updAvgPos){
				dx = ptl->avgPosition[t][0]-ptl->avgPosition[t2][0];
				dy = ptl->avgPosition[t][1]-ptl->avgPosition[t2][1];
				dz = ptl->avgPosition[t][2]-ptl->avgPosition[t2][2];
			}
			cmsAdd += pow((ptl->polys[t].cms.x-ptl->polys[t2].cms.x-dx),2);
			cmsAdd += pow((ptl->polys[t].cms.y-ptl->polys[t2].cms.y-dy),2);
			cmsAdd += pow((ptl->polys[t].cms.z-ptl->polys[t2].cms.z-dz),2);
			num++;
		}
		ptl->cmsDif[dt] += cmsAdd/num;
		dtStep = MAX(1, dt/10);
	}
}

void AddMM(SimProperties* sp, PolyTimeLapse* ptl){AddMono(sp, ptl, sp->polSize/2, ptl->mmDif);}
void AddSM(SimProperties* sp, PolyTimeLapse* ptl){AddMono(sp, ptl, 0, ptl->smDif);}
void AddEM(SimProperties* sp, PolyTimeLapse* ptl){AddMono(sp, ptl, sp->polSize-1, ptl->emDif);}

void AddSpaceRouse(SimProperties* sp, PolyTimeLapse* ptl){
	double dr[3];
	double tuv[3];
	
	for(int i=0; i<ptl->tTable.nTDT; i++){
		int dt = ptl->tTable.dt[i];
		int t1 = ptl->tTable.t[i];
		int t2 = dt+t1;
		
		dr[0] = ptl->polys[t1].cms.x-ptl->polys[t2].cms.x;
		dr[1] = ptl->polys[t1].cms.y-ptl->polys[t2].cms.y;
		dr[2] = ptl->polys[t1].cms.z-ptl->polys[t2].cms.z;
		
		for(int k=0; k<3; k++){
			for(int j=0; j<3; j++){
				tuv[j] = ptl->polys[t1].tuvCMS[j];
				for(int p=1; p<sp->polSize/2; p++){
					ptl->sAvgSpacMode[p][i][k][j] += dr[k]*sin(2*PI*p*tuv[j]/(double)sp->LT);
					ptl->cAvgSpacMode[p][i][k][j] += dr[k]*cos(2*PI*p*tuv[j]/(double)sp->LT);
				}
			}
		}
	}
}

void AddSpacDif(SimProperties* sp, PolyTimeLapse* ptl, SpacDif* sd){
	int coor[3];
	double dr[3];
	
	for(int iMono=0; iMono<sp->polSize; iMono++){
		for(int itdt=0; itdt<ptl->tTable.nTDT; itdt++){
			int dt = ptl->tTable.dt[itdt];
			int t1 = ptl->tTable.t[itdt];
			int t2 = t1+dt;
			coor[0] = ptl->polys[t1].t[iMono]; 
			coor[1] = ptl->polys[t1].u[iMono];
			coor[2] = ptl->polys[t1].v[iMono];
			RemPeriod(coor, coor);
			int tuv = coor[0]+coor[1]*sp->LT+coor[2]*sp->LT*sp->LU;
			
			dr[0] = ptl->polys[t2].x[iMono]-ptl->polys[t1].x[iMono];
			dr[1] = ptl->polys[t2].y[iMono]-ptl->polys[t1].y[iMono];
			dr[2] = ptl->polys[t2].z[iMono]-ptl->polys[t1].z[iMono];
			
// 			if(dr[0] > 500 || dr[0]<-500 || dr[1] > 500 || dr[1]<-500 || dr[2] > 500 || dr[2]<-500)
// 				printf("Uh? (%lf %lf %lf)\n", dr[0], dr[1], dr[2]);
// 			if(dt==1 && t1==11) printf("dr=%lf\n", dr[1]);
// 			if(tuv>=LT*LU*LV || tuv<0)
// 				printf("Uh oh..\n");
			
			for(LList* tLink=sd->ballList[tuv]; tLink != NULL; tLink = tLink->next){
// 				printf("hi!\n");
				int sdId = tLink->sdId;
				for(int i=0; i<3; i++){
					ptl->avgSpacDif[sdId][ptl->devId][itdt][i] += dr[i];
					if(iMono==0 && itdt==1 && sdId==0){
						printf("%lf ", ptl->avgSpacDif[sdId][ptl->devId][itdt][i]);
						if(i==2) printf("\n");
					}
				}
			}
		}
	}
}


