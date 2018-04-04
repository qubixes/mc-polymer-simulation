#include "lowm_modes.h"

void ComputeObsPoly(SimProperties* sp, PolyConfig* pcfg){
	ComputeCMS(sp, pcfg);
	ComputeTUVCMS(sp, pcfg);
	if(sp->updRouseStat || sp->updRouseDyn) ComputeModes(sp, pcfg);

	if(sp->updUnit)      ComputeUnitCor(sp, pcfg);
	if(sp->updGyr)       ComputeGyration(sp, pcfg);
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
	curPolId++;

	printf("[Compute: %.2lf ms] ", 1e3*TimerElapsed(&time)); fflush(NULL);
	TimerStart(&time);
	if(sp->updRouseStat) AddModesStat(sp,ptl);
	if(sp->updRouseDyn ) AddModesDyn(sp,ptl);
	if(sp->updPC)        AddContactProbability(sp,ptl);
	if(sp->updGenom)     AddGenom(sp,ptl);
	if(sp->updUnit)      AddUnit(sp,ptl);
	if(sp->updSPRouse)   AddSpaceRouse(sp, ptl);
	if(sp->updSpacDif)   AddSpacDif(sp,ptl,&sd);
	if(sp->updSL)        AddSL(sp,ptl);
	if(sp->updGyr)       AddRGyr(sp,ptl);
	if(sp->updDif){
		AddEM(sp,ptl); AddCMS(sp,ptl); 
		AddSM(sp,ptl); AddMM(sp,ptl); 
	}
	if(sp->updAvgPos) AddAvgPos(sp,ptl);
	if(sp->updShearMod) AddShearMod(sp,ptl);
	if(sp->updRee) AddRee(sp, ptl);
	ptl->nPolAdded++;
	printf("[Add: %.2lf ms]", 1e3*TimerElapsed(&time)); fflush(NULL);
}

void ComputeRee(SimProperties* sp, PolyConfig* pcfg){
	int monoId;
	if(pcfg->polType == POL_TYPE_LIN)
		monoId = pcfg->polSize;
	else
		monoId = pcfg->polSize/2;
	pcfg->ree[0] = (pcfg->x[0]-pcfg->x[monoId])/sqrt(2);
	pcfg->ree[1] = (pcfg->y[0]-pcfg->y[monoId])/sqrt(2);
	pcfg->ree[2] = (pcfg->z[0]-pcfg->z[monoId])/sqrt(2);
}

void ComputeShearMod(SimProperties* sp, PolyConfig* pcfg){
	
	for(int i=0; i<3; i++) pcfg->stress[i]=0;
	
	for(int i=0; i<pcfg->polSize; i++){
		double dx = pcfg->x[(i+1)%pcfg->nMono]-pcfg->x[i];
		double dy = pcfg->y[(i+1)%pcfg->nMono]-pcfg->y[i];
		double dz = pcfg->z[(i+1)%pcfg->nMono]-pcfg->z[i];
		pcfg->stress[0] += dx*dy;
		pcfg->stress[1] += dy*dz;
		pcfg->stress[2] += dx*dz;
// 		printf("[dx dy dz] %lf %lf %lf\n", dx*dy, dy*dz, dx*dz);
	}
	for(int i=0; i<3; i++)
		pcfg->stress[i] /=2;
// 	printf("%lf %lf %lf\n", pcfg->stress[0], pcfg->stress[1], pcfg->stress[2]);
}

void ComputeModes(SimProperties* sp, PolyConfig* pcfg){
	for(int ip=0; ip<pcfg->nModes; ip++){
		int p = pcfg->modeList[ip];
		pcfg->sxmode[p]=0; pcfg->symode[p]=0; pcfg->szmode[p]=0; 
		pcfg->cxmode[p]=0; pcfg->cymode[p]=0; pcfg->czmode[p]=0;
		for(int k=0; k<pcfg->polSize; k++){
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
	for(int i=0; i<pcfg->polSize; i++){
		if(pcfg->x[i] == pcfg->x[(i+1)%pcfg->nMono] && pcfg->y[i] == pcfg->y[(i+1)%pcfg->nMono] && pcfg->z[i] == pcfg->z[(i+1)%pcfg->nMono]) 
			nSl++;
	}
	
	pcfg->slRat = nSl/(double)pcfg->polSize;
}

void ComputeUnitCor(SimProperties* sp, PolyConfig* pcfg){
	int di=5;
	int dxi, dxj, dyi,dyj,dzi,dzj;
	int gMax;
	
	if(pcfg->polType == POL_TYPE_LIN)
		gMax = pcfg->polSize;
	else
		gMax = pcfg->polSize/2;
	
	for(int g=1, dg=1; g<=gMax; g+=dg){
		int num=0;
		pcfg->unitCor[g] = 0;
		di = MAX(1,MAX(pcfg->polSize/50, dg/10));
		for(int i=0; i<pcfg->polSize; i+=di){
			if(i+g+1 >= pcfg->nMono && pcfg->polType == POL_TYPE_LIN) continue;
			int j = (i+g)%pcfg->nMono;
			
			dxi = pcfg->x[i]-pcfg->x[(i+1)%pcfg->nMono];
			dyi = pcfg->y[i]-pcfg->y[(i+1)%pcfg->nMono];
			dzi = pcfg->z[i]-pcfg->z[(i+1)%pcfg->nMono];
			
			dxj = pcfg->x[j]-pcfg->x[(j+1)%pcfg->nMono];
			dyj = pcfg->y[j]-pcfg->y[(j+1)%pcfg->nMono];
			dzj = pcfg->z[j]-pcfg->z[(j+1)%pcfg->nMono];
			
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
	double invSqrt2=1/sqrt(2.0);
	for(int i=0; i<pcfg->nMono; i++){
		dx = pcfg->x[i]*invSqrt2-pcfg->cms.x;
		dy = pcfg->y[i]*invSqrt2-pcfg->cms.y;
		dz = pcfg->z[i]*invSqrt2-pcfg->cms.z;
		rg += (dx*dx+dy*dy+dz*dz)/2;
	}
	rg /= pcfg->nMono;
	pcfg->rGyr = rg;
}

void ComputeCMS(SimProperties* sp, PolyConfig* pcfg){
	pcfg->cms.x=0.0; 
	pcfg->cms.y=0.0; 
	pcfg->cms.z=0.0;
	for(int i=0; i<pcfg->nMono; i++){
		pcfg->cms.x += pcfg->x[i];
		pcfg->cms.y += pcfg->y[i];
		pcfg->cms.z += pcfg->z[i];
	}
	pcfg->cms.x /= pcfg->nMono*sqrt(2);
	pcfg->cms.y /= pcfg->nMono*sqrt(2);
	pcfg->cms.z /= pcfg->nMono*sqrt(2);
}

void ComputeTUVCMS(SimProperties* sp, PolyConfig* pcfg){
	for(int k=0; k<3; k++) pcfg->tuvCMS[k] = 0.0;
	
	for(int i=0; i<pcfg->nMono; i++){
		pcfg->tuvCMS[0] += pcfg->t[i];
		pcfg->tuvCMS[1] += pcfg->u[i];
		pcfg->tuvCMS[2] += pcfg->v[i];
	}
	for(int k=0; k<3; k++) pcfg->tuvCMS[k] /= pcfg->nMono;
}

void AddRee(SimProperties* sp, PolyTimeLapse* ptl){
	for(int dt=0; dt<ptl->nEqd; dt++){
		double addRee=0; int num=0;
		for(int t=ptl->nTherm; t<sp->nTime-dt; t++){
			int t2 = t+dt;
			for(int i=0; i<3; i++){
				addRee += ptl->polys[t].ree[i]*ptl->polys[t2].ree[i];
			}
			num++;
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
	int dx, dy, dz;
	int dr;
	for(int iGenom=0; iGenom<ptl->nGenom; iGenom++){
		int i =ptl->genomList[iGenom].x;
		int j =ptl->genomList[iGenom].y;
		int ig=ptl->genomList[iGenom].ig;
		for(int t=ptl->nTherm; t<sp->nTime; t++){
			PolyConfig* pcfg = ptl->polys+t;
			
			dx = pcfg->x[i]-pcfg->x[j];
			dy = pcfg->y[i]-pcfg->y[j];
			dz = pcfg->z[i]-pcfg->z[j];
			dr = (dx*dx+dy*dy+dz*dz);
			
			ptl->avgGenom[ig] += dr/2.0;
			double sqrtdr = sqrt(dr/2.0);
			int bin = (int) sqrtdr;
			
			bin = MIN(bin, ptl->nGenomBin);
			ptl->genomProb[ig][bin]++;
			ptl->genomR[ig][bin]+=sqrtdr;
			ptl->genomCount[ig]++;
		}
	}
}

void AddUnit(SimProperties* sp, PolyTimeLapse* ptl){
	for(int t=ptl->nTherm; t<sp->nTime; t++){
		for(int i=0; i<ptl->polys[0].polSize; i++){
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
			
			
			ptl->avgModesStat[ip][iq] += crossAdd/(double)(2*ptl->nEqd);
			if(q != p)
				ptl->avgModesStat[iq][ip] += crossAdd/(double)(2*ptl->nEqd);
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
			ptl->avgModesDyn[ip][dt] += crossAdd/(double)(2*num);
		}
		dtStep = MAX(1, dt/10);
	}
}

void GetPosImg(int t, int u, int v, PolyTimeLapse* ptl, int* pos, int* img){
	/// OK, this is confusing, but it's not actually the size of the lattice...
	/// It's a helpful lattice size to figure out if the two are on the same actual site, 
	/// or just through the boundary conditions.
	int L    = ptl->L; 
	int LIMG = ptl->LIMG; 
	
	int posT = (t+L*LIMG)%L;
	int posU = (u+L*LIMG)%L;
	int posV = (v+L*LIMG)%L;
	
	int imgT = (t-posT)/L;
	int imgU = (u-posU)/L;
	int imgV = (v-posV)/L;
	
	*pos = posT + posU*L    + posV*L*L;
	*img = imgT + imgU*LIMG + imgV*LIMG*LIMG;
}

void AddSSContactProbability(SimProperties* sp, PolyTimeLapse* ptl){
	for(int t=ptl->nTherm; t<sp->nTime; t++){
		PolyConfig* pcfg = ptl->polys+t;
		for(int iMono=0; iMono<pcfg->nMono; iMono++){
			int pos, img;
			int dt = pcfg->t[iMono]-pcfg->t[0];
			int du = pcfg->u[iMono]-pcfg->u[0];
			int dv = pcfg->v[iMono]-pcfg->v[0];
			
			GetPosImg(dt, du, dv, ptl, &pos, &img);
			ptl->ssMonoList[iMono].val = pos+img*sp->LSIZE;
			ptl->ssMonoList[iMono].id = iMono;
		}
		qsort(ptl->ssMonoList, pcfg->nMono, sizeof(IDouble), &CompareIDouble);
		
		for(int iMono=0; iMono<pcfg->nMono; ){
			int jMono;
			for(jMono=iMono+1; jMono<pcfg->nMono; jMono++){
				if(ptl->ssMonoList[iMono].val != ptl->ssMonoList[jMono].val) break;
			}
			
			for(int kMono=iMono; kMono<jMono; kMono++){
				for(int lMono=kMono+1; lMono<jMono; lMono++){
					if(lMono==kMono) continue;
					ptl->pcss[ptl->ssMonoList[kMono].id/ptl->pcBins][ptl->ssMonoList[lMono].id/ptl->pcBins]++;
					ptl->pcssAvg[abs(ptl->ssMonoList[kMono].id-ptl->ssMonoList[lMono].id)]++;
				}
			}
			iMono=jMono;
		}
	}
}

void AddContactProbability(SimProperties* sp, PolyTimeLapse* ptl){
// 	int dx, dy, dz;
// 	double dr, pc;
	
	AddSSContactProbability(sp, ptl);
	
	for(int t=ptl->nTherm; t<sp->nTime; t++){
		PolyConfig* pcfg = ptl->polys+t;
// 		int L = ptl->L;
// 		int LIMG= ptl->LIMG;
		LatPoint* lattice = ptl->lattice;
		
		int ts = pcfg->t[0];
		int us = pcfg->u[0];
		int vs = pcfg->v[0];
		
		int t,u,v;
		for(int iMono=0; iMono<pcfg->nMono; iMono++){
			t = pcfg->t[iMono]-ts;
			u = pcfg->u[iMono]-us;
			v = pcfg->v[iMono]-vs;
			
// 			int pos = t+u*L+v*L*L;
			int pos, img;
			GetPosImg(t,u,v, ptl, &pos, &img);
			
			if(lattice[pos].nOcc)
				ptl->monoList[iMono].next = lattice[pos].firstMono;
			else
				ptl->monoList[iMono].next = -1;
			ptl->monoList[iMono].img = img;
			ptl->monoList[iMono].pos = pos;
			lattice[pos].firstMono=iMono;
			lattice[pos].nOcc++;
			
			for(int curMono=ptl->monoList[iMono].next; curMono >=0; curMono=ptl->monoList[curMono].next){
				if(iMono == curMono){
					printf("????????????\n");
					exit(192);
				}
				if(ptl->monoList[curMono].img == img){
					ptl->pc[iMono/ptl->pcBins][curMono/ptl->pcBins]++;
					ptl->pcAvg[abs(iMono-curMono)]++;
				}
			}
			
			for(int i=0; i<12; i++){
				int dt = tuvRelPos[i][0];
				int du = tuvRelPos[i][1];
				int dv = tuvRelPos[i][2];
				
				int newPos, newImg;
				GetPosImg(t+dt, u+du, v+dv, ptl, &newPos, &newImg);
				
				if(lattice[newPos].nOcc){
					for(int curMono=lattice[newPos].firstMono; curMono >=0; curMono=ptl->monoList[curMono].next){
// 						printf("%i %i %i\n", iMono, pos, newPos);
						if(ptl->monoList[curMono].img == newImg){
							ptl->pc[iMono/ptl->pcBins][curMono/ptl->pcBins]++;
							ptl->pcAvg[abs(iMono-curMono)]++;
							if(iMono == curMono){
								printf("?????\n");
								exit(192);
							}
						}
					}
				}
			}
		}
// 		exit(0);
		for(int iMono=0; iMono<ptl->polys[0].nMono; iMono++){
			lattice[ptl->monoList[iMono].pos].nOcc=0;
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
		res[dt] += mAdd/(double)(2*num);
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

void AddMM(SimProperties* sp, PolyTimeLapse* ptl){AddMono(sp, ptl, ptl->polys[0].polSize/2, ptl->mmDif);}
void AddSM(SimProperties* sp, PolyTimeLapse* ptl){AddMono(sp, ptl, 0, ptl->smDif);}
void AddEM(SimProperties* sp, PolyTimeLapse* ptl){AddMono(sp, ptl, ptl->polys[0].polSize-1, ptl->emDif);}

void AddSpaceRouse(SimProperties* sp, PolyTimeLapse* ptl){
	double dr[3];
	double tuv[3];
	
	for(int i=0; i<ptl->tTable.nTDT; i++){
		int dt = ptl->tTable.tdt[i].dt;
		int t1 = ptl->tTable.tdt[i].t;
		int t2 = dt+t1;
		
		dr[0] = ptl->polys[t1].cms.x-ptl->polys[t2].cms.x;
		dr[1] = ptl->polys[t1].cms.y-ptl->polys[t2].cms.y;
		dr[2] = ptl->polys[t1].cms.z-ptl->polys[t2].cms.z;
		
		for(int k=0; k<3; k++){
			for(int j=0; j<3; j++){
				tuv[j] = ptl->polys[t1].tuvCMS[j];
				for(int p=1; p<ptl->polys[0].polSize/2; p++){
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
	double invSqrt2 = 1/sqrt(2.0);
	
	for(int iMono=0; iMono<ptl->polys[0].nMono; iMono++){
		for(int itdt=0; itdt<ptl->tTable.nTDT; itdt++){
			int dt = ptl->tTable.tdt[itdt].dt;
			int t1 = ptl->tTable.tdt[itdt].t;
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
					ptl->avgSpacDif[sdId][ptl->devId][itdt][i] += dr[i]*invSqrt2;
// 					if(iMono==0 && itdt==1 && sdId==0){
// 						printf("%lf ", ptl->avgSpacDif[sdId][ptl->devId][itdt][i]);
// 						if(i==2) printf("\n");
// 					}
				}
			}
		}
	}
}


