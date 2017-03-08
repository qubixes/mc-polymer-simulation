#include "gpupol_init.h"

int SimStateInit(SimState* ss, SimProperties* sp){
	SimState* curState;
	uint* wgSeed = (uint*) malloc(sizeof(uint)*sp->R);
	uint globalWs = sp->nwg*sp->ws;
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		curState = ss+iDev;
		curState->lattice = (char*) malloc(sp->latSize*sizeof(char));
		curState->gpuLattice = (char*) malloc(sp->latSize*sizeof(char));
// 		printf("labSize = %i\n", sp->labSize);
		curState->nPol=0;
		curState->rngState.x = sp->seed^(85127698+iDev*716576218);
		curState->rngState.y = sp->seed^(12369123+iDev*612737112);
		curState->rngState.z = sp->seed^(81261212+iDev*172635798);
		curState->rngState.w = sp->seed^(51203016+iDev*410287361);
		curState->labelHead = 0x10;
		curState->headBit = 5;
		curState->labelBit = 18;
		ClearLattice(curState, sp);
		curState->seeds = (uint*) malloc(sizeof(uint)*globalWs*2*sp->R);
		for(int wid=0; wid<sp->nwg; wid++){
			for(int i=0; i<sp->R; i++) wgSeed[i] = Rng4(&(curState->rngState));
			
			for(int lid=0; lid<sp->ws; lid++){
				uint gIndex = 2*sp->R*(lid+wid*sp->ws);
				for(int i=0; i<sp->R; i++){
					curState->seeds[gIndex+i] = Rng4(&(curState->rngState));
					curState->seeds[gIndex+i+sp->R] = wgSeed[i];
				}
			}
		}
	}
	free(wgSeed);
	return 0;
}

void ClearLattice(SimState* ss, SimProperties* sp){
	for(int i=0; i<sp->latSize; i++) ss->lattice[i]=0;
}

void SetBoxDimension(SimProperties* sp, int LT, int LU, int LV){
	if(LT%(LCELL*WST) || LU%(LCELL*WSU) || LV%(LCELL*WSV))
		printf("Warning changing box dimensions to be a multiple of 24!\n");
	
	int nwt = LT/(LCELL*WST);
	int nwu = LU/(LCELL*WSU);
	int nwv = LV/(LCELL*WSV);
	
	sp->nwt = nwt;
	sp->nwu = nwu;
	sp->nwv = nwv;
	
	sp->nwg = sp->nwt*sp->nwu*sp->nwv;
	sp->LT  = LCELL*sp->nwt*WST;
	sp->LU  = LCELL*sp->nwu*WSU;
	sp->LV  = LCELL*sp->nwv*WSV;
	sp->L   = sp->LT*sp->LU*sp->LV;
	sp->latSize = sp->L;
	sp->maxPolLength=sp->LT*sp->LU*4;
}

int SimPropDefaultInit(SimProperties* sp){
	sp->tMax=50000;
	sp->R = 4;
	sp->writeInterval = 100;
	sp->polLength = 2000;
	sp->fastEq=0;
	sp->density = 1.0;
	sp->curT= 0;
	sp->ws = WS;
	sp->seed = 239864012;
	sp->dirOut = (char*) malloc(sizeof(char)*1000);
	sp->trans = (uint*) malloc(4*sizeof(uint));
	sp->nSteps = (int)(0.8*CELL_SIZE);
	sp->equilibrated=0;
	SetGPUTrans(sp->trans);
	return 0;
}

void LoadPolymers(SimProperties* sp, SimState* ss){
	char exec[1000];
	char buf[1000];
	char fileOut[1000];
	FILE* pPipe;
	
	GPULibInit(sp, devices, &gpuContext);
	sprintf(exec, "mkdir -p %s", sp->dirOut);
	system(exec);
	
	sprintf(exec, "ls %s", sp->dirOut);
	printf("exec = %s\n", exec);
	pPipe = popen(exec, "r");
	long highestT=-1;
	int nHighestT=-1;
	char* highTFiles[N_MAX_GPU];
	long tFile;
	while(!feof(pPipe)){
		if(fscanf(pPipe, "%s", buf)>0){
			if(strncmp(buf, "t=", 2)) continue;
// 			printf("Found file %s\n", buf);
			tFile = GetTFromName(buf);
// 			printf("t = %i\n", tFile);
			if(tFile>highestT){
				for(int i=0; i<nHighestT; i++) free(highTFiles[i]);
				nHighestT=0;
			}
			if(tFile>=highestT){
				highTFiles[nHighestT++] = (char*) malloc(sizeof(char)*(strlen(buf)+1));
				highestT = tFile;
				strcpy(highTFiles[nHighestT-1], buf);
			}
		}
	}
	pclose(pPipe);
// 	printf("nHighestT = %i\n", nHighestT);
	if(nHighestT == sp->nDevices){
		printf("Starting from t=%li\n", highestT);
		sp->curT = highestT;
		sp->tMax += highestT;
		sp->seed += highestT/sp->writeInterval;
		///Simstate is initialized after reading the dimensions of the box.
		for(int i=0; i<sp->nDevices; i++){
			sprintf(buf, "%s/%s", sp->dirOut, highTFiles[i]);
			printf("file = %s\n", buf);
			ReadLatticeFile(ss+i, sp, buf);
			UpdatePolymerWithLabels(ss+i);
// 			UpdatePolymerWithLabels(ss+i);
// 			printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n\n\n\n");
		}
	}
	else{
		highestT=0;
		SetBoxDimension(sp, sp->LT, sp->LU, sp->LV);
		SimStateInit(ss, sp);
		GetBestPolymerPartitioning(sp, ss, sp->polLength, sp->density);
	}
	
	
	if(nHighestT != sp->nDevices){
		RedistribSL(ss, sp);

		for(int i=0; i<sp->nDevices; i++){
			sprintf(fileOut,"%s/t=0_dev=%i.res", sp->dirOut, i);
			WriteLatticeFile(sp, ss+ i, fileOut);
		}
	}
	GetAllRingPolymers();
}

void FindBestBoxSize(SimProperties* sp, int nPol, int polSize, int* dRes){
	int LT=sp->LT, LU=sp->LU, LV=sp->LV;
	double idealSL=0.3;
	
	double bestScore=-9999999;
	
	for(int dt=2; dt<LT; dt+=2){
		int nt = LT/dt;
		for(int du=1; du<LU; du++){
			int nu=LU/du;
			for(int dv=2; dv<LV; dv+=2){
				int nv=LV/dv;
				if(dt*du*dv*2 > polSize && nt*nu*nv>=nPol){
					double sl=(polSize-dt*du*dv)/(double)polSize;
					double score= (dt*du*dv)/(double)(dt+du+dv) - pow(dt*du*dv, 2./3.)*pow(sl-idealSL,2);
// 					printf("(%i %i %i), (%i %i %i), score=%le=x-%le*%le\n", dt, du, dv, nt, nu, nv, score, pow(dt*du*dv, 2./3.),pow(sl-idealSL,2));
					if(score>bestScore){
						dRes[0] = dt;
						dRes[1] = du;
						dRes[2] = dv;
						bestScore=score;
					}
				}
			}
		}
	}
}


void GetBestPolymerPartitioning(SimProperties* sp, SimState* ss, int polLength, double density){
	uint dt, du, dv, nt, nu, nv;
	double nDoublePol;
	double nAverage, dPol, dDensity;
	int nPol;
	int nPolAdded;
	uint sl;
	
	nDoublePol = density*sp->L/(double)polLength;
	dPol = fabs(nDoublePol-(int)(nDoublePol+0.5));
	nPol = (int)(nDoublePol+0.5);
	dDensity = dPol*polLength/(double)sp->L;
	if( dDensity/density > 1e-5 ) 
		printf("Warning: density is different from the wanted one:\n difference: %lf, density: %lf\n", dDensity, density);
	
	nAverage = pow(nPol, 1./3.)+0.5;
	dt = ((uint)(sp->LT/nAverage)/2)*2;
	du =  (uint)(sp->LU/nAverage);
	dv = ((uint)(sp->LV/nAverage)/2)*2;
	
	int dRes[3];
	FindBestBoxSize(sp, nPol, polLength, dRes);
	dt=dRes[0]; du=dRes[1]; dv=dRes[2];
	
	printf("dt,du,dv = (%i,%i,%i)\n", dt,du,dv);
	nt = sp->LT/dt;
	nu = sp->LU/du;
	nv = sp->LV/dv;
	
	printf("N=(%u,%u,%u) = %u >= %i\n", nt, nu, nv, nt*nu*nv, nPol);
	if(nt*nu*nv<nPol) {
		printf("Error: not enough space to put polymer into\n");
		exit(192);
	}
	printf("Stored length density = %lf\n", (polLength-(dt*du*dv))/(double)polLength);
	
	sl = polLength-(dt*du*dv);
	if(sl>=polLength){
		printf("Internal error: Failed to find polymer partition\n");
		exit(0);
	}
// 	if(polLength<100){
// 		for(int iDev=0; iDev<sp->nDevices; iDev++){
// 			for(int i=0; i<nPol; i++)
// 				ConstructBarRingPolymer(ss+iDev, sp, i, polLength);
// 		}
// 	}
// 	else{
		nPolAdded=0;
		for(uint t=0; t<nt*dt && nPolAdded<nPol; t += sp->LT/nt){
			for(uint u=0; u<nu*du  && nPolAdded<nPol; u += sp->LU/nu){
				for(uint v=0; v<nv*dv && nPolAdded<nPol; v += sp->LV/nv){
					for(int iDev=0; iDev<sp->nDevices; iDev++){
// 						printf("Construction at (%i,%i,%i)\n", t, u, v);
						ConstructCubeRingPolymer(ss+iDev, sp, t,u,v, dt,du,dv,sl);
					}
					nPolAdded++;
				}
			}
		}
// 	}
	sp->maxPolLength = polLength+1;
	GetAllRingPolymers();
	for(int i=0; i<sp->nDevices; i++) UpdatePolymerWithLabels(ss+i);
	printf("------------------------------------------\n");
}

void ConstructBarRingPolymer(SimState* ss, SimProperties* sp, int polId, int length){
	int nt = sp->LT/((length+1)/2-1)+1;
	int nu = sp->LU/2;
	int nv = sp->LV;
	
	int dt = sp->LT/nt;
	int du = 2;
	int dv = 1;
	
	if(polId >= nt*nu*nv){
		printf("Error constructing polymers: Too little space!\n");
		exit(120);
	}
	
// 	if(2*(dt+1)>length) dt = length/2-1;
	
	int sl = length - dt*du*dv;
	int t  = (polId%nt)     *dt;
	int u  = ((polId/nt)%nu)*du;
	int v  = (polId/(nt*nu))*dv;
	
	uint* bonds, *slBonds;
	bonds = (uint*) malloc(sizeof(uint)*length);
	slBonds = (uint*) malloc(sizeof(uint)*length);
	
	int iMono=0;
	for(int i=0; i<dt-1; i++) bonds[iMono++] = 0x1;
	bonds[iMono++] = 0x2;
	for(int i=0; i<dt-1; i++) bonds[iMono++] = (~0x1)&0xf;
#if POL_TYPE == POL_RING
	bonds[iMono++] = (~0x2)&0xf;
#elif POL_TYPE == POL_LIN
	bonds[iMono++] = 0xf;
#endif
	
	for(int i=0; i<sl; i++){
		slBonds[2*i] = 0x0;
		slBonds[2*i+1] = bonds[i];
	}
// 	for(int i=0; i<length; i++) printf("%x", bonds[i]);
// 	printf("\n"); printf("nt=%i\n", nt); exit(0);
	
	for(int i=2*sl; i<length; i++)
		slBonds[i] = bonds[i-sl];
	AddPolymerToSim(ss, sp, t,u,v, length, slBonds);
	free(bonds); free(slBonds);
}

void ConstructCubeRingPolymer(SimState* ss, SimProperties* sp, uint tStart, uint uStart, uint vStart, uint dt, uint du, uint dv, uint sl){
	
	uint dtMax=dt/2;
	uint duMax=du;
	uint dvMax=dv;
	uint tUnit=0x1, uUnit=0x2, vUnit=0x4;
	uint length=0, slLength=0;
	uint t,u,v,i,iDest,iSrc;
// 	uint nFreq[16];
	uint* bonds = (uint*) malloc(2*dt*du*dv*sizeof(uint));
	uint* slBonds = (uint*) malloc(2*dt*du*dv*sizeof(uint));
	for(u=0; u<duMax; u++){
		for(t=0; t<dtMax-1; t++){
			bonds[length++] = tUnit;
		}
		if(u<duMax-1){
			bonds[length++] = uUnit;
			tUnit = (~tUnit)&0xf;
		}
	}
	
		
	bonds[length++] = vUnit;
	
	for(i=0; i<length-1; i++){
		bonds[length+i] = (~bonds[length-i-2])&0xf;
	}
	length += length-1;
	bonds[length++] = vUnit;
	
	for(v=0; v<dvMax; v+=2){
		memcpy(bonds+(v/2)*length, bonds, sizeof(uint)*length);
	}
	length = length*(dvMax/2)-1;
	bonds[length++] = 0xe;

	
	///Keep t, reflect u, v;
	for(i=0; i<length-1; i++){
		iSrc = length-i-2;
		iDest = length+i;
		if(bonds[iSrc] == 0x1 || bonds[iSrc] == 0xe)
			bonds[iDest] = bonds[iSrc];
		else
			bonds[iDest] = (~bonds[iSrc])&0xf;
	}
	length += length-1;
#if POL_TYPE == POL_RING
	bonds[length++] = 0x1;
#elif POL_TYPE == POL_LIN
	bonds[length++] = 0xf;
#endif

	for(i=0; i<sl; i++){
		slBonds[slLength++] = 0x0;
		slBonds[slLength++] = bonds[i];
	}
	for(;i<length; i++)
		slBonds[slLength++] = bonds[i];
	
	if(sl+dt*du*dv != slLength) {
		printf("Internal error: Failed to create ring polymer, length not equal to its desired length (%i vs %i vs %i vs %i)\n", dt*du*dv+sl, slLength, length, sl);
		exit(0);
	}
	
	AddPolymerToSim(ss, sp, tStart, uStart, vStart, slLength, slBonds);
	free(bonds); free(slBonds);
}

int ConstructRingPolymer(SimState* ss, SimProperties* sp, uint t, uint u, uint v){
	int length=0, i;
	uint* bonds, *slBonds;
	
	uint uMax, tMax;
	uint maxPolLength=sp->LT*sp->LU*4;
	bonds = (uint*) malloc(sizeof(uint)*maxPolLength);
	slBonds = (uint*) malloc(sizeof(uint)*maxPolLength);
	uint ib=0;
	uMax = sp->LU-3; tMax=sp->LT-2;
	
	uint tMin=0;
	while(u<uMax){
		while(t>0+tMin && t<tMax+tMin){
			if(tMin)
				bonds[ib++] = 0xe;
			else
				bonds[ib++] = 0x1;
			AddUnitTUV(bonds[ib-1], &t,&u,&v);
		}
		bonds[ib++] = 0x2;
		AddUnitTUV(bonds[ib-1], &t, &u, &v);
		tMin ^= 0x1;
	}
	bonds[ib++] = 0x4;
	while(u>1){
		while(t>0+tMin && t<tMax+tMin){
			if(tMin)
				bonds[ib++] = 0xe;
			else 
				bonds[ib++] = 0x1;
			AddUnitTUV(bonds[ib-1], &t,&u,&v);
		}
		bonds[ib++]=0xd;
		AddUnitTUV(bonds[ib-1], &t,&u,&v);
		tMin ^= 0x1;
	}
	bonds[ib++]=0xb;
	AddUnitTUV(bonds[ib-1], &t,&u,&v);
	for(i=0; i<ib; i++){
		slBonds[length++]=bonds[i];
	}
	printf("length=%i\n", length);
	AddPolymerToSim(ss,sp, t,u,v,length, slBonds);
	free(bonds); free(slBonds);
	return length;
}

int UniqueLabel(int startLabel, int nLabelBit, int head, int headBit){
	int totLabel, headMask;
	
	headMask = (1<<headBit)-1;
	totLabel = head | (startLabel<<headBit) | (head<<(headBit+nLabelBit));
	
	for(int i=1; i<headBit+nLabelBit; i++){
		if(((totLabel>>i)&headMask)==head) return 0;
	}
	return 1;
}

void AddPolymerToSim(SimState* ss, SimProperties *sp, uint t, uint u, uint v, int length, uint* bonds){
	Polymer* pol;
	if(!ss->nPol){
		ss->pol = (Polymer*) malloc(sizeof(Polymer));
		ss->polNumTranslate = (int*) malloc(sizeof(int));
		ss->label2Id = (int*) malloc((1<<ss->labelBit)*sizeof(int));
		ss->id2Label = (int*) malloc((1<<ss->labelBit)*sizeof(int));
		for(int i=0; i<(1<<18); i++){
			ss->label2Id[i]=-1;
			ss->id2Label[i]=-1;
		}
		pol = ss->pol;
	}
	else{
		ss->pol = (Polymer*) realloc(ss->pol, sizeof(Polymer)*(ss->nPol+1));
		ss->polNumTranslate = (int*) realloc(ss->polNumTranslate, sizeof(int)*(ss->nPol+1));
		pol = ss->pol+ss->nPol;
	}
	ss->nPol++;
	
	pol->label = (uint*) malloc(sizeof(uint)*sp->maxPolLength);
	pol->bonds = (uint*) malloc(sizeof(uint)*sp->maxPolLength);
// 	pol->polCoor = (Coor*) malloc(sizeof(Coor)*sp->maxPolLength);
	pol->labelStart = 0;
	
	int startLabel;
	if(!(ss->nPol-1)){
		startLabel=0;
	}
	else{
		startLabel=ss->id2Label[ss->nPol-2]+1;
	}
	while(!UniqueLabel(startLabel, ss->labelBit,ss->labelHead, ss->headBit)) startLabel++;
	if(startLabel> (1<<ss->labelBit)){
		printf("Ran out of labels, please make a smaller box/longer polymers!\n");
	}
	
// 		pol->label[0]=1;

// 		for(int i=1; i<ss->headBit; i++) pol->label[i] = 0;
// 	for(int i=ss->headBit; i<ss->headBit+ss->labelBit; i++) 
// 		pol->label[i] = (startLabel>>(ss->headBit+ss->labelBit-1-i))&0x1;
	
	for(int i=0; i<ss->labelBit; i++) 
		pol->label[i] = (startLabel>>(ss->labelBit-1-i))&0x1;
	for(int i=ss->labelBit, curLab=1; i<length-ss->headBit; i++, curLab^=0x1) 
		pol->label[i] = curLab;
	pol->label[length-ss->headBit]=1;
	for(int i=length-ss->headBit+1; i<length; i++) pol->label[i] = 0;
	ss->id2Label[ss->nPol-1] = startLabel;
	ss->label2Id[startLabel] = ss->nPol-1;
// 	printf("polId=%i, label=0x%x\n", ss->nPol-1, startLabel);
	
// 	for(int i=ss->labelBit+ss->headBit, curLab=1; i<length; i++, curLab^=0x1) pol->label[i] = curLab;
	
// 	pol->label[0]=1;
// 	polId = ss->nPol-1;
// 	for(int i=1; i<sp->maxPolLength; i++) pol->label[i]=0;
// 	
// 	for(int i=2; polId>0; i+=2){
// 		pol->label[i] = polId&0x1;
// 		polId >>= 1;
// 	}
	
	SetBondVecs(ss->lattice, t, u, v, bonds, length, pol->label, sp);
	memcpy(pol->bonds, bonds, sizeof(uint)*length);
	pol->length = length;
	ss->polNumTranslate[ss->nPol-1] = ss->nPol-1;
	pol->startTUV.t = t;
	pol->startTUV.u = u;
	pol->startTUV.v = v;
}

void RedistribSL(SimState* ss, SimProperties* sp){
	printf("Redistributing stored length distribution\n");
	uint* bonds = (uint*) malloc(sizeof(uint)*sp->maxPolLength);
	uint* slDist = (uint*) malloc(sizeof(uint)*sp->maxPolLength);
	int nBond, nSl;
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		ClearLattice(ss+iDev, sp);
		for(int iPol=0; iPol<ss[iDev].nPol; iPol++){
			nBond=0; nSl=0;
			Polymer* pol = ss[iDev].pol+iPol;
// 			PrintBonds(pol->bonds, pol->length);
			for(int i=0; i<pol->length; i++){
				if(pol->bonds[i] == 0x0){
					nSl++;
				}
				else if (pol->bonds[i] != 0xf){
					nBond++;
				}
			}
			
			for (int i=0; i<nBond; i++){
				int j = Rng4(&(ss->rngState)) % (i+1);
				slDist[i] = slDist[j];
				slDist[j] = (i<nSl)?1:0;
			}
			
// 			for(int i=0; i<nBond; i++){printf("%u", slDist[i]);} ;printf("\n");
			
			for(int i=0, id=0, is=0; i<nBond; i++, id++,is++){
				if(pol->bonds[is] == 0x0) is++;
				if(slDist[i] == 1){
					bonds[id++] = 0x0;
				}
				bonds[id] = pol->bonds[is];
			}
			
			memcpy(pol->bonds, bonds, sizeof(uint)*(nBond+nSl));
			uint t = pol->startTUV.t;
			uint u = pol->startTUV.u;
			uint v = pol->startTUV.v;
// 			PrintBonds(pol->bonds, pol->length);
// 			printf("(%u %u %u)\n", t,u,v);
// 			printf("iPol = %i\n", iPol);
			SetBondVecs(ss[iDev].lattice, t,u,v,pol->bonds,pol->length,pol->label, sp);
// 			exit(0);
		}
	}
	printf("Finished redistributing stored length\n");
}

int SetGPUTrans(uint* trans){
	
	uint unit1, unit2;
	uint sum;
	uint nFound[16];
	uint i, bit, block;
	
	for(int i=0; i<16; i++) nFound[i]=0;
	for(int i=0; i<4;  i++) trans[i]=0;
	
	for(unit1=1; unit1<15; unit1++){
		if(!IsValid(unit1)) continue;
		for(unit2=unit1+1; unit2<15; unit2++){
			if(!IsValid(unit2)) continue;
			if(!ValidateAddUnitVectors(unit1, unit2, &sum))
				continue;
			bit = sum%4;
			block = sum/4;
			trans[block] |= unit1<<(8*bit+4*nFound[sum]);
			nFound[sum]++;
			if(nFound[sum]>2){ 
				printf("Internal error: setting trans reactions failed: nFound[0x%x] = %i\n", sum, nFound[sum]);
				exit(0);
			}
		}
	}
	for(i=1; i<15; i++){
		if(i != 0x3 && i != 0xc && nFound[i] != 2){
			printf("Error %i (%i)\n", i, nFound[i]);
			exit(0);
		}
	}
	return 0;
}
