#include "gpupol_cuda.h"
#include "gpupol.h"
#include "gpupol_cuda_def.h"
#include "ringpol_cu.h"

int GPULibInit(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* cudaContext){
	printf("Initializing GPU\n");
	gpuErrchk( cudaGetDeviceCount(&(sp->nDevices)) );
	
	double* perf = (double*)malloc(sizeof(double)*sp->nDevices);
	cudaDeviceProp prop;
	double maxPerf=0;
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		gpuErrchk( cudaGetDeviceProperties(&prop, iDev) );
		perf[iDev] = prop.clockRate*prop.multiProcessorCount;
		if(perf[iDev]>maxPerf) 
			maxPerf=perf[iDev];
	}
	
	int nDev=0;
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		if(perf[iDev]/maxPerf > 0.5)
			cudaContext->devIds[nDev++]=iDev;
	}
	sp->nDevices = nDev;
	free(perf);
	printf("Finished initializing GPU\n");
	return 0;
}

int GPULoadLattices(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext){
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		GPUDeviceState* curDev = devStates+iDev;
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaMemcpy(curDev->latBuf[0],   ss[iDev].lattice, sizeof(char)*sp->latSize,   cudaMemcpyHostToDevice) );
		
		curDev->dt = WLT*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
		curDev->du = WLU*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
		curDev->dv = WLV*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
		
		curDev->cumDt = curDev->dt;
		curDev->cumDu = curDev->du;
		curDev->cumDv = curDev->dv;
		curDev->curBuf=0;
		
		int tuvOff = curDev->dt + curDev->du*WLT+curDev->dv*WLT*WLU;
		polmove <<< sp->nwg, sp->ws, 0 >>> (0, curDev->seedBuf, curDev->latBuf[curDev->curBuf], curDev->latBuf[curDev->curBuf^0x1], curDev->transBuf, 0, tuvOff, sp->nwt, sp->nwu, sp->nwv);
		curDev->curBuf ^= 0x1;		
	}
	return 0;
}

int GPULatticeToCPU(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext){
	
	gpuErrchk( cudaPeekAtLastError() );
	
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		GPUDeviceState* curDev = devStates+iDev;
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaMemcpy(ss[iDev].gpuLattice, curDev->latBuf[curDev->curBuf], sizeof(char)*sp->latSize, cudaMemcpyDeviceToHost) );
	}
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaDeviceSynchronize() );
	}
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		CopyGPUToCPULattice(ss[iDev].gpuLattice, ss[iDev].lattice, devStates[iDev].cumDt, devStates[iDev].cumDu, devStates[iDev].cumDv, devStates[iDev].dt, devStates[iDev].du, devStates[iDev].dv, sp);
	}

	GetAllRingPolymers();
	for(int iDev=0; iDev<sp->nDevices; iDev++)
		UpdatePolymerWithLabels(ss+iDev);
	return 0;	
}

int GPULibLoadBuffers(SimProperties* sp, SimState* ss, GPUDeviceState * devStates, GPULibContext* cudaContext){
	uint globalWs;
	GPUDeviceState* curDev;
	printf("Loading GPU buffers\n");
	globalWs = sp->nwg*sp->ws;
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		curDev = devStates + iDev;
		curDev->latBuf = (char**) malloc(sizeof(char*)*2);
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		
		gpuErrchk( cudaMalloc(&(curDev->seedBuf),  sizeof(uint)*globalWs*2*sp->R) );
		gpuErrchk( cudaMalloc(&(curDev->latBuf[0]),   sizeof(char)*sp->latSize) );
		gpuErrchk( cudaMalloc(&(curDev->latBuf[1]),   sizeof(char)*sp->latSize) );
		gpuErrchk( cudaMalloc(&(curDev->transBuf), sizeof(uint)*4) );
		
		gpuErrchk( cudaMemcpy(curDev->seedBuf,  ss[iDev].seeds,      sizeof(uint)*globalWs*2*sp->R, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(curDev->transBuf, sp->trans,           sizeof(uint)*4,                cudaMemcpyHostToDevice) );
	}
	GPULoadLattices(sp, ss, devStates, cudaContext);
	printf("Finished loading buffers\n");
	return 0;
}

int GPULibRelease(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* cudaContext){
	GPUDeviceState* curDev;
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		curDev = devStates+iDev;
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaFree(curDev->seedBuf) );
		gpuErrchk( cudaFree(curDev->latBuf[0]) );
		gpuErrchk( cudaFree(curDev->latBuf[1]) );
		gpuErrchk( cudaFree(curDev->transBuf) );
	}
	return 0;
}

int GPULibRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext, int nTime){
	int iDev;
	int* tuvOff, *prevTuvOff;
	uint NWT, NWU, NWV;
	GPUDeviceState* curDev;
	
	tuvOff = (int*) malloc(sizeof(int)*sp->nDevices);
	prevTuvOff = (int*) malloc(sizeof(int)*sp->nDevices);
	NWT=sp->nwt; NWU=sp->nwu; NWV=sp->nwv;
	
	for(iDev=0; iDev<sp->nDevices; iDev++){
		tuvOff[iDev] = devStates[iDev].dt+devStates[iDev].du*WLT+devStates[iDev].dv*WLT*WLU;
	}
	
	for(int i=0; i<nTime; i++){
		for(iDev=0; iDev<sp->nDevices; iDev++){
			curDev = devStates+iDev;
			cudaSetDevice(cudaContext->devIds[iDev]);
			
			curDev->dt = WLT*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			curDev->du = WLU*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			curDev->dv = WLV*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			
			curDev->cumDt += curDev->dt;
			curDev->cumDu += curDev->du;
			curDev->cumDv += curDev->dv;
			
			prevTuvOff[iDev] = tuvOff[iDev];
			tuvOff[iDev] = curDev->dt + curDev->du*LCELL*WST + curDev->dv*LCELL*WSU*LCELL*WST;
			
			polmove <<< sp->nwg, sp->ws, 0 >>> (sp->nSteps, curDev->seedBuf, curDev->latBuf[curDev->curBuf], curDev->latBuf[curDev->curBuf^0x1], curDev->transBuf, prevTuvOff[iDev], tuvOff[iDev], NWT, NWU, NWV);
		}
		curDev->curBuf ^= 0x1;
	}
	
	GPULatticeToCPU(sp, ss, devStates, cudaContext);
	
	return 0;
}

int GPULibEqRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext, int nTime){
	int iDev;
	int* tuvOff, *prevTuvOff;
	int nEq = nTime/sp->fastEq;
	int dt = sp->fastEq;
	uint NWT, NWU, NWV;
	GPUDeviceState* curDev;
	
	tuvOff = (int*) malloc(sizeof(int)*sp->nDevices);
	prevTuvOff = (int*) malloc(sizeof(int)*sp->nDevices);
	NWT=sp->nwt; NWU=sp->nwu; NWV=sp->nwv;
	
	for(long t=0; t<nTime; t+=dt){
		long curDt = dt + (((t/dt)<(nEq%dt))?1:0);
		for(iDev=0; iDev<sp->nDevices; iDev++){
			tuvOff[iDev] = devStates[iDev].dt+devStates[iDev].du*WLT+devStates[iDev].dv*WLT*WLU;
		}

		for(int i=0; i<curDt; i++){
			for(iDev=0; iDev<sp->nDevices; iDev++){
				curDev = devStates+iDev;
				cudaSetDevice(cudaContext->devIds[iDev]);
				
				curDev->dt = WLT*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
				curDev->du = WLU*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
				curDev->dv = WLV*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
				
				curDev->cumDt += curDev->dt;
				curDev->cumDu += curDev->du;
				curDev->cumDv += curDev->dv;
				
				prevTuvOff[iDev] = tuvOff[iDev];
				tuvOff[iDev] = curDev->dt + curDev->du*LCELL*WST + curDev->dv*LCELL*WSU*LCELL*WST;
				
				polmove <<< sp->nwg, sp->ws, 0 >>> (sp->nSteps, curDev->seedBuf, curDev->latBuf[curDev->curBuf], curDev->latBuf[curDev->curBuf^0x1], curDev->transBuf, prevTuvOff[iDev], tuvOff[iDev], NWT, NWU, NWV);
			}
			curDev->curBuf ^= 0x1;
		}
		GPULatticeToCPU(sp, ss, devStates, cudaContext);
		RedistribSL(ss, sp);
		GPULoadLattices(sp, ss, devStates, cudaContext);
	}
	
	return 0;
}

