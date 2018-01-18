#include "gpupol_cuda.h"
#include "gpupol.h"
#include "gpupol_cuda_def.h"
#include "ringpol_cu.h"

int GPULibInit(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* cudaContext){
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
	return 0;
}

int GPULibLoadBuffers(SimProperties* sp, SimState* ss, GPUDeviceState * devStates, GPULibContext* cudaContext){
	uint globalWs;
	GPUDeviceState* curDev;
	
	globalWs = sp->nwg*sp->ws;
	CreateGPULattice(ss, sp);
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		curDev = devStates + iDev;
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		
		gpuErrchk( cudaMalloc(&(curDev->seedBuf),  sizeof(uint)*globalWs*2*sp->R) );
		gpuErrchk( cudaMalloc(&(curDev->latBuf),   sizeof(uint)*sp->gpuLatSize) );
		gpuErrchk( cudaMalloc(&(curDev->transBuf), sizeof(uint)*4) );
		
		gpuErrchk( cudaMemcpy(curDev->seedBuf,  ss[iDev].seeds,      sizeof(uint)*globalWs*2*sp->R, cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(curDev->latBuf,   ss[iDev].gpuLattice, sizeof(uint)*sp->gpuLatSize,   cudaMemcpyHostToDevice) );
		gpuErrchk( cudaMemcpy(curDev->transBuf, sp->trans,           sizeof(uint)*4,                cudaMemcpyHostToDevice) );
	}
	return 0;
}

int GPULibRelease(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* cudaContext){
	GPUDeviceState* curDev;
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		curDev = devStates+iDev;
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaFree(curDev->seedBuf) );
		gpuErrchk( cudaFree(curDev->latBuf) );
		gpuErrchk( cudaFree(curDev->transBuf) );
	}
	return 0;
}

int GPULatticeToCPU(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext){
	
	gpuErrchk( cudaPeekAtLastError() );
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		GPUDeviceState* curDev = devStates+iDev;
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaMemcpy(ss[iDev].gpuLattice, curDev->latBuf, sizeof(uint)*sp->gpuLatSize, cudaMemcpyDeviceToHost) );
	}
	
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		gpuErrchk( cudaSetDevice(cudaContext->devIds[iDev]) );
		gpuErrchk( cudaDeviceSynchronize() );
	}
	
	CopyGPUToCPULattice();
	GetAllRingPolymers();
	for(int iDev=0; iDev<sp->nDevices; iDev++)
		UpdatePolymerWithLabels(ss+iDev);
	return 0;	
}

int GPULibRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext, int nTime){
	int iDev;
	int* tuvOff;
	uint NWT, NWU, NWV;
	int dt, du, dv;
	GPUDeviceState* curDev;
	
	tuvOff = (int*) malloc(sizeof(int)*sp->nDevices);
	NWT=sp->nwt; NWU=sp->nwu; NWV=sp->nwv;
	
	for(int i=0; i<nTime; i++){
		for(iDev=0; iDev<sp->nDevices; iDev++){
			curDev = devStates+iDev;
			cudaSetDevice(cudaContext->devIds[iDev]);
			
			dt = WLT*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			du = WLU*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			dv = WLV*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			
			tuvOff[iDev] = dt + du*LCELL*WST + dv*LCELL*WSU*LCELL*WST;
			
			polmove <<< sp->nwg, sp->ws, 0 >>> (sp->nSteps, curDev->seedBuf, curDev->latBuf, curDev->transBuf, tuvOff[iDev], NWT, NWU, NWV);
		}
	}
	
	GPULatticeToCPU(sp, ss, devStates, cudaContext);
	
	return 0;
}




/*
int GPULibRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext, int nTime){
	int iDev;
	uint tOff, uOff, vOff, tuvOff;
	uint NWT, NWU, NWV;
	GPUDeviceState* curDev;
	
	NWT=sp->nwt; NWU=sp->nwu; NWV=sp->nwv;
	
	for(int i=0; i<nTime; i++){
		for(iDev=0; iDev<sp->nDevices; iDev++){
			curDev = devStates+iDev;
			cudaSetDevice(cudaContext->devIds[iDev]);
			
			tOff = (LCELL*WST)*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			uOff = (LCELL*WSU)*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			vOff = (LCELL*WSV)*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			
			tuvOff = tOff+uOff*LCELL*WST+vOff*LCELL*LCELL*WST*WSU;
// 			printf("nStep=%i\n", sp->nSteps);
			polmove <<< sp->nwg, sp->ws, 0 >>> (sp->nSteps, curDev->seedBuf, curDev->latBuf, curDev->transBuf, tuvOff, NWT, NWU, NWV);
		}
	}
	GPULatticeToCPU(sp, ss, devStates, cudaContext);
	return 0;
}*/
