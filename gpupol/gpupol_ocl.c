#include "gpupol_ocl.h"
#include "gpupol.h"
#include "gpupol_ocl_def.h"

int GPULibInit(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* oclContext){
	
	char binFile[100];
	char kernelName[100];
	sprintf(kernelName, "polmove");
#if POL_TYPE == POL_LIN
	sprintf(binFile, "linpol.gpu");
#elif POL_TYPE == POL_RING
	sprintf(binFile, "ringpol.gpu");
#endif
	
	cl_int ret;
	cl_platform_id platformId;
	cl_uint nPlatforms;
	cl_device_id deviceIds[32];
	cl_device_type deviceType;
	cl_uint nDevices;
	char tString[1000];
	int i;
	sp->nDevices=0;
	unsigned char* fftSrc[32];
	size_t binSizes[32];
	GPUDeviceState* curDev;
	
	clGetPlatformIDs(1, &platformId, &nPlatforms);
	clGetPlatformInfo(platformId, CL_PLATFORM_VENDOR, 1000, (void*) tString, NULL);
	clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL, 32, deviceIds, &nDevices);
	
	for(i=0; i<nDevices; i++){
		curDev = devStates+sp->nDevices;
		ret=clGetDeviceInfo(deviceIds[i], CL_DEVICE_TYPE, sizeof(deviceType), (void*) &deviceType, NULL);
		if(deviceType==CL_DEVICE_TYPE_CPU) {
			printf("Found CPU\n");
		}
		else if(deviceType==CL_DEVICE_TYPE_GPU) {
// 			printf("Found GPU\n");
			curDev->devId = deviceIds[i];
			ret=clGetDeviceInfo(deviceIds[i], CL_DEVICE_VENDOR, 1000, (void*) tString, NULL);
			if(!strcmp(tString, "Advanced Micro Devices, Inc."))
				curDev->devType = __ATI_GPU__;
			else if(!strncmp(tString, "NVIDIA", 6))
				curDev->devType = __NVIDIA_GPU__;
			else{
				printf("Unknown GPU device (%s) found, continuing..\n", tString);
				continue;
			}
			sp->nDevices++;
		}
		else sprintf(tString, "other");
	}
	size_t binSize;
	fftSrc[0] = BinaryFromFile(binFile, &binSize);
	if(!binSize || !fftSrc[0]){
		printf("Error reading binary file %s\n", binFile);
		return 196;
	}
	for(i=0; i<sp->nDevices; i++){
		fftSrc[i] = fftSrc[0];
		binSizes[i] = binSize;
	}
	for(i=0; i<sp->nDevices; i++)
		deviceIds[i] = devStates[i].devId;
	
	oclContext->context = clCreateContext(NULL, sp->nDevices, deviceIds, NULL, NULL, &ret);
	if(ret==CL_SUCCESS);
	else{
		printf("Failed to create context\n");
		return 192;
	}
	
	for(i=0; i<sp->nDevices; i++){
		curDev = devStates+i;
		curDev->cQueue = clCreateCommandQueue(oclContext->context, curDev->devId, 0, &ret);
		if(ret==CL_SUCCESS);
		else{
			printf("Failed to create command queue: error %i\n", ret);
			return 193;	
		}
	}

	oclContext->program=clCreateProgramWithBinary(oclContext->context, sp->nDevices, deviceIds, binSizes, (const unsigned char**) fftSrc, NULL, &ret);
	if(ret==CL_SUCCESS);
	else{
		printf("Failed to create program from binary: ret = %i\n", ret);
		printf("Size of binary = %lu bytes\n", binSize);
		clGetProgramBuildInfo(oclContext->program, devStates[0].devId, CL_PROGRAM_BUILD_LOG, 1000, tString, NULL);
		printf("\n%s\n", tString);
		return 194;
	}

	ret = clBuildProgram(oclContext->program, sp->nDevices, deviceIds, NULL, NULL, NULL);
	if(ret==CL_SUCCESS);
	else{
		printf("Failed to build program\n");
		return 195;	
	}
	
	for(i=0; i<sp->nDevices; i++){
		curDev = devStates+i;
		curDev->kernel = clCreateKernel(oclContext->program, kernelName, &ret);
		assert(ret == CL_SUCCESS);
	}
	if(ret==CL_SUCCESS);
	else
		return 196;	
	free(fftSrc[0]);
	return CL_SUCCESS;
}

int GPULibLoadBuffers(SimProperties* sp, SimState* ss, GPUDeviceState * devStates, GPULibContext* oclContext){
	uint globalWs;
	uint i;
	cl_int ret;
	GPUDeviceState* curDev;
	
	globalWs = sp->nwg*sp->ws;
	CreateGPULattice(ss, sp);
	for(int iDev=0; iDev<sp->nDevices; iDev++){
		curDev = devStates + iDev;
		curDev->seedBuf = clCreateBuffer(oclContext->context, CL_MEM_READ_WRITE, sizeof(cl_uint)*globalWs*2*sp->R, NULL, &ret);
		if(ret != CL_SUCCESS) printf("ret = %i\n", ret); 
		curDev->latBuf  = clCreateBuffer(oclContext->context, CL_MEM_READ_WRITE, sizeof(cl_uint)*sp->gpuLatSize,   NULL, &ret);
		if(ret != CL_SUCCESS) printf("ret = %i (%i)\n", ret, sp->gpuLatSize);
		assert(ret==CL_SUCCESS);
		curDev->transBuf = clCreateBuffer(oclContext->context, CL_MEM_READ_ONLY, sizeof(cl_uint)*4, NULL, &ret);
		assert(ret==CL_SUCCESS);
		
		
		ret = clEnqueueWriteBuffer(curDev->cQueue, curDev->seedBuf, CL_FALSE, 0, sizeof(cl_uint)*globalWs*2*sp->R, (void*) (ss[iDev].seeds), 0, NULL, NULL);
		if(ret != CL_SUCCESS) printf("ret = %i\n", ret);
		assert(ret == CL_SUCCESS);
		ret = clEnqueueWriteBuffer(curDev->cQueue, curDev->latBuf, CL_FALSE, 0, sizeof(cl_uint)*sp->gpuLatSize, (void*) ss[iDev].gpuLattice, 0, NULL, NULL);
		assert(ret == CL_SUCCESS);
		clEnqueueWriteBuffer(curDev->cQueue, curDev->transBuf, CL_FALSE, 0, sizeof(cl_uint)*4, (void*) (sp->trans), 0, NULL, NULL);
		assert(ret==CL_SUCCESS);
	}
	for(i=0; i<sp->nDevices; i++)
		clFlush(devStates[i].cQueue);
	for(i=0; i<sp->nDevices; i++)
		clFinish(devStates[i].cQueue);
	
	for(i=0; i<sp->nDevices; i++){
		curDev = devStates+i;
		assert( clSetKernelArg(curDev->kernel, 0, sizeof(cl_uint), (void*)&sp->nSteps)       == CL_SUCCESS);
		assert( clSetKernelArg(curDev->kernel, 1, sizeof(void*),   (void*)&curDev->seedBuf)  == CL_SUCCESS);
		assert( clSetKernelArg(curDev->kernel, 2, sizeof(void*),   (void*)&curDev->latBuf)   == CL_SUCCESS);
		assert( clSetKernelArg(curDev->kernel, 3, sizeof(void*),   (void*)&curDev->transBuf) == CL_SUCCESS);
	}
	
	return 0;
}




int GPULibRelease(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* oclContext){
	int i;
	for(i=0; i<sp->nDevices; i++){
		clReleaseMemObject(devStates[i].seedBuf);
		clReleaseMemObject(devStates[i].latBuf);
		clReleaseKernel(devStates[i].kernel);
	}
	clReleaseProgram(oclContext->program);
	for(i=0; i<sp->nDevices; i++)
		clReleaseCommandQueue(devStates[i].cQueue);
	clReleaseContext(oclContext->context); 
	return 0;
}


int GPULibRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* oclContext, int nTime){
	int i,j,iDev;
	size_t localWs, globalWs;
	cl_int ret;
	cl_uint tOff, uOff, vOff, tuvOff;
	cl_uint nwt, nwu, nwv;
	localWs = sp->ws;
	globalWs = localWs*sp->nwg;
	cl_event **kerEvent;
	GPUDeviceState* curDev;
	
	nwt = sp->nwt; nwu = sp->nwu; nwv = sp->nwv;
	kerEvent = (cl_event**) malloc(sizeof(cl_event*)*sp->nDevices);
	for(iDev=0; iDev<sp->nDevices; iDev++){
		kerEvent[iDev] = (cl_event*) malloc(sizeof(cl_event)*nTime);
		assert( clSetKernelArg(devStates[iDev].kernel, 5, sizeof(cl_uint), (void*)&nwt) );
		assert( clSetKernelArg(devStates[iDev].kernel, 6, sizeof(cl_uint), (void*)&nwu) );
		assert( clSetKernelArg(devStates[iDev].kernel, 7, sizeof(cl_uint), (void*)&nwv) );
		
	}
	j=0;
	for(i=0; i<nTime; i++){
		for(iDev=0; iDev<sp->nDevices; iDev++){
			curDev = devStates+iDev;
			tOff = (LCELL*WST)*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			uOff = (LCELL*WSU)*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			vOff = (LCELL*WSV)*(double)RNG_FAC*Rng4(&ss[iDev].rngState);
			
			tuvOff = tOff+uOff*LCELL*WST+vOff*LCELL*LCELL*WST*WSU;
			
			assert( clSetKernelArg(curDev->kernel, 4, sizeof(cl_uint), (void*)&tuvOff) == CL_SUCCESS);
			if(i==0)
				ret = clEnqueueNDRangeKernel(curDev->cQueue, curDev->kernel, 1, NULL, &globalWs, &localWs, 0, NULL, kerEvent[iDev]+j);
			else
				ret = clEnqueueNDRangeKernel(curDev->cQueue, curDev->kernel, 1, NULL, &globalWs, &localWs, 1, kerEvent[iDev]+j-1, kerEvent[iDev]+j);
			if(ret != CL_SUCCESS) printf("Failed enqueing polmove kernel (ret=%i)\n", ret);
			assert( ret == CL_SUCCESS );
			
		}
		j++;
	}
	for(iDev=0; iDev<sp->nDevices; iDev++){
		ret = clFinish(devStates[iDev].cQueue);
		if(ret != CL_SUCCESS){
			printf("Failed execution!\n");
			printf("ret = %i\n", ret);
			assert(ret == CL_SUCCESS);
		}
	}
	for(iDev=0; iDev<sp->nDevices; iDev++){
		for(j=0; j<nTime; j++){
			clReleaseEvent(kerEvent[iDev][j]);
		}
		free(kerEvent[iDev]);
	}
	
// 	memset(sp->gpuLattice, 0, sizeof(cl_uint)*sp->gpuLatSize);
	for(iDev=0; iDev<sp->nDevices; iDev++){
		curDev = devStates+iDev;
		ret = clEnqueueReadBuffer(curDev->cQueue, curDev->latBuf, CL_FALSE, 0, sizeof(cl_uint)*sp->gpuLatSize, (void*)ss[iDev].gpuLattice, 0, NULL, NULL);
		if(ret != CL_SUCCESS){
			printf("Failed to read buffer to CPU host (ret = %i)\n", ret);
			assert(ret == CL_SUCCESS);
		}
	}

	for(iDev=0; iDev<sp->nDevices; iDev++){
		ret = clFinish(devStates[iDev].cQueue);
		if(ret != CL_SUCCESS){
			printf("Error reading data back (ret = %i)\n", ret);
			return 192;
		}
	}
	free(kerEvent);
	CopyGPUToCPULattice();
	GetAllRingPolymers();
	for(iDev=0; iDev<sp->nDevices; iDev++)
		UpdatePolymerWithLabels(ss+iDev);
	
	return 0;
}
