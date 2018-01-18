#ifndef GPUPOL_OCL_DEF_H_INC
#define GPUPOL_OCL_DEF_H_INC

#if __APPLE__
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif



typedef struct GPULibContext{
	cl_context context;
	cl_program program;
}GPULibContext;

typedef struct GPUDeviceState{
	cl_command_queue cQueue;
	cl_device_id devId;
	uint devType;
	cl_kernel kernel;
	
	cl_mem latBuf, seedBuf, transBuf;
}GPUDeviceState;
#endif
