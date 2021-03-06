#ifndef GPUPOL_CUDA_DEF_H_INC
#define GPUPOL_CUDA_DEF_H_INC

#include <cuda.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct GPULibContext{
	int devIds[8];
}GPULibContext;

typedef struct GPUDeviceState{
	int devId;
	char **latBuf;
	unsigned int *transBuf; 
	uint4 *seedBuf;
	int curBuf;
	long cumDt, cumDu, cumDv;
	int dt, du, dv;
}GPUDeviceState;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
	if (code != cudaSuccess){
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}
#endif
