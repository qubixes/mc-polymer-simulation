#include "gpupol_def.h"
#include "gpupol_cuda_def.h"



int GPULibInit(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* oclContext);
int GPULibLoadBuffers(SimProperties* sp, SimState* ss, GPUDeviceState * devStates, GPULibContext* oclContext);
int GPULibRelease(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* oclContext);
int GPULibRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* oclContext, int time);
int GPULoadLattices(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext);
int GPULatticeToCPU(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext);
int GPULibEqRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* cudaContext, int nTime);