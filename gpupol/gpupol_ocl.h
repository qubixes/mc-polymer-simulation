#include "gpupol_def.h"
#include "gpupol_ocl_def.h"

int GPULibInit(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* oclContext);
int GPULibLoadBuffers(SimProperties* sp, SimState* ss, GPUDeviceState * devStates, GPULibContext* oclContext);
int GPULibRelease(SimProperties* sp, GPUDeviceState* devStates, GPULibContext* oclContext);
int GPULibRun(SimProperties* sp, SimState* ss, GPUDeviceState* devStates, GPULibContext* oclContext, int time);
