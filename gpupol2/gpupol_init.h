#include "gpupol.h"
#include "util.h"
#include "oclpol.h"

void LoadPolymers(SimProperties* sp, SimState* ss);
int SimStateInit(SimState* ss, SimProperties* sp);
int SimPropDefaultInit(SimProperties* sp);
void GetBestPolymerPartitioning(SimProperties* sp, SimState* ss, int polLength, double density);
void ConstructCubeRingPolymer(SimState* ss, SimProperties* sp, uint tStart, uint uStart, uint vStart, uint dt, uint du, uint dv, uint sl);
int ConstructRingPolymer(SimState* ss, SimProperties* sp, uint t, uint u, uint v);
void AddPolymerToSim(SimState* ss, SimProperties* sp, uint t, uint u, uint v, int length, uint* bonds);
int SetGPUTrans(uint* trans);
void LoadPolymers(SimProperties* sp, SimState* ss);
void ClearLattice(SimState* ss, SimProperties* sp);
void RedistribSL(SimState* ss, SimProperties* sp);
void ConstructBarRingPolymer(SimState* ss, SimProperties* sp, int polId, int length);
void SetBoxDimension(SimProperties* sp, int LT, int LU, int LV);
