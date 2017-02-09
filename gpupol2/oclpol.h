#include "gpupol.h"


uint GetGpuSite(uint t, uint u, uint v, SimProperties* sp);
uint ValidateAddUnitVectors(uint a, uint b, uint* c);
uint AddUnitVectors(uint a, uint b);
uint IsValid(uint a);
void AddUnitTUV(uint unit, uint* t, uint* u, uint* v);

void SetBondVecs(char* lattice, uint t, uint u, uint v, uint* bonds, uint polSize, uint* labels, SimProperties* sp);
void GetAllRingPolymers();
int GetRingPolymer(SimState* ss, uint t, uint u, uint v, Polymer* pol, char* lattice);
void FindPolymerStart(SimState* ss, uint *t, uint *u, uint *v, Polymer* pol, char* lattice);
void UpdatePolymerModes(SimState* ss);
void UpdatePolymerWithLabels(SimState* ss);


char GetSl(char* lattice, uint t, uint u, uint v, SimProperties* sp);
char GetSlSite(char* lattice, uint site);
void SetSl(char* lattice, uint t, uint u, uint v, int sl, SimProperties* sp);
void SetSlSite(char* lattice, uint site, int sl);
uint GetBond(char* lattice, uint t, uint u, uint v, int sl, SimProperties* sp);
uint GetBondSite(char* lattice, uint site);
void SetBond(char* lattice, uint t, uint u, uint v, int bond, SimProperties* sp);
void SetBondSite(char* lattice, uint site, int bond);
uint GetLabel(uint* labelLat, uint t, uint u, uint v, SimProperties* sp);
uint GetLabelSite(char* lattice, uint site);
void SetLabel(char* lattice, uint t, uint u, uint v, uint label, uint left, SimProperties* sp);
void SetLabelSite(char* lattice, uint site, uint label, uint left);
void CopyGPUToCPULattice(char* gpuLattice, char* cpuLattice, int tOff, int uOff, int vOff, int dt, int du, int dv, SimProperties* sp);
