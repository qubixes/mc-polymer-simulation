#include "gpupol.h"


uint GetGpuSite(uint t, uint u, uint v, SimProperties* sp);
void CopyCPUToGPULattice();
void CopyGPUToCPULattice();
uint ValidateAddUnitVectors(uint a, uint b, uint* c);
uint AddUnitVectors(uint a, uint b);
uint IsValid(uint a);
uint GetSl(uint* sl, uint t, uint u, uint v);
uint GetLabel(uint* labelLat, uint t, uint u, uint v);
void AddUnitTUV(uint unit, uint* t, uint* u, uint* v);
void SetBondVecs(SimState* ss, uint t, uint u, uint v, uint* bonds, uint polSize, uint* labels);
void SetPrevNextVecs(uint* gLattice, uint* gSl, uint* t, uint* u, uint* v, uint unit, uint sl);
void GetSiteSlot(uint t, uint u, uint v, uint* site, uint* slot, SimProperties* sp, uint mSpin);
void GetAllRingPolymers();
int GetRingPolymer(SimState* ss, uint t, uint u, uint v, Polymer* pol, uint* lattice, uint* slLattice, uint* labLattice);
void FindPolymerStart(SimState* ss, uint *t, uint *u, uint *v, Polymer* pol, uint* lattice, uint* slLattice, uint* labLattice);
Coor GetCoor(uint t, uint u, uint v);
Coor DGetCoor(double t, double u, double v);
void SetSite(uint* gLattice, uint t, uint u, uint v, uint val);
uint GetNextSite(uint* gLattice, uint* t, uint* u, uint* v, uint dirMode);
int RSQA(int* d);
int PolDistance(uint* start, uint* end, uint nMono);
int MinDistance(uint* start, uint* end);
double DRSQ(Coor* a, Coor* b);
double Distance(Coor* a, Coor* b);
Coor AddUnitToCoor(uint unit, Coor coor);
Coor CoorToBaseBlock(Coor coor);
Coor GetCMS(Coor* coor, uint polSize);
int CEQ(Coor a, Coor b);
void SetSecAttPolymer(uint t, uint u, uint v, Polymer* pol);
void UpdatePolymerModes(SimState* ss);
void UpdatePolymerWithLabels(SimState* ss);
