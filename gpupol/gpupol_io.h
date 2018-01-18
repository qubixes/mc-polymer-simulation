#include "gpupol.h"

void PrintPol(Polymer* pol);
void PrintBinary(uint num);
void PrintSite(uint* gLattice, uint t, uint u, uint v);
void PrintSummary(SimProperties* sp, SimState* ss);
int ReadLatticeFile(SimState* ss, SimProperties* sp, char* file);
int WriteLatticeFile(SimProperties* sp, SimState* ss, char* file);
void WriteRingPolymer(FILE* pFile, Polymer* pol);
void WriteSimulationSettings(SimProperties* sp, SimState* ss);
