#ifndef _UTIL_H_INCLUDED_
#define _UTIL_H_INCLUDED_
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>

#ifndef UINT_DEF
#define UINT_DEF
typedef unsigned int uint;
#endif

// void PrintPolymer();
void PrintCoor(uint coor);
void PrintUnitVector(uint vec);
int BitCount(uint vec);
double TimeElapsed(struct timeval tStart, struct timeval tEnd);
char* CStringFromFile(char* filename);
char unsigned* BinaryFromFile(char* filename, size_t* n);
int CountBits(uint* ar, int nElem);
long GetTFromName(char* file);
void PrintBonds(uint* bonds, int nBonds);

#endif
