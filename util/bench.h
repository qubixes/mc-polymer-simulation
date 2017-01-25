#include "timer.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define L 96
#define LAT_SIZE (L*L*L)

void CleanDir(char* dir);
double ComputeAverageDif(double** timings, int nRemove, int nMeas);
