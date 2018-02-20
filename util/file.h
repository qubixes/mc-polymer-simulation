#ifndef __FILE_H_INC__
#define __FILE_H_INC__
#define _XOPEN_SOURCE
#include <sys/stat.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

int NeedsUpdate(char* src, char* dst);
int NeedsUpdatePath(char* src, char* dst, char* path);
int GetLastT(char* dir);
#endif