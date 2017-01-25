#ifndef __TIMER_INCLUDED__
#define __TIMER_INCLUDED__
#include <stdlib.h>
#include <sys/time.h>
#include <stdio.h>

typedef struct Timer{
	struct timeval start, end;
	double elapsed;
	int running;
}Timer;

void TimerStart(Timer* t);
void TimerStop(Timer* t);
double TimerElapsed(Timer* t);
void TimerContinue(Timer* t);
#endif