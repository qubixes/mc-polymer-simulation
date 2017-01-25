#include "timer.h"

void TimerStart(Timer* t){
	t->elapsed = 0;
	t->running = 1;
	gettimeofday(&(t->start), NULL);
}

void TimerStop(Timer* t){
	t->elapsed += TimerElapsed(t);
	t->running = 0;
}

double TimerElapsed(Timer* t){
	if(t->running)
		gettimeofday(&(t->end), NULL);
	return ((t->end.tv_sec-t->start.tv_sec)+1e-6*(t->end.tv_usec-t->start.tv_usec));
}

void TimerContinue(Timer* t){
	if(t->running) return;
	t->running = 1;
	gettimeofday(&(t->start), NULL);
}
