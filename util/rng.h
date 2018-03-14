#ifndef __RNG_H_INCLUDED__
#define __RNG_H_INCLUDED__
#ifndef IS_MAIN
unsigned int Rng(unsigned int* state);
void Seed(unsigned int* state, unsigned int seed);
double DRng(unsigned int* state);
#else
/** Tausworthe generator 
 **/
unsigned int Rng(unsigned int* state){
  unsigned int b;
  b = (((state[0] << 6 )^state[0]) >> 13);
  state[0] = (((state[0] & 4294967294) << 18) ^b);
  b = (((state[1] << 2 )^state[1]) >> 27);
  state[1] = (((state[1] & 4294967288) << 2 ) ^b);
  b = (((state[2] << 13)^state[2]) >> 21);
  state[2] = (((state[2] & 4294967280) << 7 ) ^b);
  b = (((state[3] << 3 )^state[3]) >> 12);
  state[3] = (((state[3] & 4294967168) << 13) ^b);
  return (state[0]^state[1]^state[2]^state[3]);
}

#define RNG_FAC_W1 2.328306437080797375431469961868E-10
#define RNG_FAC    2.3283064365386962890625E-10

void Seed(unsigned int* state, unsigned int seed){
	state[0] = seed;
	state[1] = seed^0xa741f123;
	state[2] = seed^0x213910ab;
	state[3] = seed^0xcd123e14;
}

double DRng(unsigned int* state){
	return (double)(RNG_FAC*Rng(state));
}
#endif
#endif
