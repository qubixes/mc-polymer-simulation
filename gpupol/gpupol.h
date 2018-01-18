#ifndef _GPUPOL_H_INCLUDED_
#define _GPUPOL_H_INCLUDED_

#define TUD 2
#define LCELL 4
#define WST 8
#define WSU 8
#define WSV 8
#define WS (WST*WSU*WSV)

#define WLT (LCELL*WST)
#define WLU (LCELL*WSU)
#define WLV (LCELL*WSV)

#define MSPIN 4
#define SL_MSPIN 32
#define LAB_MSPIN 16
#define CELL_SIZE (LCELL*LCELL*LCELL)

#define PREV_MASK 0x00555555
#define NEXT_MASK 0x00aaaaaa

#define PREV1_MASK 0x00041041
#define PREV2_MASK (PREV1_MASK<<2)
#define PREV3_MASK (PREV2_MASK<<2)

#define NEXT1_MASK (PREV1_MASK<<1)
#define NEXT2_MASK (NEXT1_MASK<<2)
#define NEXT3_MASK (NEXT2_MASK<<2)

#define SL1_MASK 0x01000000
#define SL2_MASK (SL1_MASK<<2)
#define SL3_MASK (SL2_MASK<<2)
#define SL_MASK (SL1_MASK|SL2_MASK|SL3_MASK)

#define ED1_MASK 0x02000000
#define ED2_MASK (ED1_MASK<<2)
#define ED3_MASK (ED2_MASK<<2)
#define ED_MASK (ED1_MASK|ED2_MASK|ED3_MASK)

#define SUB1_MASK (PREV1_MASK|NEXT1_MASK|SL1_MASK|ED1_MASK)
#define SUB2_MASK (PREV2_MASK|NEXT2_MASK|SL2_MASK|ED2_MASK)
#define SUB3_MASK (PREV3_MASK|NEXT3_MASK|SL3_MASK|ED3_MASK)

#define ENT_PIXEL_ENABLE (SUB1_MASK|SUB2_MASK|SUB3_MASK)
#define TOT_PIX ENT_PIXEL_ENABLE

#define SUBV1_MASK (PREV1_MASK|NEXT1_MASK|SL1_MASK)
#define SUBV2_MASK (PREV2_MASK|NEXT2_MASK|SL2_MASK)
#define SUBV3_MASK (PREV3_MASK|NEXT3_MASK|SL3_MASK)
#define SUBV_MASK (SUBV1_MASK|SUBV2_MASK|SUBV3_MASK)

#define RNG_FAC_W1 2.328306437080797375431469961868E-10
#define RNG_FAC    2.3283064365386962890625E-10

#define __NVIDIA_GPU__ 1
#define __ATI_GPU__ 2
#define __CPU__ 4

#define TRUE 1
#define FALSE 0


#define POL_LIN 1
#define POL_RING 2
// #define POL_TYPE POL_RING

#define N_MAX_GPU 16

#define PRINT_MONO 0x1
#define PRINT_LENGTH 0x2
#define DIR_FORWARD 0x4
#define DIR_BACKWARD 0x8
#define PRINTED_BIT 0x80000000


#define TOSTR2(x) #x
#define TOSTR(x) TOSTR2(x)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "gpupol_def.h"
#ifdef OPENCL_GPULIB
	#include "gpupol_ocl_def.h"
	#include "gpupol_ocl.h"
#endif

#ifdef CUDA_GPULIB
	#include "gpupol_cuda_def.h"
	#include "gpupol_cuda.h"
#endif

#include "ocl_util.h"
#include "gpupol_init.h"
#include "gpupol_io.h"


extern SimProperties sp;
extern GPULibContext gpuContext;
extern SimState ss[N_MAX_GPU];
extern GPUDeviceState devices[N_MAX_GPU];

uint DistribOneBits(uint sum);
uint AddUnitVecs(uint a, uint b, uint* c);
uint Rng4(uint4_t* state);


#ifdef IS_MAIN


inline uint DistribOneBits(uint sum){
	sum = sum | (sum >> 6 );
	sum = sum | (sum >> 12);
	sum = sum | (sum << 6 );
	sum = sum | (sum << 12);
	sum = sum | ((sum&0x15)<<24);
	sum = sum | (sum<<1);
	return sum;
}

uint Rng4(uint4_t* state){
	uint b;
	b = ((((*state).x << 6)^(*state).x)>>13);
	(*state).x = ((((*state).x & 4294967294) <<18) ^b);
	b = ((((*state).y << 2)^(*state).y) >> 27);
	(*state).y = ((((*state).y & 4294967288) << 2) ^b);
	b = ((((*state).z << 13)^(*state).z) >>21);
	(*state).z = ((((*state).z & 4294967280) << 7) ^b);
	b = ((((*state).w << 3)^(*state).w) >> 12);
	(*state).w = ((((*state).w & 4294967168) << 13)^b);
	return ((*state).x^(*state).y^(*state).z^(*state).w);
}

#endif
#endif