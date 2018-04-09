.PHONY: clean all install

CFLAGS=-O2 -Wall -std=c99 -DRELEASE=$(RELEASE)
NVFLAGS=-O2 --compiler-options "-DCUDA_GPULIB -fPIC" -x cu -DRELEASE=$(RELEASE)
#-lineinfo
LDFLAGS=-lm
CC=gcc
CXX=g++
NVCC=nvcc -m64
BIN_DIR=../bin

EXEC=label_test
CUDA_EXECS= gpupol2_cuda_ring gpupol3_cuda_ring
RING_EXECS= gpupol2_cuda_ring gpupol3_cuda_ring

DEFINES=
CUDA_PATH=/usr/local/cuda/lib
OCL_INC=-I$(CUDA_PATH)/include
NV_OBJ=ringpol.o

TARGETS=$(EXECS) $(OCL_EXECS) $(OCL_BINS) $(CUDA_EXECS)

$(RING_EXECS): EXTRA_CFLAGS= -DPOL_TYPE=POL_RING
OS=$(shell uname)
$(CUDA_EXECS): override LDFLAGS+= -lcudart
gpupol3_cuda_ring: override EXTRA_CFLAGS+= -DLCELL=4
gpupol2_cuda_ring: override EXTRA_CFLAGS+= -DLCELL=3

ifeq ($(OS), Darwin)
 $(CUDA_EXECS): override CXX=/usr/bin/clang
 $(CUDA_EXECS): override LDFLAGS+= -L/usr/local/cuda/lib
 override NVCC+= -ccbin=/usr/bin/clang
else
 override TARGETS=$(CUDA_EXECS)
 $(CUDA_EXECS): override LDFLAGS+= -L$(CUDA_HOME)/lib64
 NVCC += -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35 -gencode=arch=compute_35,code=compute_35
 NVCC += -ccbin clang-3.8
endif



# %.gpu: %.ocl
# 	./oclcc gpu $< $@

# %.cpu: %.ocl
# 	./oclcc cpu $< $@

# %_ring.o: %.c
# 	$(CC) -c  $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS)

# %_lin.o: %.c
# 	$(CC) -c  $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS)

%3_cu_ring.o: %.c
	$(NVCC) -c  $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)
%4_cu_ring.o: %.c
	$(NVCC) -c  $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)

# %_cu_lin.o: %.c
# 	$(NVCC) -c  $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)

%3_cu_ring.o: %.cu
	$(NVCC) -c $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)
%4_cu_ring.o: %.cu
	$(NVCC) -c $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)

# %_cu_lin.o: %.cu
# 	$(NVCC) -c $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)


all: $(TARGETS)



gpupol3_cuda_ring: gpupol24_cu_ring.o oclpol4_cu_ring.o gpupol_init4_cu_ring.o gpupol_cuda4_cu_ring.o gpupol_io4_cu_ring.o util4_cu_ring.o ringpol34_cu_ring.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gpupol2_cuda_ring: gpupol23_cu_ring.o oclpol3_cu_ring.o gpupol_init3_cu_ring.o gpupol_cuda3_cu_ring.o gpupol_io3_cu_ring.o util3_cu_ring.o ringpol23_cu_ring.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gpupol.o: gpupol.c oclpol.h gpupol.h
oclpol.o: oclpol.c oclpol.h gpupol.h

install:
	cp $(TARGETS) $(BIN_DIR)

clean: 
	rm -f $(EXECS) *.o *~  $(OCL_BINS) $(OCL_EXECS) $(CUDA_EXECS) 
