.PHONY: clean all install

CFLAGS=-O2 -Wall -std=c99
NVFLAGS=-O2 --compiler-options "-DCUDA_GPULIB -fPIC" -x cu -DRELEASE=$(RELEASE)
LDFLAGS=-lm
CC=gcc
CXX=g++
NVCC=nvcc -m64
BIN_DIR=../bin

EXEC=label_test
OCL_EXECS= oclcc gpupol_ocl_ring gpupol_ocl_lin
CUDA_EXECS= gpupol_cuda_ring gpupol_cuda_lin
RING_EXECS= gpupol_ocl_ring gpupol_cuda_ring
LIN_EXECS= gpupol_ocl_lin gpupol_cuda_lin

OCL_BINS=linpol.gpu ringpol.gpu
DEFINES=
CUDA_PATH=/usr/local/cuda/lib
OCL_INC=-I$(CUDA_PATH)/include
NV_OBJ=ringpol.o

TARGETS=$(EXECS) $(OCL_EXECS) $(OCL_BINS) $(CUDA_EXECS)

$(OCL_EXECS): override CFLAGS+= $(OCL_INC)
$(OCL_EXECS): override CFLAGS+= -DOPENCL_GPULIB
$(RING_EXECS): EXTRA_CFLAGS= -DPOL_TYPE=POL_RING
$(LIN_EXECS): EXTRA_CFLAGS= -DPOL_TYPE=POL_LIN
OS=$(shell uname)
$(CUDA_EXECS): override LDFLAGS+= -lcudart
ifeq ($(OS), Darwin)
 $(OCL_EXECS): override LDFLAGS+= -framework OpenCL
 $(CUDA_EXECS): override CXX=/usr/bin/clang
 $(CUDA_EXECS): override LDFLAGS+= -L/usr/local/cuda/lib
 override NVCC+= -ccbin=/usr/bin/clang
else
 override TARGETS=$(CUDA_EXECS)
 $(OCL_EXECS): override LDFLAGS+= -lOpenCL
 $(CUDA_EXECS): override LDFLAGS+= -L$(CUDA_HOME)/lib64
 override NVCC+= -ccbin clang-3.8
endif



%.gpu: %.ocl
	./oclcc gpu $< $@

%.cpu: %.ocl
	./oclcc cpu $< $@

%_ring.o: %.c
	$(CC) -c  $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS)

%_lin.o: %.c
	$(CC) -c  $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS)

%_cu_ring.o: %.c
	$(NVCC) -c  $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)

%_cu_lin.o: %.c
	$(NVCC) -c  $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)

%_cu_ring.o: %.cu
	$(NVCC) -c $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)

%_cu_lin.o: %.cu
	$(NVCC) -c $^ -o $@ $(NVFLAGS) $(EXTRA_CFLAGS)


all: $(TARGETS)
$(OCL_BINS): oclcc


oclcc: oclcc.o util.o
oclcc.o: oclcc.c util.h

gpupol_ocl_ring: gpupol_ring.o oclpol_ring.o gpupol_init_ring.o gpupol_ocl_ring.o gpupol_io_ring.o util_ring.o
	$(CC) $^ -o $@ $(LDFLAGS)
	
gpupol_ocl_lin: gpupol_lin.o oclpol_lin.o gpupol_init_lin.o gpupol_ocl_lin.o gpupol_io_lin.o util_lin.o
	$(CC) $^ -o $@ $(LDFLAGS)

gpupol_cuda_ring: gpupol_cu_ring.o oclpol_cu_ring.o gpupol_init_cu_ring.o gpupol_cuda_cu_ring.o gpupol_io_cu_ring.o util_cu_ring.o ringpol_cu_ring.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gpupol_cuda_lin: gpupol_cu_lin.o oclpol_cu_lin.o gpupol_init_cu_lin.o gpupol_cuda_cu_lin.o gpupol_io_cu_lin.o util_cu_lin.o linpol_cu_lin.o
	$(CXX) $^ -o $@ $(LDFLAGS)

gpupol.o: gpupol.c oclpol.h gpupol.h
oclpol.o: oclpol.c oclpol.h gpupol.h

install:
	cp $(TARGETS) $(BIN_DIR)

clean: 
	rm -f $(EXECS) *.o *~  $(OCL_BINS) $(OCL_EXECS) $(CUDA_EXECS) 
