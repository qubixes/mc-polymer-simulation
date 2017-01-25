.PHONY: clean install

LIN_EXECS=efpol_lin 
RING_EXECS=efpol_ring  
EXECS=$(LIN_EXECS) $(RING_EXECS)  
BIN_DIR=../bin
VPATH=../util
CFLAGS= -Wall -O2 -std=c99
# CFLAGS=-g -std=c99 -Wall
LDFLAGS=-lpthread -lm
EXTRA_INC=../util
$(EXECS): override CFLAGS+= -I$(EXTRA_INC)
$(RING_EXECS): EXTRA_CFLAGS= -DPOL_TYPE=POL_RING
$(LIN_EXECS): EXTRA_CFLAGS= -DPOL_TYPE=POL_LIN

%_ring.o: %.c
	$(CC) -c  $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS)

%_lin.o: %.c
	$(CC) -c  $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS)

all: $(EXECS)

cpupol: cpupol.o timer.o
bitpol: bitpol.o timer.o
bitpol2_lin: bitpol2_lin.o timer.o
bitpol2_ring: bitpol2_ring.o timer.o
efpol_lin: efpol_lin.o efpol_init_lin.o efpol_io_lin.o timer.o efpol_lib_lin.o efpol_meas_lin.o
efpol_ring: efpol_ring.o efpol_init_ring.o efpol_io_ring.o timer.o efpol_lib_ring.o efpol_meas_ring.o
sspol_lin: sspol_lin.o efpol_init_lin.o efpol_io_lin.o timer.o efpol_lib_lin.o
simpol: simpol.o timer.o
mtpol: mtpol.o timer.o sync.o
sem_test: sem_test.o timer.o sync.o

clean:
	rm -f $(EXECS) *.o *~ 

install:       
	cp $(EXECS) $(BIN_DIR)