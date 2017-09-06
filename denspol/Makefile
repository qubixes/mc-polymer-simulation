.PHONY: clean install

EXECS= ee denspol_ring denspol_lin
BIN_DIR=../bin
VPATH=../util
CFLAGS= -Wall -O2 
# CFLAGS=-g -std=c99 -Wall
LDFLAGS= -lm
EXTRA_INC=../util
$(EXECS): override CFLAGS+= -I$(EXTRA_INC) -std=c99

%_ring.o: %.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -DPOL_TYPE=POL_TYPE_RING

%_lin.o: %.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -DPOL_TYPE=POL_TYPE_LIN
	



all: $(EXECS)

ee: ee.o lattice.o

# denspol: denspol.o denspol_init.o denspol_mc.o  denspol_lib.o denspol_io.o timer.o denspol_dbg.o file.o

denspol_ring: denspol_ring.o denspol_init_ring.o denspol_lib_ring.o denspol_io_ring.o timer.o denspol_dbg_ring.o denspol_mc_ring.o file.o

denspol_lin: denspol_lin.o denspol_init_lin.o denspol_lib_lin.o denspol_io_lin.o denspol_dbg_lin.o denspol_mc_lin.o timer.o file.o

clean:
	rm -f $(EXECS) *.o *~ 

install:       
	cp $(EXECS) $(BIN_DIR)