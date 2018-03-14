.PHONY: clean install

EXECS= ee denspol denspol_hp contact_map denspol_scaleup
BIN_DIR=../bin
VPATH=../util
CFLAGS= -Wall -O2 
# CFLAGS=-g -Wall
LDFLAGS= -lm
EXTRA_INC=../util
RING_OBJ=denspol.o denspol_init.o denspol_lib.o denspol_io.o timer.o denspol_dbg.o file.o
$(EXECS): override CFLAGS+= -I$(EXTRA_INC) -std=c99 -DRELEASE=$(RELEASE)

all: $(EXECS)

%_hp.o: %.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__HP_ENABLED__


ee: ee.o lattice.o


denspol: $(RING_OBJ) denspol_mc.o denspol_mc_lin.o denspol_mc_ring.o

denspol_hp: $(RING_OBJ) denspol_mc_hp.o denspol_mc_lin_hp.o denspol_mc_ring_hp.o
	$(CC) $(LDFLAGS) $^ -o $@ 

contact_map: contact_map.o denspol_lib.o

denspol_scaleup: denspol_scaleup.o denspol_lib.o

clean:
	rm -f $(EXECS) *.o *~ 

install:       
	cp $(EXECS) $(BIN_DIR)
