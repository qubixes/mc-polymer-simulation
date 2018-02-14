.PHONY: clean install

EXECS= ee denspol_ring denspol_lin contact_map denspol_ring_hp denspol_lin_hp denspol_ring_topo denspol_lin_topo denspol_ring_hp_topo denspol_lin_hp_topo
BIN_DIR=../bin
VPATH=../util
CFLAGS= -Wall -O2 
# CFLAGS=-g -std=c99 -Wall
LDFLAGS= -lm
EXTRA_INC=../util
RING_OBJ=denspol_ring.o denspol_init_ring.o denspol_lib_ring.o denspol_io_ring.o timer.o denspol_dbg_ring.o file.o
LIN_OBJ=denspol_lin.o denspol_init_lin.o denspol_lib_lin.o denspol_io_lin.o denspol_dbg_lin.o timer.o file.o
$(EXECS): override CFLAGS+= -I$(EXTRA_INC) -std=c99 -DRELEASE=$(RELEASE)

all: $(EXECS)

%_ring.o: %.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -DPOL_TYPE=POL_TYPE_RING

%_lin.o: %.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -DPOL_TYPE=POL_TYPE_LIN

denspol_mc_ring_hp.o: denspol_mc_ring.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__HP_ENABLED__

denspol_mc_lin_hp.o: denspol_mc_lin.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__HP_ENABLED__

denspol_mc_ring_topo.o: denspol_mc_ring.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__TOPO_MOVE_ENABLED__ 

denspol_mc_lin_topo.o: denspol_mc_lin.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__TOPO_MOVE_ENABLED__

denspol_mc_ring_hp_topo.o: denspol_mc_ring.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__HP_ENABLED__ -D__TOPO_MOVE_ENABLED__

denspol_mc_lin_hp_topo.o: denspol_mc_lin.c
	$(CC) -c $^ -o $@ $(CFLAGS) $(EXTRA_CFLAGS) -D__HP_ENABLED__ -D__TOPO_MOVE_ENABLED__



ee: ee.o lattice.o

# denspol: denspol.o denspol_init.o denspol_mc.o  denspol_lib.o denspol_io.o timer.o denspol_dbg.o file.o

denspol_ring: $(RING_OBJ) denspol_mc_ring.o

denspol_lin: $(LIN_OBJ) denspol_mc_lin.o 

denspol_ring_hp: $(RING_OBJ) denspol_mc_ring_hp.o
	$(CC) $(LDFLAGS) $^ -o $@ 

denspol_lin_hp: $(LIN_OBJ) denspol_mc_lin_hp.o
	$(CC) $(LDFLAGS) $^ -o $@ 

denspol_ring_topo: $(RING_OBJ) denspol_mc_ring_topo.o 
	$(CC) $(LDFLAGS) $^ -o $@ 

denspol_lin_topo: $(LIN_OBJ) denspol_mc_lin_topo.o 
	$(CC) $(LDFLAGS) $^ -o $@ 

denspol_ring_hp_topo: $(RING_OBJ) denspol_mc_ring_hp_topo.o
	$(CC) $(LDFLAGS) $^ -o $@ 

denspol_lin_hp_topo: $(LIN_OBJ) denspol_mc_lin_hp_topo.o
	$(CC) $(LDFLAGS) $^ -o $@ 

contact_map: contact_map.o denspol_lib_ring.o

clean:
	rm -f $(EXECS) *.o *~ 

install:       
	cp $(EXECS) $(BIN_DIR)