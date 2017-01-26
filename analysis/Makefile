.PHONY: clean all install

VPATH=../util

EXECS=lowm_modes create_ptl
EXTRA_INC=../util
BIN_DIR=../bin
CFLAGS= -Wall -O2 -std=c99
LDFLAGS= -lm

$(EXECS): override CFLAGS+= -I$(EXTRA_INC)
all: $(EXECS)

lowm_modes: lowm_modes.o mode_compute.o mode_io.o timer.o mode_init.o mode_util.o

create_ptl: create_ptl.o timer.o



install:
	cp $(EXECS) $(BIN_DIR)

clean:
	rm -f $(EXECS) *.o *~ 
