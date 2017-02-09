.PHONY: clean all install

VPATH=../util

EXECS=lowm_modes create_ptl avg_data localization merge_tseries get_diff
EXTRA_INC=../util
BIN_DIR=../bin
CFLAGS= -Wall -O2
LDFLAGS= -lm


$(EXECS): override CFLAGS+= -I$(EXTRA_INC) -std=c99
all: $(EXECS)

lowm_modes: lowm_modes.o mode_compute.o mode_io.o timer.o mode_init.o mode_util.o file.o

create_ptl: create_ptl.o timer.o file.o



install:
	cp $(EXECS) $(BIN_DIR)

clean:
	rm -f $(EXECS) *.o *~ 
