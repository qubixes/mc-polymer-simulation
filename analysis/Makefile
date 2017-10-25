.PHONY: clean all install

VPATH=../util

EXECS=lowm_modes create_ptl avg_data localization merge_tseries get_diff create_cms cms_cor ptl2t scaleup denspol_scaleup
EXTRA_INC=../util
BIN_DIR=../bin
CFLAGS= -Wall -O2
LDFLAGS= -lm


$(EXECS): override CFLAGS+= -I$(EXTRA_INC) -std=c99 -DRELEASE=$(RELEASE)
all: $(EXECS)

lowm_modes: lowm_modes.o mode_compute.o mode_io.o timer.o mode_init.o mode_util.o file.o

create_ptl: create_ptl.o timer.o file.o

create_cms: create_cms.o timer.o file.o mode_init.o mode_io.o mode_util.o

cms_cor: cms_cor.o timer.o file.o mode_init.o mode_util.o mode_io.o

ptl2t: ptl2t.o timer.o file.o mode_init.o mode_util.o mode_io.o

scaleup: scaleup.o timer.o file.o mode_util.o

denspol_scaleup: denspol_scaleup.o mode_util.o

install:
	cp $(EXECS) $(BIN_DIR)

clean:
	rm -f $(EXECS) *.o *~ 
