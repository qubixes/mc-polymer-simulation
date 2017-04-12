.PHONY: clean install

EXECS= ee
BIN_DIR=../bin
VPATH=../util
CFLAGS= -Wall -O2 
# CFLAGS=-g -std=c99 -Wall
LDFLAGS=-lpthread -lm
EXTRA_INC=../util
$(EXECS): override CFLAGS+= -I$(EXTRA_INC) -std=c99

all: $(EXECS)

ee: ee.o lattice.o


clean:
	rm -f $(EXECS) *.o *~ 

install:       
	cp $(EXECS) $(BIN_DIR)