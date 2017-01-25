.PHONY: clean install

EXECS=bench
BIN_DIR=../bin
CFLAGS= -Wall -O2 -std=c99
# CFLAGS=-g -std=c99 -Wall
all: $(EXECS)

bench: bench.o timer.o

clean:
	rm -f $(EXECS) *.o *~ 

install: 
	cp $(EXECS) $(BIN_DIR)