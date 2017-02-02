.PHONY: clean install

EXECS=bench
BIN_DIR=../bin
CFLAGS= -Wall -O2 

$(EXECS): override CFLAGS+= -std=c99
all: $(EXECS)

bench: bench.o timer.o

clean:
	rm -f $(EXECS) *.o *~ 

install: 
	cp $(EXECS) $(BIN_DIR)