.PHONY: all clean install
FC = gfortran-mp-4.9
FFLAGS=-O2 
CFLAGS=-O2 -Wall -std=c99
CC = gcc
EXTRA_INC= ../analysis ../util
VPATH= $(EXTRA_INC) 

BIN_DIR = ../bin
SRCS = main.f90 readparam.f90 readpos.f90 ener.f90 tW.f90 tV.f90
OBJS = main.o readparam.o readpos.o ener.o tW.o tV.o
PRGM = dJost.exe 
EXECS= $(PRGM) secstr_wrap

override CFLAGS += -I../analysis -I../util

all: $(EXECS)

$(PRGM): $(OBJS)
	$(FC) -o $(PRGM) $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

secstr_wrap: secstr_wrap.o mode_util.o
# %.o: %.c
# 	gcc -c $<

install:
	cp $(EXECS) input_template.dat $(BIN_DIR)

clean:
	rm -f $(PRGM) *.o
