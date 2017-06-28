#!/bin/bash

SIM=$*

NPOL=100
LINE=(`grep ^"$SIM" "ne_list.dat"`)
NE=${LINE[2]}
TF=${LINE[4]}
MS=${LINE[5]}

# NE=$1
# TF=$2


gnuplot -persist <<EOFGNU
trouse(Ne, tFac, N) = 1.4*(N/Ne)**2.2*tFac*Ne**2.2
trept(Ne, tFac, N) = 0.261*(N/Ne)**3.08*tFac*Ne**2.2

set log x
set log y
set xlabel "N"
set ylabel "t_{eq}"
max(x,y)=(x>y)?x:y
print $NE, $TF
plot [10:10000] max(trouse($NE, $TF, x), trept($NE, $TF, x))
EOFGNU

gnuplot -persist <<EOFGNU
trouse(Ne, tFac, N) = 1.4*(N/Ne)**2.2*tFac*Ne**2.2
trept(Ne, tFac, N) = 0.261*(N/Ne)**3.08*tFac*Ne**2.2

set log x
set log y
set xlabel "N/N_e"
set ylabel "t_{eq} (hours)"
max(x,y)=(x>y)?x:y
print $NE, $TF
real_time(x) = max(trouse($NE, $TF, x), trept($NE, $TF, x))*$NPOL*x/($MS*1e6*3600)
plot [0.1:100] real_time(x*$NE)
EOFGNU