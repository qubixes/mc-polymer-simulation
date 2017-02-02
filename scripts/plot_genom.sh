#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $*`) || { echo ${DIRS[*]}; exit $?; }

PLOT="plot "

for DIR in ${DIRS[*]}; do
	PLOT="$PLOT \"$DIR/rgyr_time.dat\" w l, "
done

PLOT=${PLOT:0:${#PLOT}-2}

gnuplot << EOFGNU
set log x
# set log y
$PLOT
EOFGNU