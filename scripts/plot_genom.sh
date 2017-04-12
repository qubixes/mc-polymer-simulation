#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $*`) || { echo ${DIRS[*]}; exit $?; }

# PLOT="plot "
# 
# for DIR in ${DIRS[*]}; do
# 	PLOT="$PLOT \"$DIR/rgyr_time.dat\" w l, "
# done

PLOT="plot [:][:20]"
for DIR in ${DIRS[*]}; do
	PLOT="$PLOT \"$DIR/genom.dat\" u 1:(2**0.5*1.2*\$2**1.5/(\$1)) w l notitle, "
done


PLOT=${PLOT:0:${#PLOT}-2}

gnuplot << EOFGNU
set grid
set log x
# set log y
$PLOT, 11.8
EOFGNU