#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $*`) || { echo ${DIRS[*]}; exit $?; }

# PLOT="plot "
# 
# for DIR in ${DIRS[*]}; do
# 	PLOT="$PLOT \"$DIR/rgyr_time.dat\" w l, "
# done

PLOT2="plot [][-2:2]"
for DIR in ${DIRS[*]}; do
	N=`get_attr 'Length' $DIR`
	PLOT2="$PLOT2 \"$DIR/ucor_avg.dat\" u 1:(\$2*$N) w l notitle, " 
done


PLOT2=${PLOT2:0:${#PLOT2}-2}

# gnuplot << EOFGNU
# set grid
# set log x
# # set log y
# $PLOT, 11.8
# EOFGNU

gnuplot << EOFGNU
set grid
set log x
# set log y
$PLOT2
EOFGNU