#!/bin/bash

. ./util.sh
DIRS=(`get_dirs $*`) || { echo ${DIRS[*]}; exit $?; }


PLOT="plot "

for DIR in ${DIRS[*]}; do
	TITLE=`get_title $DIR "N=%n"`
	PLOT="$PLOT \"$DIR/slrat.dat\" w l title '$TITLE', "
done

PLOT=${PLOT:0:${#PLOT}-2}

gnuplot << EOFGNU
set log x
set ylabel 'N_{sl}/N'
set xlabel 't'
$PLOT
EOFGNU