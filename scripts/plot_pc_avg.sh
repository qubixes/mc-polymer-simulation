#!/bin/bash

. util.sh

DIRS=`get_dirs $*`
BETA=1.17
PLOT="plot [][0.5:3]"

for DIR in ${DIRS[*]}; do
	PLOT="$PLOT '$DIR/pc_avg.dat' u 1:(\$2*\$1**$BETA) w l , "
done

PLOT=${PLOT:0:${#PLOT}-2}

# echo $PLOT
gnuplot -persist <<EOFGNU
set log x
set log y
$PLOT, 1.02*x**-(1.3-$BETA)
EOFGNU