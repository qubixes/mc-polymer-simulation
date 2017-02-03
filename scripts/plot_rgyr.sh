#!/bin/bash

. ./util.sh
DIRS=(`get_dirs $* -e`) || { echo ${DIRS[*]}; exit $?; }

DATA_DIR="../data"
RGYR_FILE="$DATA_DIR/rgyr`echo "$*" | sed 's/ //g' | sed 's/-/_/g'`.dat"
rm -f $RGYR_FILE

I=0


for DIR in ${DIRS[*]}; do
	RGYR=(`cat $DIR/rgyr.dat`)
	N=`get_attr 'Length' $DIR`
	NP=`ls $DIR/ptl | wc -w`
	echo "$N $RGYR" >> $RGYR_FILE
done

gnuplot << EOF
set term aqua enhanced font 'Helvetica, 22'
set log x
set log y
set xlabel "N"
set ylabel "r^2"
set key bottom right
plot "$RGYR_FILE" u 1:(\$2/\$1) w l lt 1 lw 2
EOF