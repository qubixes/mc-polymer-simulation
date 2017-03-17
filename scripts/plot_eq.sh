#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }
PLOT="plot [:6][] "
for DIR in ${DIRS[*]}; do
	DIR="$DIR/long/"
	F1="$DIR/rgyr_time.dat"
	F2="$DIR/cmsdif.dat"
	RG=`cat "$DIR/rgyr.dat"`
	FOUT="$DIR/cms_v_rgyr.dat"
	cross_sect_files $F1 $F2 > $FOUT
	PLOT="$PLOT \"$FOUT\" u (\$2/$RG):(1-\$3/$RG) w l title \"`get_title $DIR "N=%n"`\", "
	
# 	cat $FOUT
done

PLOT=${PLOT:0:${#PLOT}-2}

gnuplot -persist <<EOFGNU
# set log x
set key bottom right
set log y
$PLOT, 0.1*exp(-x)
EOFGNU