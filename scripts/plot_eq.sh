#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }
PLOT="plot [:6][] "
PLOT2="plot "
for DIR in ${DIRS[*]}; do
	DIR="$DIR/long/"
	F1="$DIR/rgyr_time.dat"
	F2="$DIR/cmsdif.dat"
	RG=`cat "$DIR/rgyr.dat"`
	FOUT="$DIR/cms_v_rgyr.dat"
	cross_sect_files $F1 $F2 > $FOUT
	TITLE=`get_title $DIR "N=%n"`
	PLOT="$PLOT \"$FOUT\" u (\$2/$RG):(1-\$3/$RG) w l title \"$TITLE\", "
	PLOT2="$PLOT2 \"$F1\" u 1:2 w l title \"$TITLE\", "
# 	cat $FOUT
done

PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}


gnuplot -persist <<EOFGNU
# set log x
set key bottom right
set log y
$PLOT, 0.1*exp(-x)
EOFGNU

gnuplot -persist <<EOFGNU
set log x
set key bottom right
set log y
$PLOT2
EOFGNU