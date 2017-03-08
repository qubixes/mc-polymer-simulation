#!/bin/bash

. util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }

COL=2
PLOT="plot "

for DIR in ${DIRS[*]}; do
	DIR_SH=`echo $DIR/short*`
	DIR_LO=`echo $DIR/long`
	PLOT="$PLOT '$DIR_SH/cms_cor.dat' u 1:$COL w l, "
done

PLOT=${PLOT:0:${#PLOT}-2}

gnuplot -persist <<EOFGNU
set log x
# set log y
$PLOT
EOFGNU