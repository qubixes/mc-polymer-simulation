#!/bin/bash

. util.sh

DIRS=`get_dirs $*`
BETA=1
PLOT="plot [][]"
PLOT2="plot [][]"
FNAME="pc_vs_genom.dat"

for DIR in ${DIRS[*]}; do
	FOUT="$DIR/$FNAME"
	cross_sect_files $DIR/pc_avg.dat $DIR/genom.dat > $FOUT
	PLOT="$PLOT '$DIR/pc_avg.dat' u 1:(\$2*\$1**$BETA) w l title '`get_title $DIR "N=%n"`', "
	PLOT2="$PLOT2 '$FOUT' u 2:(\$3/\$2**-1.5) w l title '`get_title $DIR "N=%n"`', "
done

PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}

# echo $PLOT
gnuplot -persist <<EOFGNU
set grid
set log x
set log y
set ylabel 'p_c'
set xlabel '|i-j|'
$PLOT, 1.02*x**-(1.3-$BETA)
EOFGNU

gnuplot -persist <<EOFGNU
set log x
set log y
set ylabel 'p_c'
set xlabel 'r^2'
$PLOT2, 1.5*x**0.02
EOFGNU
