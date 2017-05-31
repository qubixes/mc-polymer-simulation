#!/bin/bash

EXECS=("efpol" "gpupol2" "gpupol3" "denspol")
PLOT="plot "
PLOT2="plot "

for EXEC in ${EXECS[*]}; do
	FILE=${EXEC}_ne_list.dat
	grep "$EXEC" ne_list.dat > $FILE
	PLOT="$PLOT \"$FILE\" u 2:(\$3**-3.2/\$5) title \"$EXEC\", "
	PLOT2="$PLOT2 \"$FILE\" u 2:(\$3**-3.2/\$5*\$6) title \"$EXEC\", "
done

PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}


gnuplot -persist <<EOFGNU
set log y
$PLOT
EOFGNU

gnuplot -persist <<EOFGNU
set log y
set key bottom right
set term aqua enhanced font 'Times Roman, 20'
set xlabel '{/Symbol r}'
set ylabel 'Efficiency= N_e^{-3.2} t_{fac}^{-1} M/s 
$PLOT2
EOFGNU