#!/bin/bash

EXECS=("efpol" "gpupol2" "gpupol3" "denspol")
PLOT="plot "
PLOT2="plot "

for EXEC in ${EXECS[*]}; do
	FILE=${EXEC}_ne_list.dat
	grep "$EXEC" ne_list.dat > $FILE
	PLOT="$PLOT \"$FILE\" u 2:(\$3**3.0/\$5) title \"$EXEC\", "
	PLOT2="$PLOT2 \"$FILE\" u 2:(\$3**3.0/\$5*\$6) title \"$EXEC\", "
done

PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}


gnuplot -persist <<EOFGNU
set log y
$PLOT
EOFGNU

gnuplot -persist <<EOFGNU
set log y
$PLOT2
EOFGNU