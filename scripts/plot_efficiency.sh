#!/bin/bash
#  "denspol"
EXECS=("efpol" "gpupol2" "gpupol3" "denspol")
TITLES=("CPU - standard" "GPU" "GPU" "CPU - denspol")
PLOT="plot "
PLOT2="plot "

I=0
for EXEC in ${EXECS[*]}; do
	FILE=${EXEC}_ne_list.tmp
	grep "$EXEC" ne_list.dat > $FILE
	if [ $EXEC == "gpupol2" -o $EXEC == "gpupol3" ]; then
		PLOT="$PLOT \"$FILE\" u 2:(\$3**-3.2/\$5) title \"${TITLES[I]}\", "
		PLOT2="$PLOT2 \"$FILE\" u 2:(\$3**-3.2/\$5*\$6*\$2*1e6) title \"${TITLES[I]}\", "
	else
		PLOT="$PLOT \"$FILE\" u 2:(\$3**-3.2/\$5) title \"${TITLES[I]}\", "
		PLOT2="$PLOT2 \"$FILE\" u 2:(\$3**-3.2/\$5*\$6*1e6) title \"${TITLES[I]}\", "
	fi
	let "I=I+1"
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
set term aqua enhanced font 'Times Roman, 24'
set xlabel '{/Symbol r}'
set ylabel 'Efficiency= N_e^{-3.2} t_{fac}^{-1} M/s 
$PLOT2
EOFGNU