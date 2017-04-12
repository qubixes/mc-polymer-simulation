#!/bin/bash
. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }
TYPE="ring"

if [ $TYPE == "ring" ]; then
	EX_PLOT="2.6*x**-0.3 title '\$5 t^{-0.3}\$' lc 1 lw 2, 0.25*x**-0.13 title '\$0.48 t^{-0.13}\$' lc 2 lw 2"
else
	EX_PLOT="0.6*x**-0.16 title '\$0.6 t^{-0.16}\$' lc 2 lw 2"
fi
PLOT="plot [:][0.001:1]"

I=0
for DIR in ${DIRS[*]}; do
	if [ ! -f "$DIR/cmsdif.dat" ]; then continue; fi
	N=`get_attr 'Length' $DIR`
	PLOT="${PLOT} \"$DIR/cmsdif.dat\" u (\$1/f($N)):(\$2/g($N)/h(\$1)) w l notitle, "
	let "I=I+1"
done

gnuplot << EOFGNU
set term aqua enhanced font 'Helvetica, 22'
set log x
set log y
set ylabel '<{r^2}_{cms}> N t^{-1}'
set xlabel 't'
f(x) = x**0
g(x) = x**-1
h(x) = x**1.0
$PLOT $EX_PLOT
EOFGNU
