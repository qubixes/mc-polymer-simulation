#!/bin/bash
. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }
TYPE="ring"
EXEC=`get_attr 'Executable' $DIRS`
ODIR=`get_paper_dir`
OFILE=cmsdif_$TYPE


if [ $TYPE == "ring" ]; then
	if [ "$EXEC" == "denspol" ]; then
		EX_FUNC="i(x) = (x>1e5 && x<3e7)?1*x**-0.45:1/0; j(x) = (x>1e2 && x<7e4)?0.15*x**-0.3:1/0"
		EX_PLOT="i(x) title 't^{-0.45}' lc 1 lw 2, j(x) title 't^{-0.3}' lc 2 lw 2"
	else
		EX_FUNC="i(x) = (x>2e6 && x<3e8)?4.0*x**-0.35:1/0; j(x) = (x>1e2 && x<7e4)?0.21*x**-0.13:1/0"
		EX_PLOT="i(x) title 't^{-0.35}' lc 1 lw 2, j(x) title 't^{-0.13}' lc 2 lw 2"
	fi
	EX_PLOT_FILE="i(x) lc rgb \"black\" lw 2 dt 2 notitle, j(x) lc rgb \"black\" lw 2 dt 2 notitle"
else
	EX_PLOT="0.6*x**-0.16 title '\$0.6 t^{-0.16}\$' lc 2 lw 2"
fi
PLOT="plot [:][0.0005:0.1]"

I=0
for DIR in ${DIRS[*]}; do
	if [ ! -f "$DIR/cmsdif.dat" ]; then continue; fi
	N=`get_attr 'Length' $DIR`
	PLOT="${PLOT} \"$DIR/cmsdif.dat\" u (\$1/f($N)):(\$2/g($N)/h(\$1)) w l notitle, "
	let "I=I+1"
done

gnuplot << EOFGNU
set term aqua enhanced font 'Times Roman, 22'
set log x
set log y
set yrange [:0.1]
set ylabel '<{r^2}_{cms}> N t^{-1}'
set xlabel 't'
f(x) = x**0
g(x) = x**-1
h(x) = x**1.0
$EX_FUNC
$PLOT $EX_PLOT
EOFGNU


# gnuplot << EOFGNU
# set term epslatex color standalone colortext 14
# set log x 
# set log y
# set format x "\$10^{%T}\$"
# set label '{\\tiny \$t^{-0.13}\$ }' at 800,0.063
# set label '{\\tiny \$t^{-0.35}\$ }' at 3e6,0.01
# set ylabel '\$\\left<{r^2}_{\\mathrm{cms}}\\right> N t^{-1}\$'
# set xlabel '\$t\$'
# set key spacing 1.5
# set out "$OFILE.tex"
# f(x) = x**0
# g(x) = x**-1
# h(x) = x**1.0
# $EX_FUNC
# $PLOT $EX_PLOT_FILE
# EOFGNU
# 
# make $OFILE.eps && cp $OFILE.eps $ODIR
