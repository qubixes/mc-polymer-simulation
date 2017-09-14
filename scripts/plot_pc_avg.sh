#!/bin/bash

. util.sh

DIRS=`get_dirs $*`
BETA=1
PLOT="plot [][]"
PLOT2="plot [][0.3:3]"
FNAME="pc_vs_genom.dat"
ODIR=`get_paper_dir`
OFILE="pc_resc"


for DIR in ${DIRS[*]}; do
	FOUT="$DIR/$FNAME"
	cross_sect_files $DIR/pc_avg.dat $DIR/genom.dat > $FOUT
	PLOT="$PLOT '$DIR/pc_avg.dat' u 1:(\$2*\$1**$BETA) w l notitle, "
	PLOT2="$PLOT2 '$FOUT' u 2:(\$3/\$2**-1.5) w l notitle, "
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
set terminal aqua enhanced dashed font 'Times Roman, 26'
set log x
set log y
set ylabel 'p_c <|r|^2>^{3/2}'
set xlabel '<r^2>'
set label "x^{0.01}" at 100,1.85
$PLOT2, (x>20)?(1.65*x**0.01):1/0 lw 2 dt 2 lc rgb "black" notitle 
EOFGNU


# gnuplot << EOFGNU
# set term epslatex color standalone colortext 14
# set out "$OFILE.tex"
# 
# set log x 
# set log y
# set ylabel '\$ p_c \left<|\\mathbf{r}|^2\right>^{3/2} \$'
# set xlabel '\$\left<|\\mathbf{r}|^2\right>\$'
# set label '{\\tiny \$x^{0.01}\$ }' at 100,1.85
# $PLOT2, (x>20)?(1.65*x**0.01):1/0 lw 2 dt 2 lc rgb "black" notitle 
# EOFGNU
# 
# make $OFILE.eps && cp $OFILE.eps $ODIR
