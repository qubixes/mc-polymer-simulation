#!/bin/bash
. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }
TYPE="ring"
EXEC=`get_attr 'Executable' $DIRS`
ODIR=`get_paper_dir`
OFILE=cmsperp_$TYPE


PLOT="plot [:][0.002:0.1]"

I=0

LT=0
for DIR in ${DIRS[*]}; do
	if [ ! -f "$DIR/cmsdif.dat" ]; then continue; fi
	N=`get_attr 'Length' $DIR`
	
	if [ $LT -eq 0 ]; then 
		TITLE1="title 'cms_{perp}'"
		TITLE2="title 'cms_{par}'"
		TITLE3="title 'cms\_v\_mag_{perp}'"
		TITLE4="title 'cms\_v\_mag_{par}'"
	else
		TITLE1="notitle"
		TITLE2="notitle"
		TITLE3="notitle"
		TITLE4="notitle"
	fi
	
	PLOT="${PLOT} \"$DIR/cms_dir_mag.dat\" u (\$1/f($N)):(\$2/g($N)/h(\$1)) w l lt 1 $TITLE1, "
	PLOT="${PLOT} \"$DIR/cms_dir_mag.dat\" u (\$1/f($N)):(\$3/g($N)/h(\$1)) w l lt 2 $TITLE2, "
	PLOT="${PLOT} \"$DIR/cms_dir_mag.dat\" u (\$1/f($N)):(\$5/g($N)/h(\$1)) w l lt 3 $TITLE3, "
	PLOT="${PLOT} \"$DIR/cms_dir_mag.dat\" u (\$1/f($N)):(\$6/g($N)/h(\$1)) w l lt 4 $TITLE4, "
	let "I=I+1"
	let "LT++"
done

PLOT=${PLOT:0:${#PLOT}-2}


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
