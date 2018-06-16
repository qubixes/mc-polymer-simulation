#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $*`) || { echo ${DIRS[*]}; exit $?; }

# PLOT="plot "
# 
# for DIR in ${DIRS[*]}; do
# 	PLOT="$PLOT \"$DIR/rgyr_time.dat\" w l, "
# done

TMP_FILE="ne_est_`echo $* | sed 's/ /_/g' | sed 's/-/_/g'`.tmp"
rm -f $TMP_FILE

function ne_estimate {
	DIR=$1
	
	FILE="$DIR/genom.dat"
	
	awk 'function pow(x,y) {return exp(y*log(x))};
		BEGIN{
			max_gauss=0; 
			max_compact=0;
		}{
			g=$(1); rsq=$(2);
			if(rsq/g > max_gauss)
				max_gauss = rsq/g;
			if(rsq/pow(g,2./3.) > max_compact)
				max_compact = rsq/pow(g,2./3.);
		}END{ print max_compact, max_gauss, pow(max_compact/max_gauss, 3);}
	' $FILE
}


PLOT="plot [:][:20]"
PLOT2="plot "
for DIR in ${DIRS[*]}; do
	DENSITY=`get_attr 'Density' $DIR`
	echo `get_attr 'Length' $DIR` `ne_estimate $DIR` >> $TMP_FILE
	PLOT="$PLOT \"$DIR/genom.dat\" u 1:(2**0.5*$DENSITY*\$2**1.5/(\$1)) w l notitle, "
	PLOT2="$PLOT2 \"$DIR/genom.dat\" u 1:2 w l notitle, " 
done


PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}

gnuplot << EOFGNU
set grid
set log x
# set log y
$PLOT, 11.8
EOFGNU

gnuplot << EOFGNU
set grid
set log x
set log y
$PLOT2, 0.59*x, 1.296*x**(2./3.)
EOFGNU

gnuplot << EOFGNU

# set log x
# set log y
plot "$TMP_FILE" u (\$1**-0.3):4 w l, 16.4-35*x

EOFGNU
# 10.8
# 10.4