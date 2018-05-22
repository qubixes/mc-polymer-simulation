#!/bin/bash

LK="1.31"
LE="5.6"
NE="10.1"
PI="3.1415926536"


MAG_FILE="magdip.tmp"
rm -f $MAG_FILE

. ./util.sh
DIRS=(`get_dirs $* -e`) || { echo ${DIRS[*]}; exit $?; }

I=0

PLOT="plot [:0.03]"
PLOT2="plot "

for DIR in ${DIRS[*]}; do
	MAGDIP=(`head -1 $DIR/magdip.dat`)
	MAGDIP=${MAGDIP[1]}
	RGYR=`cat $DIR/rgyr.dat`
	N=`get_attr 'Length' $DIR`
	echo "$N $MAGDIP $RGYR" >> $MAG_FILE
	PLOT="$PLOT \"$DIR/magdip.dat\" u (\$1*$N**-3.5):(\$2*$N**-(4./3.)) w l title \"N=$N\", "
	PLOT2="$PLOT2 \"$DIR/magdip.dat\" u 1:(\$2*$N**-(4./3.)) w l title \"N=$N\", "

done

PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}

gnuplot << EOF
set term aqua enhanced font 'Helvetica, 22'
set log x
set log y
set xlabel "N/N_e"
set ylabel "r_m^2/(L_e*l_k/12)"
set key top right
plot [1:100] "$MAG_FILE" u (\$1/$NE):(\$2/($LE*$LK)*(12./(2*$PI))) title "r_m" , 0.29*x**(2./3.), "$MAG_FILE" u (\$1/$NE):(\$3/($LE*$LK)*(12.)) title "r_g" , 0.8*x**(2./3.)
EOF
# 
# 
# gnuplot << EOF
# set term aqua enhanced font 'Helvetica, 22'
# set log x
# set log y
# set xlabel "r_g^2"
# set ylabel "<m>"
# set key top right
# plot "$MAG_FILE" u 2:3 w l, 0.46*x
# EOF
# 
# gnuplot <<EOFGNU
# set term aqua enhanced font 'Helvetica, 22'
# set log y
# set xlabel "t/N^{3.5}"
# set ylabel "<m(t)*m(t+dt)>/N^{4/3}"
# $PLOT
# EOFGNU
# 
# gnuplot <<EOFGNU
# set term aqua enhanced font 'Helvetica, 22'
# set log y
# set log x
# set key bottom left
# set xlabel "t"
# set ylabel "<m(t)*m(t+dt)>"
# $PLOT2
# EOFGNU