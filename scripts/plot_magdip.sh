#!/bin/bash

# NOTE These values are for the denspol algorithm. 

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
	MAGDIP_SIG=${MAGDIP[2]}
	RGYR=(`cat $DIR/rgyr.dat`)
	RGYR_SIG=${RGYR[1]}
	RGYR=${RGYR[0]}
	
	N=`get_attr 'Length' $DIR`
	echo "$N $MAGDIP $MAGDIP_SIG $RGYR $RGYR_SIG" >> $MAG_FILE
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
set key top left
plot [1:] "$MAG_FILE" u (\$1/$NE):(\$2/($LE*$LK)*(12./(2*$PI))):(\$3/($LE*$LK)*(12./(2*$PI))) with errorbars title "r_m" , 2.5*x**(2./3.), "$MAG_FILE" u (\$1/$NE):(\$4/($LE*$LK)*(12.)):(\$5/(($LE*$LK)*(12.))) with errorbars title "r_g" , 0.95*x**(2./3.)
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