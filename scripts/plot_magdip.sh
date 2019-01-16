#!/bin/bash

# NOTE These values are for the denspol algorithm. 

LK="1.31"
LE="5.6"
NE="10.1"
PI="3.1415926536"
EXP="3.5"
MIN_N="1000"


MAG_FILE="magdip.tmp"
rm -f $MAG_FILE

. ./util.sh
DIRS=(`get_dirs $* -e`) || { echo ${DIRS[*]}; exit $?; }

I=0

PLOT="plot [:2e7]"
PLOT2="plot [:2e7]"

FRAC=0.5

for DIR in ${DIRS[*]}; do
	MAGDIP=(`head -1 $DIR/magdip.dat`)
	MAGDIP_SQ=(`head -2 $DIR/magdip.dat | tail -1`)
	MAGDIP_SQ=${MAGDIP_SQ[1]}
	MAGDIP=${MAGDIP[1]}
	MAGDIP_SIG=${MAGDIP[2]}
	RGYR=(`cat $DIR/rgyr.dat`)
	RGYR_SIG=${RGYR[1]}
	RGYR=${RGYR[0]}
	
	N=`get_attr 'Length' $DIR`
	if [ $N -lt $MIN_N ]; then continue; fi;
	echo "$N $MAGDIP $MAGDIP_SIG $RGYR $RGYR_SIG" >> $MAG_FILE
	LAST_T=`awk '{if(!found && $(2)<0){found=1; last_t=$(1);}}END{print last_t;}' "$DIR/magdip.dat"`
	PLOT="$PLOT \"$DIR/magdip.dat\" u (\$1*($N/1000.0)**-$EXP):((\$1>$LAST_T)?(1/0):(\$2/$MAGDIP_SQ)) title \"N=$N\", "
	PLOT2="$PLOT2 \"$DIR/magdip.dat\" u (\$1*($N/1000.0)**-3.5):(\$2/$MAGDIP_SQ - 0.92 * ((1.0-$FRAC)*exp(-\$1/(8*$N)) + $FRAC*exp(-\$1/(0.55*$N**2)))) w l title \"N=$N\", "

done

PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}

# gnuplot << EOF
# set term aqua enhanced font 'Helvetica, 22'
# set log x
# set log y
# set xlabel "N/N_e"
# set ylabel "r_m^2/(L_e*l_k/12)"
# set key top left
# plot [1:] "$MAG_FILE" u (\$1/$NE):(\$2/($LE*$LK)*(12./(2*$PI))):(\$3/($LE*$LK)*(12./(2*$PI))) with errorbars title "r_m" , 2.5*x**(2./3.), "$MAG_FILE" u (\$1/$NE):(\$4/($LE*$LK)*(12.)):(\$5/(($LE*$LK)*(12.))) with errorbars title "r_g" , 0.95*x**(2./3.)
# EOF

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
gnuplot <<EOFGNU
set term aqua enhanced font 'Helvetica, 15'
set log y
set xlabel "t/N^{3.5}"
set ylabel "<m(t)*m(t+dt)>/N^{4/3}"
remove_lt0(x,zero) = (zero>0)?x:(1/0)
$PLOT, 0.06*exp(-x/2.5e6)
EOFGNU
# 
# gnuplot <<EOFGNU
# set term aqua enhanced font 'Helvetica, 22'
# set log y
# # set log x
# set key bottom left
# set xlabel "t"
# set ylabel "<m(t)*m(t+dt)>"
# $PLOT2
# EOFGNU