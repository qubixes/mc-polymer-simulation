#!/bin/bash

. util.sh

DIR=$1
ALPHA=1.0

LENGTH=$1
DIR=`get_dirs $*`
TYPE="ring"

if [ ! -d "$DIR" ] || [ $# -eq 0 ]; then echo "Need valid directory ($DIR)!"; exit 0; fi
FILE="$DIR/spac_dif.dat"

N=(`head -1 $FILE`)
echo "${N[*]}"
SCHIL=4
# PEELN=("30" "100" "150" "300" "400" "500" "1000" "1500" "2000" "3000" "4000")
# PEELD=("1.5" "2.5" "2.8" "3.7" "4.2" "4.4"  "6.5" "9.6"  "12.4" "15.5" "17.4")

I=0
while read X Y; do
PEELN[I]=$X
PEELD[I]=$Y
(( I=I+1 ))
done < spac_dif_${TYPE}.dat

# for I in `seq 0 1 ${#PEELD[*]}`; do
# 	echo ${PEELN[I]} ${PEELD[I]};
# done

L=`get_attr 'Length' $DIR`

I=0
while [ $I -lt ${#PEELN[*]} ]; do
	if [ $L -eq ${PEELN[I]} ]; then
		SCHIL=${PEELD[I]}
		break;
	fi
	let "I=I+1"
done

let "NCOL=${#N[*]}"
let "BS_MAX=N[1]"

PLOT="plot [][:1.0]"
PLOT2="plot "

for i in `seq 2 $NCOL`; do 
	R=`echo "print (3*${N[i-1]}/(4*pi))**(1./3.)" | gnuplot 2> /dev/stdout`
	if [ $BS_MAX -eq ${N[i-1]} ]; then
		COR_FAC="(1.0)"
	else
		COR_FAC="(($BS_MAX)/(1.0*($BS_MAX-${N[i-1]})))"
	fi
	
	PLOT="$PLOT \"$FILE\" u 1:(\$$i*$COR_FAC*((3*${N[i-1]}/(4*pi))**(1./3.)+$SCHIL)**3/\$1) w l title 'r=$R , d=$SCHIL'"
	PLOT2="$PLOT2 \"$FILE\" u 1:(\$$i/\$1) w l "
	if [ $i -ne $NCOL ]; then 
		PLOT="$PLOT, "; 
		PLOT2="$PLOT2, "; 
	fi
done

gnuplot -persist <<EOFGNU
set log x
set log y
set title "N=${L}"
$PLOT
EOFGNU


# gnuplot -persist <<EOFGNU
# set log x
# set log y
# 
# $PLOT2
# EOFGNU
# 
