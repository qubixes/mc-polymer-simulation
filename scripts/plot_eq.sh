#!/bin/bash

. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }
PLOT="plot [:6][] "
for DIR in ${DIRS[*]}; do
	DIR="$DIR/long/"
	F1="$DIR/rgyr_time.dat"
	F2="$DIR/cmsdif.dat"
	RG=`cat "$DIR/rgyr.dat"`
	FOUT="$DIR/cms_v_rgyr.dat"
	rm -f $FOUT
	if [ ! -f "$F1" -o ! -f "$F2" ]; then continue; fi
	
	exec 3<$F2
	exec 4<$F1
	
# 	head $F1
	read T1 R1 <&3
	read T2 R2 <&4
	while true; do
# 		echo "$T1 $T2"
		if [ "$T1" == "$T2" ]; then 
			echo "$T2 $R1 $R2" >> $FOUT
			read T1 R1 <&3;
		elif [ "$T1" == "" -o "$T2" == "" ]; then break;
		fi
		while [ "$T1" != "" ] && [ $T1 -lt $T2 ]; do
			read T1 R1 <&3;
		done
# 		echo "$T1 $T2"
		if [ "$T1" == "$T2" ]; then 
			echo "$T2 $R1 $R2" >> $FOUT
			read T2 R2 <&4;
		elif [ "$T1" == "" -o "$T2" == "" ]; then break;
		fi
		
		while [ "$T2" != "" ] && [ $T2 -lt $T1 ]; do
			read T2 R2 <&4;
		done
	done
	exec 3<&-
	exec 4<&-
	PLOT="$PLOT \"$FOUT\" u (\$2/$RG):(1-\$3/$RG) w l title \"`get_title $DIR "N=%n"`\", "
	
# 	cat $FOUT
done

PLOT=${PLOT:0:${#PLOT}-2}

gnuplot -persist <<EOFGNU
# set log x
set key bottom right
set log y
$PLOT, 0.1*exp(-x)
EOFGNU