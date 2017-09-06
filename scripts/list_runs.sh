#!/bin/bash

. util.sh

DIRS=(`get_dirs $*`)

I=0
printf "%6s %12s %12s %8s %12s\n" "N" "Ttot" "Teq" "Interval" "TMeasured"
echo "------------------------------------------------------"
for DIR in ${DIRS[*]}; do
	N=`get_attr 'Length' $DIR`
	INTER=`get_attr 'Interval' $DIR`
	TTOT=(`tail -1 $DIR/slrat.dat`)
	TEQED=(`tail -1 $DIR/ree.dat`)
	let "TEQ=TTOT-TEQED"
	printf "%6s %12s %12s %8s %12s\n" $N $TTOT $TEQ $INTER $TEQED
# 	echo $N $TTOT $TEQ $INTER $TEQED
done

