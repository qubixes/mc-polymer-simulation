#!/bin/bash

. util.sh

DIRS=(`get_dirs $*`)
TFILE="latex_table.tmp"
rm -f $TFILE
I=0
printf "%6s %6s %12s %12s %8s %12s\n" "N" "M" "Ttot" "Teq" "Interval" "TMeasured"
echo "-------------------------------------------------------------"
for DIR in ${DIRS[*]}; do
	N=`get_attr 'Length' $DIR`
	NPOL=`get_attr 'Npol' $DIR`
	INTER=`get_attr 'Interval' $DIR`
	TTOT=(`tail -1 $DIR/slrat.dat`)
	TEQED=(`tail -1 $DIR/ree.dat`)
	let "TEQ=TTOT-TEQED"
	printf "%6s %6s %12s %12s %8s %12s\n" $N $NPOL $TTOT $TEQ $INTER $TEQED
	Z=`echo "print $N/121." | gnuplot 2> /dev/stdout` 
	Z=`printf "%.2lf " "$Z"`
	NTAU=`echo "print $TTOT.0/($TEQ.0)" | gnuplot 2> /dev/stdout`
	NTAU=`printf "%.1lf" "$NTAU"`
	
# 	echo $TEQ `get_teq $DIR` >> $TFILE
	TEQ_THEO=`get_teq $DIR`
	TEQ_THEO=`printf "%le" "$TEQ_THEO"`
	NTAU_THEO=`echo "print $TTOT.0/(1.0*$TEQ_THEO)" | gnuplot 2> /dev/stdout`
	NTAU_THEO=`printf "%.1le" "$NTAU_THEO"`
	DECADES=`echo "print $NTAU_THEO*$NPOL" | gnuplot 2> /dev/stdout`
	DECADES=`printf "%.1le" $DECADES`
	echo '{\footnotesize '$Z'} &  {\footnotesize $'$N'\times '$NPOL'$} & {\footnotesize $'$TTOT'$} & {\footnotesize $'$TEQ_THEO'$} & {\footnotesize $'$NTAU_THEO'$} & {\footnotesize $'$DECADES'$}\\' >> $TFILE
# 	echo $N $TTOT $TEQ $INTER $TEQED
done

cat $TFILE
