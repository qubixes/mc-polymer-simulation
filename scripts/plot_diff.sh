#!/bin/bash
. util.sh
TYPE=ring
DIRS=(`get_dirs $* -e`)
# LENGTHS=(`./get_lengths.sh`)
# TYPE=`./get_type.sh`
# ODIR=`./paper_dir.sh`
OUT_FILE="./dif_coef_$TYPE.dat"
# OFILE="dif_coef_$TYPE"

if [ $# -eq 1 ] && [ $1 == "-np" ]; then
	NOPLOT=1
else	
	NOPLOT=0
fi

rm -f $OUT_FILE

if [ $TYPE == "lin" ]; then
	LINES="0.06*x**-1.2 title \"\$0.06 N^{-1.2}\$\" lw 2, 58*x**-2.5 title \"\$32 N^{-2.4}\$\" lw 2"
	XSTART=50
else
	LINES='0.072*x**-1.2 title "\$0.072 N^{-1.2}\$" lw 2, 210*x**-2.4 title "\$55 N^{-2.2}\$" lw 2'
	XSTART=30
fi

I=0
for DIR in ${DIRS[*]}; do
	N=`get_attr 'Length' $DIR`
	D1=`
	awk '
		BEGIN{minim=1e300}{
		x=$(2)/$(1); 
		minim=(minim<x?minim:x);
		}
		END{
			print '$N', minim;
		}' "$DIR/cmsdif.dat"`
	D2=`../bin/get_diff $DIR/cmsdif.dat`
	echo $D1 $D2 >> $OUT_FILE
	let "I=I+1"
done


if [ $NOPLOT -eq 0 ]; then
gnuplot << EOF
set term aqua enhanced font 'Helvetica, 22'
set log x
set log y
set format y "10^{%T}"
set ylabel "D"
set xlabel "N"
plot [$XSTART:][:0.01] "$OUT_FILE" u 1:(\$3) w l title "MC Data" lw 2, $LINES
EOF

fi