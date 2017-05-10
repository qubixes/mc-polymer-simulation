#!/bin/bash
. util.sh
TYPE=ring
DIRS=(`get_dirs $* -e`)
# LENGTHS=(`./get_lengths.sh`)
# TYPE=`./get_type.sh`
# ODIR=`./paper_dir.sh`
# OUT_FILE="./dif_coef_$TYPE.dat"
# OFILE="dif_coef_$TYPE"

NFILES=0
FILES=()
EXECS=()
DENSITIES=()

if [ $# -eq 1 ] && [ $1 == "-np" ]; then
	NOPLOT=1
else	
	NOPLOT=0
fi

rm -f ./dif_coef*

if [ $TYPE == "lin" ]; then
	LINES="0.06*x**-1.2 title \"\$0.06 N^{-1.2}\$\" lw 2, 58*x**-2.5 title \"\$32 N^{-2.4}\$\" lw 2"
	XSTART=50
else
	LINES='0.0328*x**-1.2 title "\$0.072 N^{-1.2}\$" lw 2, 88*x**-2.4 title "\$55 N^{-2.2}\$" lw 2'
	XSTART=30
fi

I=0
for DIR in ${DIRS[*]}; do
	N=`get_attr 'Length' $DIR`
	RHO=`get_attr 'Density' $DIR`
	EXEC=`get_attr 'Executable' $DIR`
	FILE="./dif_coef_${EXEC}_${RHO}"
	NE_LINE=(`grep '^'"$EXEC ${RHO:0:3}" "ne_list.dat"`)
	NE_FAC=${NE_LINE[2]}
	R_FAC=${NE_LINE[3]}
	T_FAC=${NE_LINE[4]}
	if [ ! -f $FILE ]; then
		FILES[NFILES]=$FILE
		DENSITIES[NFILES]=$RHO
		EXECS[NFILES]=$EXEC
		let "NFILES++"
	fi
	
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
	RG=`cat "$DIR/rgyr.dat"`
	echo $D1 $D2 $RG $NE_FAC $R_FAC $T_FAC >> $FILE
	let "I=I+1"
done

PLOT="plot [:][:0.01]"
PLOT2="plot [1:50]"

I=0
for FILE in ${FILES[*]}; do
	BTITLE="{/Symbol r}=${DENSITIES[I]}, exec=${EXECS[I]}"
	PLOT="$PLOT \"$FILE\" u (\$1/\$7):(\$3*\$9*\$7**1.2) w l title \"$BTITLE\" lw 2, "
	PLOT2="$PLOT2 \"$FILE\" u (\$1/\$7):(7.8*\$6/\$3/6.0/\$9/\$7**2.2) w l title \"$BTITLE (Diffuse RG)\" lw 2, "
# 	PLOT2="$PLOT2 \"$FILE\" u 1:4 w l title \"$BTITLE (CMS Decay)\" lw 2, "
	let "I++"
done

if [ $NOPLOT -eq 0 ]; then

gnuplot << EOF
set term aqua enhanced font 'Helvetica, 22'
set log x
set log y
set format y "10^{%T}"
set ylabel "D"
set xlabel "N"
$PLOT $LINES
EOF

gnuplot << EOF
set term aqua enhanced font 'Helvetica, 18'
set log x
set log y
set format y "10^{%T}"
set ylabel "t"
set xlabel "N"
set key bottom right
alpha=5
max(x,y) = (x>y)?x:y;
f(x) = (alpha/exp(1.32))**(1./0.85)*(x)**2.21
g(x) = (alpha/exp(5.02))**(1./0.61)*(x)**3.08*70
h(x) = 0.042*x**2.5
i(x) = 1.5e-3*x**2.97
# print g(3000)
$PLOT2 max(f(x),g(x)) w l lw 2
EOF
# cat $OUT_FILE

fi