#!/bin/bash
. util.sh
TYPE=ring
DIRS=(`get_dirs $* -e`)
# LENGTHS=(`./get_lengths.sh`)
# TYPE=`./get_type.sh`
PAPER_DIR=`get_paper_dir`
# OUT_FILE="./dif_coef_$TYPE.dat"
OFILE="dif_coef_$TYPE.tmp"

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
	FUNCS=""
else
	FUNCS="f(x) = (x>0.2 && x<1)?1.5e-4*x**-1.2:1/0; g(x) = (x>10 && x<40)?1.1e-3*x**-2.4:1/0;"
	LINES='f(x) title "N^{-1.2}" lw 2, g(x) title "N^{-2.4}" lw 2'
	EXPLOT_TEX="f(x) dt 2 lw 2 lc rgb 'black' notitle, g(x) dt 2 lw 2 lc rgb 'black' notitle"
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
	RG=(`cat "$DIR/rgyr.dat"`)
	RG=${RG[0]}
	echo $D1 $D2 $RG $NE_FAC $R_FAC $T_FAC >> $FILE
	let "I=I+1"
done

PLOT="plot [:][:0.5]"
PLOT_TEX="plot [:][:0.5]"
PLOT2="plot [1:50]"

I=0
for FILE in ${FILES[*]}; do
	BTITLE="{/Symbol r}=${DENSITIES[I]}, exec=${EXECS[I]}"
# 	PLOT="$PLOT \"$FILE\" u (\$1/\$7):(\$3*\$9*\$7**1.2) w l title \"$BTITLE\" lw 2, "
	PLOT="$PLOT \"$FILE\" u (\$1/\$7):(\$3) w l title \"$BTITLE\" lw 2, "
	PLOT_TEX="$PLOT_TEX \"$FILE\" u (\$1/\$7):(\$3) w l title \"MC Data\" lw 2, "
	PLOT2="$PLOT2 \"$FILE\" u (\$1/\$7):(3.9*\$6/\$3/6.0/\$9/\$7**2.2) w l title \"$BTITLE (Diffuse RG)\" lw 2, "
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
set xlabel "N/N_e"
$FUNCS
$PLOT $LINES
EOF

gnuplot << EOF
set term aqua enhanced font 'Helvetica, 18'
set log x
set log y
set format y "10^{%T}"
set ylabel "t"
set xlabel "N/N_e"
set key bottom right
alpha=5
max(x,y) = (x>y)?x:y;
f(x) = 1.4*x**2.21
g(x) = 0.261*x**3.08
h(x) = 0.042*x**2.5
i(x) = 1.5e-3*x**2.97
# print g(3000)
$PLOT2 max(f(x),g(x)) w l lw 2
EOF
# cat $OUT_FILE

fi
# 
# gnuplot << EOFGNU
# set term epslatex color standalone colortext 14
# set log x
# set log y
# set format y "\$10^{%T}\$"
# set label '{\\tiny \$ N^{-1.2} \$}' at 0.6, 1.1e-3  
# set label '{\\tiny \$ N^{-2.4} \$}' at 25,3e-6
# set ylabel '\$D\$'
# set xlabel '\$N/N_e\$'
# $FUNCS
# set out "$OFILE.tex"
# $PLOT_TEX $EXPLOT_TEX
# EOFGNU
# make $OFILE.eps
# cp $OFILE.eps $PAPER_DIR

# exit 0
