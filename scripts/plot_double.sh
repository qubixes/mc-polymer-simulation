#!/bin/bash
T=3500000

. ./util.sh
DIRS=(`get_dirs $* -e`) || { echo ${DIRS[*]}; exit $?; }

DATA_DIR="../data"
RGYR_FILE="$DATA_DIR/rgyr`echo "$*" | sed 's/ //g' | sed 's/-/_/g'`.dat"
RGYR_DBL_FILE="$DATA_DIR/rgyr_dbl`echo "$*" | sed 's/ //g' | sed 's/-/_/g'`.dat"
rm -f $RGYR_FILE $RGYR_DBL_FILE

I=0


for DIR in ${DIRS[*]}; do
	RGYR=(`cat $DIR/rgyr.dat`)
		
	N=`get_attr 'Length' $DIR`
# 	NP=`ls $DIR/ptl | wc -w`
	echo "$N $RGYR" >> $RGYR_FILE
	TFILE="$DIR/double_1/rgyr_time.dat"
	if [ -f "$TFILE" ]; then
		LINE=(`grep ^$T "$TFILE"`)
		RGYR_DBL=${LINE[1]}
		let "N_DBL=N*8"
		echo "$N_DBL $RGYR_DBL" >> $RGYR_DBL_FILE
	fi
done

# cat $RGYR_DBL_FILE

gnuplot << EOF
set term aqua enhanced font 'Helvetica, 22'
set log x
set log y
set xlabel "N"
set ylabel "r_g^2/N^{2/3}"
set key top right
plot "$RGYR_FILE" u 1:(\$2/\$1**(2./3.)) w l lt 1 lw 2 title "Brute force", "$RGYR_DBL_FILE" u 1: (\$2/\$1**(2./3.)) title "Hierarchical building with {/Symbol t}=$T"
EOF