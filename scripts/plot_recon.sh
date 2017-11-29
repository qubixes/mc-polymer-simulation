#!/bin/bash

. util.sh

BASE_DIR=$1

ALL_DIRS=(${BASE_DIR}/rec_*/)

# echo ${ALL_DIRS[*]}

PLOT_FILES=()
NP_FILES=0

for DIR in ${ALL_DIRS[*]}; do
# 	echo $DIR
	if [[ ! $DIR =~ _b[1,2] ]]; then
		if [[ $# -gt 1 ]] && [[ ! $DIR =~ $2 ]]; then continue; fi
		FIRST_DIR=${DIR%*/}
		BASE_TITLE=${FIRST_DIR##*_hp}
		
		OUT_FILE="$BASE_DIR/similarity_${BASE_TITLE}.dat"
		
		cp $FIRST_DIR/similarity_* $OUT_FILE
		
		DT=`get_last_t $FIRST_DIR`
		CURT=$DT
		
		for I in `seq 2`; do
			NEW_DIR=${FIRST_DIR}_b${I}
			if [ ! -f $NEW_DIR/simulation_settings.txt ]; then break; fi
			awk '{print $(1)+'$CURT', $(2), $(3), $(4), $(5);}' $NEW_DIR/simila* >> $OUT_FILE
			let "CURT += DT"
		done
		PLOT_FILES[NP_FILES]=$OUT_FILE
		let "NP_FILES++"
	fi
done

PLOT="plot "

for FILE in ${PLOT_FILES[*]}; do 
	PLOT="${PLOT} \"$FILE\" u 1:5 w l title '`basename $FILE | sed 's/_/ /g'`', "
done 

PLOT=${PLOT:0:${#PLOT}-2}


gnuplot -persist <<EOFGNU

$PLOT

EOFGNU