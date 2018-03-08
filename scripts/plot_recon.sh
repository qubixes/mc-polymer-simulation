#!/bin/bash

. util.sh

function pattern_match {

local DIR=$1; shift;

while [[ $# -gt 0 ]]; do
	local PAT=$1
	if [[ ${PAT:0:1} == "~" ]]; then
		if [[ $DIR =~ ${PAT:1} ]]; then
			return 1;
		fi
	else
		if [[ ! $DIR =~ $PAT ]]; then
			return 1;
		fi
	fi
	shift
done
return 0;

}

BASE_DIR=$1
shift

ALL_DIRS=(${BASE_DIR}/rec_*/)

# echo ${ALL_DIRS[*]}

PLOT_FILES=()
NP_FILES=0


for DIR in ${ALL_DIRS[*]}; do
	if ! pattern_match `basename $DIR` $*; then 
		continue
	fi
	if [[ ! $DIR =~ _b[1,2] ]]; then
		FIRST_DIR=${DIR%*/}
		BASE_TITLE=${FIRST_DIR##*_hp}
		
		OUT_FILE="$BASE_DIR/similarity_${BASE_TITLE}.dat"
		
		cp $FIRST_DIR/similarity.dat $OUT_FILE
		
		DT=`get_last_t $FIRST_DIR`
		CURT=$DT
		
		for I in `seq 2`; do
			NEW_DIR=${FIRST_DIR}_b${I}
			if [ ! -f $NEW_DIR/simulation_settings.txt ]; then break; fi
			awk '{print $(1)+'$CURT', $(2), $(3), $(4), $(5), $(6), $(7);}' $NEW_DIR/simila* >> $OUT_FILE
			DT=`get_last_t $NEW_DIR`
			let "CURT += DT"
		done
		PLOT_FILES[NP_FILES]=$OUT_FILE
		let "NP_FILES++"
	fi
done

PLOT="plot "
PLOT_CONTACT="plot "

for FILE in ${PLOT_FILES[*]}; do 
	TITLE=`basename $FILE | sed 's/_/ /g'`
	PLOT="${PLOT} \"$FILE\" u 1:5 w l title '$TITLE', "
	PLOT_CONTACT="${PLOT_CONTACT} \"$FILE\" u 1:7 w l title '$TITLE', "
done 

PLOT=${PLOT:0:${#PLOT}-2}
PLOT_CONTACT=${PLOT_CONTACT:0:${#PLOT_CONTACT}-2}

gnuplot -persist <<EOFGNU
set grid
$PLOT

EOFGNU

gnuplot -persist <<EOFGNU
set grid
$PLOT_CONTACT
EOFGNU