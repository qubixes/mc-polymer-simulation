#!/bin/bash

. ./scripts/util.sh

DIR_ORIG=$1
DIR=$2
BOUNDARY_COND=$3

# echo $DIR $DIR_ORIG

ORIG_FILES=$DIR_ORIG/`get_last_tfile $DIR_ORIG`
if [ `get_last_tfile $DIR_ORIG` == "NOT_FOUND" ]; then
	SIM_EXEC=./bin/similarity_pascal
	BIN_FILE=$DIR_ORIG/bins*
	CONFIG_FILE=$DIR_ORIG/config*
	ORIG_FILES="$BIN_FILE $CONFIG_FILE"
else
	SIM_EXEC=./bin/similarity
fi
INTERVAL=`get_attr 'Interval' $DIR`
CONTACT_MAP="$DIR/contact_map.dat"

T=0
LAST_T=$(get_last_t $DIR)

while [ $T -le $LAST_T ]; do
	CUR_FILE="$DIR/t=${T}_dev=0.res"
	echo $T `$SIM_EXEC $ORIG_FILES $CUR_FILE $BOUNDARY_COND` `./bin/sim_contact $CUR_FILE $CONTACT_MAP`
	let "T=T+INTERVAL"
done
