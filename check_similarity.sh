#!/bin/bash

. ./scripts/util.sh

DIR_ORIG=$1
DIR=$2

# echo $DIR $DIR_ORIG

ORIG_FILE=$DIR_ORIG/`get_last_tfile $DIR_ORIG`
INTERVAL=`get_attr 'Interval' $DIR`

T=0
LAST_T=$(get_last_t $DIR)

while [ $T -le $LAST_T ]; do
	echo $T `./bin/similarity $ORIG_FILE "$DIR/t=${T}_dev=0.res"`
	let "T=T+INTERVAL"
done
