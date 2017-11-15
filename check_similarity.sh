#!/bin/bash

. ./scripts/util.sh

DIR_ORIG=$1
DIR=$2

# echo $DIR $DIR_ORIG

ORIG_FILE=$DIR_ORIG/`get_last_tfile $DIR_ORIG`
INTERVAL=`get_attr 'Interval' $DIR`
# let "START=$INTERVAL*450"

for T in `seq 0 $INTERVAL $(get_last_t $DIR)`; do 
	echo $T `./bin/similarity $ORIG_FILE "$DIR/t=${T}_dev=0.res"`
done
