#!/bin/bash

. ./scripts/util.sh

DIR_ORIG=$1
DIR=$2
BOUNDARY_COND=$3

# echo $DIR $DIR_ORIG

ORIG_FILE=$DIR_ORIG/`get_last_tfile $DIR_ORIG`
INTERVAL=`get_attr 'Interval' $DIR`
CONTACT_MAP="$DIR/contact_map.dat"

T=0
LAST_T=$(get_last_t $DIR)

while [ $T -le $LAST_T ]; do
	CUR_FILE="$DIR/t=${T}_dev=0.res"
	echo $T `./bin/similarity $ORIG_FILE $CUR_FILE $BOUNDARY_COND` `./bin/sim_contact $CUR_FILE $CONTACT_MAP`
	let "T=T+INTERVAL"
done
