#!/bin/bash

. util.sh

DIRS=`./get_dirs.sh -x gpupol3`

for DIR in ${DIRS[*]}; do
	FILE=`get_last_tfile $DIR/long`
	FILE="$DIR/long/$FILE"
# 	echo $FILE
	
	DEST_FILE="/scratch/schram/gpupol_data/`echo $FILE | sed 's/..\/data//g'`"
# 	echo $DEST_FILE
	scp $FILE k40.cbp.ens-lyon.fr:$DEST_FILE
done
