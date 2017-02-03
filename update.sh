#!/bin/bash

NPROC=4

if [ $# -gt 0 ]; then
	NPROC=$1
fi

DIRS=(./data/*/*/);

for DIR in ${DIRS[*]}; do 
	./bin/create_ptl $DIR || exit $?
done

parallel -j $NPROC ./bin/lowm_modes ::: ${DIRS[*]}
