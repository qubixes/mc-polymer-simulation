#!/bin/bash

DIRS=(`ls -d ../data/*/`)

for DIR in ${DIRS[*]}; do
	if [ ! -f ${DIR}/simulation_settings.txt ]; then 
		echo "${DIR}"
	fi
done
