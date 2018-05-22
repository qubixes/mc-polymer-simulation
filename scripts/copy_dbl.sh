#!/bin/bash

. util.sh

DIRS=(`get_dirs -b 2 --ring -x gpupol3`)

# echo ${DIRS[*]}

for DIR in ${DIRS[*]}; do
	N=`get_attr 'Length' $DIR/../long`
	
	DIR=${DIR:0:${#DIR}-1}
	cp -r $DIR ../../doc/Files4Raoul/gpupol/N${N}/
done