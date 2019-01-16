#!/bin/bash

. util.sh

DIRS=(`get_dirs -b 0 --ring -x gpupol3`)

# echo ${DIRS[*]}

CURDIR=`pwd`
cd ../data

rsync -vaz --update --include '*.dat' --include '*/' --exclude '*' k40.cbp.ens-lyon.fr:/scratch/FromBackup/schram/gpupol_data/ .

cd $CURDIR

for DIR in ${DIRS[*]}; do
	for SUB in long double_1 double_2; do
# 		echo $SUB
		N=`get_attr 'Length' $DIR/long`
		SRC_DIR=$DIR/$SUB
		DEST_DIR=../../doc/Files4Raoul/gpupol/N${N}/$SUB
# 		ls $SRC_DIR/*.dat
# 		ls $DEST_DIR/*.dat
		cp $SRC_DIR/*.dat $DEST_DIR
	done
done
# 
# 	
# 	DIR=${DIR:0:${#DIR}-1}
# 	cp -r $DIR ../../doc/Files4Raoul/gpupol/N${N}/
# done
