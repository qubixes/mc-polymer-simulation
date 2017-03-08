#!/bin/bash

NPROC=4

if [ $# -gt 0 ]; then
	NPROC=$1
fi

if [ $# -gt 1 ]; then
	BDIR=$2
else
	BDIR="./data"
fi

DIRS=($BDIR/*/*/);

for DIR in ${DIRS[*]}; do 
	./bin/create_ptl $DIR || exit $?
done

parallel -j $NPROC ./bin/lowm_modes ::: ${DIRS[*]}

for DIR in ${DIRS[*]}; do 
	if [ -f $DIR/shearmod.dat -a ! $DIR/shearmod_avg.dat -nt $DIR/shearmod.dat ]; then
		./bin/avg_data $DIR/shearmod.dat > $DIR/shearmod_avg.dat || echo "Failed: $DIR"
	fi
	if [ -f $DIR/pc.dat -a ! $DIR/pc_point.dat -nt $DIR/pc.dat ]; then
		./scripts/plot_2d_pc.sh $DIR/pc.dat || echo "Failed: $DIR"
	fi
	if [ -f $DIR/pc_point.dat -a ! $DIR/dr.dat -nt $DIR/pc_point.dat ]; then
		echo "./bin/localization $DIR/pc_point.dat > $DIR/dr.dat"
		./bin/localization $DIR/pc_point.dat > $DIR/dr.dat || echo "Failed: $DIR"
	fi
	
	if [ -f $DIR/pc_point.dat -a ! $DIR/pc_avg.dat -nt $DIR/pc_point.dat ]; then
		./scripts/get_pc_avg.sh $DIR/pc_point.dat > $DIR/pc_avg.dat || echo "Failed: $DIR, pc avging"
	fi
done

for DIR in ${DIRS[*]}; do 
	./bin/create_cms $DIR || exit $?
done

parallel -j $NPROC ./bin/cms_cor ::: ${DIRS[*]}

BASE_DIRS=(`echo $BDIR/*gpupol*`);

MERGE_FILES=("cmsdif.dat" "emdif.dat" "mmdif.dat" "smdif.dat" )
LONG_FILES=("slrat.dat" "rgyr.dat" "genom.dat" "ucor.dat" "simulation_settings.txt" "pc_avg.dat" "rouse_stat.dat" "rgyr_time.dat")
ROUSE_FILES=("ree.dat" "rouse_dyn.dat")
ALL_FILES=(${MERGE_FILES[*]} ${LONG_FILES[*]} ${ROUSE_FILES[*]})

for DIR in ${BASE_DIRS[*]}; do
	for FILE in ${ALL_FILES[*]}; do
		if [ -d $DIR/long ]; then
			cp $DIR/long/$FILE $DIR
		fi
	done
done

for DIR in ${BASE_DIRS[*]}; do
	SH_DIR=($DIR/short*)
	LN_DIR=($DIR/long*)
# 	echo $SH_DIR $LN_DIR
	if [ ! -d $SH_DIR -o ! -d $LN_DIR ]; then continue; fi
	
	for FILE in ${MERGE_FILES[*]}; do
		./bin/merge_tseries {$SH_DIR,$LN_DIR}/$FILE > $DIR/$FILE || echo "Failed $FILE"
	done
	
	for FILE in ${ROUSE_FILES[*]}; do
		./bin/merge_tseries {$SH_DIR,$LN_DIR}/$FILE 1 > $DIR/$FILE || echo "Failed $FILE"
	done
done