#!/bin/bash

SUBDIR=(ring)
TDIRS=(`ls -d ../data/${SUBDIR}*/*/`)
NMONO="-1"
DENSITY="-1"
EQUILIBRATED="0"

while (( "$#" )); do
	case $1 in 
		-n|--nmono)
			if [ $# -lt 2 ]; then
				echo "Need number after -n/--nmono option"
				exit 1
			fi
# 			echo "Set number of monomers: $2"
			NMONO=$2
			shift ;;
		-d|--density)
			if [ $# -lt 2 ]; then
				echo "Need number after -n/--nmono option"
				exit 1
			fi
			DENSITY=$2
			shift ;;
		-e|--equilibrated)
			EQUILIBRATED=1 ;;
		*)
			echo "Unknown option: $1"
			exit 2 ;;
	esac
	shift
done

I=0
for DIR in ${TDIRS[*]}; do
	if [ ! -f "$DIR/rgyr.dat" ]; then continue; fi
	if [ "`cat "$DIR/genom.dat"`" == "" -a $EQUILIBRATED == "1" ]; then continue; fi
	LLINE=(`grep 'Length' $DIR/simulation_settings.txt`)
	LENGTHS[I]=${LLINE[2]}
	DLINE=(`grep 'Density' $DIR/simulation_settings.txt`)
	DENSITIES[I]=${DLINE[2]}
	if [ $NMONO -ne -1 -a $NMONO -ne ${LENGTHS[I]} ]; then continue; fi
	if [ "$DENSITY" != "-1" -a "${DENSITY:0:3}" != "${DENSITIES[I]:0:3}" ]; then continue; fi
	DIRS[I]=$DIR
	let "I=I+1"
done

# echo ${DIRS[*]}
# echo ${LENGTHS[*]}
# exit 0

NDIR=$I

I=0; J=0;

while [ $I -lt $NDIR ]; do
	BEST=999999;
	BESTJ=0;
	J=0;
	while [ $J -lt $NDIR ]; do
		if [ "${DIRS[J]}" != "" ]; then
			if [ ${LENGTHS[J]} -lt $BEST ]; then
				BEST=${LENGTHS[J]};
				BESTJ=$J;
			fi
		fi
		let "J=J+1"
	done
	NEW_DIRS[I]=${DIRS[BESTJ]}
	NEW_LENGTHS[I]=$BEST
	DIRS[BESTJ]=""
	let "I=I+1"
done

I=0
# NDIR=9
while [ $I -lt $NDIR ]; do
	printf "${NEW_DIRS[I]}\n"
	let "I=I+1"
done
