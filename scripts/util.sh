#!/bin/bash

function get_attr {
	LINE=(`grep "$1" "$2/simulation_settings.txt"`)
	echo "${LINE[2]}"
}

function contains () { 
	[[ "$1" =~ (^|[[:space:]])"$2"($|[[:space:]]) ]]; 
}

function get_dirs {

SUBDIR=(ring)
TDIRS=(`ls -d ../data/${SUBDIR}*/*/`)
N_LIST="-1"
N_COUNT="0"
DENS_LIST="-1"
DENS_COUNT="0"
EQUILIBRATED="0"
DOUBLE_LIST="-1"
DOUBLE_COUNT="0"
EXEC_LIST="-1"
EXEC_COUNT="0"

while (( "$#" )); do
	case $1 in 
		-n|--nmono)
			if [ $# -lt 2 ]; then
				echo "Need number after -n/--nmono option"
				exit 1
			fi
			N_LIST[N_COUNT]=$2
			let "N_COUNT=N_COUNT+1"
			shift ;;
		-d|--density)
			if [ $# -lt 2 ]; then
				echo "Need number after -d/--density option"
				exit 1
			fi
			DENS_LIST[DENS_COUNT]=$2
			let "DENS_COUNT++"
			shift ;;
		-b|--double)
			if [ $# -lt 2 ]; then
				echo "Need number after -b/--double option"
				exit 1
			fi			
			DOUBLE_LIST[DOUBLE_COUNT]=$2
			let "DOUBLE_COUNT++"
			shift ;;
		-e|--equilibrated)
			EQUILIBRATED=1 ;;
		-x|--exec)
			if [ $# -lt 2 ]; then
				echo "Need exec after -x/--exec option"
				exit 1
			fi
			EXEC_LIST[EXEC_COUNT]=$2
			let "EXEC_COUNT++"
			shift ;;
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
	LENGTHS[I]=`get_attr 'Length' $DIR`
	DENSITIES[I]=`get_attr 'Density' $DIR`; DENSITIES[I]=${DENSITIES[I]:0:3};
	DOUBLES[I]=`get_attr 'Double_step' $DIR`
	EXECS[I]=`get_attr 'Executable' $DIR`
	if [[ $N_COUNT -gt 0 ]] && ! contains "${N_LIST[*]}" "${LENGTHS[I]}"; then continue; fi
	if [[ $DENS_COUNT -gt 0 ]] && ! contains "${DENS_LIST[*]}" "${DENSITIES[I]}"; then continue; fi
	if [[ $DOUBLE_COUNT -gt 0 ]] && ! contains "${DOUBLE_LIST[*]}" "${DOUBLES[I]}"; then continue; fi
	if [[ $EXEC_COUNT -gt 0 ]] && ! contains "${EXEC_LIST[*]}" "${EXECS[I]}"; then continue; fi
	DIRS[I]=$DIR
	let "I=I+1"
done

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

}