#!/bin/bash

function get_paper_dir {
	echo "../denspol_paper/pics/"
}

function get_attr {
	LINE=(`grep "$1" "$2/simulation_settings.txt"`)
	echo "${LINE[2]}"
}

function contains () { 
	[[ "$1" =~ (^|[[:space:]])"$2"($|[[:space:]]) ]]; 
}

function get_title {
	DIR=$1
	STR=$2
	LENGTH=`get_attr 'Length' $DIR`
	echo "$STR" | sed 's/%n/'$LENGTH'/g'
}

function needs_update {
	SRC=$1
	DST=$2
	
	if [ ! -f $SRC ]; then return 1; fi
	if [ ! -f $DST ]; then return 0; fi 
	if [ $SRC -nt $DST ]; then return 0; fi
	return 1;
}

function get_teq {
local DIR=$1
local DENSITY=`get_attr 'Density' $DIR`
local EXEC=`get_attr 'Executable' $DIR`
local N=`get_attr 'Length' $DIR`
LINE=(`grep '^'"$EXEC ${DENSITY:0:3}" "ne_list.dat"`)
NE=${LINE[2]}
TFAC=${LINE[4]}

gnuplot << EOFGNU 2> /dev/stdout

max(x,y)=(x>y)?x:y;
tConst(z) = 128*128*2;
tRouse(z,tfac,Ne) = 1.4*z**2.2*tfac*Ne**2.2;
tRept(z,tfac,Ne) = 0.261*z**3.08*tfac*Ne**2.2;

z=$N/(1.0*$NE);
print max(tConst(z), max(tRouse(z,$TFAC,$NE), tRept(z,$TFAC,$NE)))/3.9;
EOFGNU

}

function cross_sect_files {
	F1=$1
	F2=$2
	
	if [ ! -f "$F1" -o ! -f "$F2" ]; then return; fi
	
	exec 3<$F2
	exec 4<$F1
	
# 	head $F1
	read T1 R1 <&3
	read T2 R2 <&4
	while true; do
# 		echo "$T1 $T2"
		if [ "$T1" == "$T2" ]; then 
			echo "$T2 $R1 $R2"
			read T1 R1 <&3;
		elif [ "$T1" == "" -o "$T2" == "" ]; then break;
		fi
		while [ "$T1" != "" ] && [ $T1 -lt $T2 ]; do
			read T1 R1 <&3;
		done
# 		echo "$T1 $T2"
		if [ "$T1" == "$T2" ]; then 
			echo "$T2 $R1 $R2"
			read T2 R2 <&4;
		elif [ "$T1" == "" -o "$T2" == "" ]; then break;
		fi
		
		while [ "$T2" != "" ] && [ $T2 -lt $T1 ]; do
			read T2 R2 <&4;
		done
	done
	exec 3<&-
	exec 4<&-
}

# function get_last_tfile {
# DIR=$1
# FILES=(`echo $1/t*.res`)
# MAX_T=0
# MAX_FILE="NOT_FOUND"
# for FILE in ${FILES[*]}; do
# 	if [ ! -f $FILE ]; then continue; fi
# 	TSTR=${FILE%_dev*}
# 	T=${TSTR#*t=}
# 	if [ $T -gt $MAX_T ]; then
# 		MAX_T=$T
# 		MAX_FILE=$FILE
# 	fi
# done
# echo ${MAX_FILE##*/}
# if [ $MAX_FILE == "NOT_FOUND" ]; then
# 	return 192
# else 
# 	return 0
# fi
# }

function get_last_tfile {
	DIR=$1
	LAST_T=`get_last_t $DIR`
	if [ $LAST_T == "NOT_FOUND" ]; then
		echo $LAST_T
		return 192;
	fi
	
	echo "t=${LAST_T}_dev=0.res"
	return 0
}

function get_last_savfile {
	DIR=$1
	LAST_T=`get_last_t $DIR`
	if [ $LAST_T == "NOT_FOUND" ]; then
		echo $LAST_T
		return 192;
	fi
	
	echo "sav_t${LAST_T}_dev=0.res"
	return 0
}

function get_last_t {
DIR=$1
FILES=(`echo $1/t*.res`)
MAX_T=0
MAX_FILE="NOT_FOUND"
for FILE in ${FILES[*]}; do
	if [ ! -f $FILE ]; then continue; fi
	TSTR=${FILE%_dev*}
	T=${TSTR#*t=}
	if [ $T -gt $MAX_T ]; then
		MAX_T=$T
		MAX_FILE=$FILE
	fi
done
# echo ${MAX_T}
if [ $MAX_FILE == "NOT_FOUND" ]; then
	echo "NOT_FOUND"
	return 192
else 
	echo "$MAX_T"
	return 0
fi
}


function get_dirs {

SUBDIR=(ring)
# echo ${TDIRS[*]}
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
		-l|--linear)
			SUBDIR=(lin) ;;
		-r|--ring)
			SUBDIR=(ring) ;;
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

# echo "${EXEC_LIST[*]}"
TDIRS=(../data/${SUBDIR}*/N*/ ../data/${SUBDIR}_{gpupol,denspol}*/{,double_*/} )
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