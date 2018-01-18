#!/bin/bash

. ./scripts/util.sh

BIN_DIR=./bin

function abs_path {
	X=`cd $(dirname $1); pwd -P`
	echo $X
}

if [ $# -lt 2 ]; then
	echo "Need at least two arguments: directory and HP strength"
	exit 192;
fi

TOPO="FALSE"
POL_TYPE="ring"
NPC_SAMPLES=""
EE_TOPO_FILE=./denspol/ee_topo_comp.dat

DIR=$1
HP_STRENGTH=$2

SEED=12938173
TIME=2e5
INTERVAL=1e3
FIRST_TIME=2e5
FIRST_INTERVAL=1e4

BEND_ENERGY=0.3
SHUFFLE=1
HP_STRENGTH_2=`echo "$HP_STRENGTH+0.1" | bc -l`
HP_STRENGTH_3=`echo "$HP_STRENGTH+0.2" | bc -l`


if [ $# -ge 3 -a "$3" == "TOPO" ]; then
	TOPO="TRUE"
fi

if [ $# -ge 4 -a "$4" == "lin" ]; then
	POL_TYPE="lin"
fi

if [ $# -ge 5 ]; then
	NPC_SAMPLES=$5
fi

if [ $# -ge 6 ]; then
	SEED=$6
fi

if [ $TOPO == "TRUE" ]; then
	EXEC=$BIN_DIR/denspol_${POL_TYPE}_hp_topo
else
	EXEC=$BIN_DIR/denspol_${POL_TYPE}_hp
fi

if [ $POL_TYPE == "lin" ]; then
	FIRST_EXEC=$BIN_DIR/denspol_${POL_TYPE}_hp_topo
else
	FIRST_EXEC=$EXEC
fi


if [ ! -f "$DIR/simulation_settings.txt" ]; then
	echo "Directory $DIR is not a valid directory (missing simulation_settings.txt)"
	exit 193;
fi

make all install &> /dev/null || { echo "Error in compilation.\n"; exit 192; }

SCRIPT_DIR=`pwd`
cd $DIR/..
BASE_DIR=`pwd`
cd $SCRIPT_DIR
LENGTH=`get_attr 'Length' $DIR`
L=`get_attr 'Latsize' $DIR`
DENSITY=`get_attr 'Density' $DIR`
NPOL=`get_attr 'Npol' $DIR`


DENSITY=`echo "$DENSITY/2.0" | bc -l`
# echo $DENSITY; exit 0
# echo "get_attr 'Npol' $BASE_DIR"
# echo "NPOL=$NPOL"
# exit 0
NSTEP=0
CUR_LENGTH=$LENGTH

CUR_LENGTH=7;
MIN_L=4;
BEST_L=-1;
BEST_DL=100;
BEST_LENGTH=-1;
MAX_VALID_LENGTH=9999999
while true; do
	LVAL=`echo "$DENSITY $NPOL $CUR_LENGTH" | awk 'function pow(x,y) {return exp(y*log(x))}; {print pow($(2)*$(3)/$(1), 1./3.);}'`
	LINT=`printf "%.0f" $LVAL`
	DLVAL=`echo "$LINT $LVAL" | awk 'function abs(v) {return v < 0 ? -v : v}; {print abs($(1)-$(2))}'`
	echo $LVAL $LINT $DLVAL $CUR_LENGTH $DENSITY $NPOL
	if [ `echo "$LVAL >= $MIN_L" | bc ` == "1" ]; then
		if [ $BEST_LENGTH == "-1" -o `echo "$DLVAL*$CUR_LENGTH < $BEST_DL*$BEST_LENGTH" | bc` == "1" ]; then
			if [ `echo "$NPOL < $LINT*$LINT*$LINT" | bc` == "1" ]; then
				if [ $MAX_VALID_LENGTH == "9999999" ]; then
					let "MAX_VALID_LENGTH=2*CUR_LENGTH"
				fi
				BEST_DL=$DLVAL
				BEST_L=$LINT
				BEST_LENGTH=$CUR_LENGTH
			fi
		fi
	fi
	
	if [ $CUR_LENGTH -gt $MAX_VALID_LENGTH -a $BEST_DL != "100" ]; then break; fi
	if [ `echo "$BEST_DL<0.1" | bc` == "1" ]; then break; fi	
	let "CUR_LENGTH++"
done 

# echo $BEST_DL $BEST_L $BEST_LENGTH 
# exit 01
# while [ `echo "$CUR_LENGTH>28" | bc -l` == "1" ]; do
# 	CUR_LENGTH=`echo "scale=12; $CUR_LENGTH/8" | bc -l`
# 	let  "NSTEP++"
# done

START_LENGTH=$BEST_LENGTH
START_L=$BEST_L
DENSITY=`echo "$NPOL $START_L $START_LENGTH" | awk '{print $(1)*$(3)/($(2)*$(2)*$(2));}'`

# echo "$START_L $START_LENGTH $DENSITY $NPOL"
# exit 0

# HP_STRENGTH=0.5

ISTEP=0

BASE_DEST_DIR="$BASE_DIR/rec_l${START_LENGTH}_hp${HP_STRENGTH}_t${TIME}_${POL_TYPE}_s${SEED}"

if [ $TOPO == "TRUE" ]; then
	BASE_DEST_DIR="${BASE_DEST_DIR}_TOPO"
fi

if [ $# -ge 5 ]; then
	BASE_DEST_DIR="${BASE_DEST_DIR}_ns${NPC_SAMPLES}"
fi

FIRST_DIR=$BASE_DEST_DIR
mkdir -p $FIRST_DIR
SRC_FILE=$DIR/`get_last_tfile $DIR`
$BIN_DIR/contact_map $SRC_FILE $FIRST_DIR/contact_map.dat $NPC_SAMPLES || exit $?
# echo "$BIN_DIR/contact_map $SRC_FILE $FIRST_DIR/contact_map.dat $NPC_SAMPLES"; exit 0
LAST_TFILE=`get_last_tfile $FIRST_DIR`
LAST_SAVFILE=`get_last_savfile $FIRST_DIR`

# echo "$LAST_TFILE"; exit 192;

if [ $LAST_TFILE == "NOT_FOUND" -o $LAST_TFILE == "t=0_dev=0.res" ]; then
# 	echo "$EXEC $SEED $FIRST_DIR $DENSITY $TIME $INTERVAL $START_LENGTH $START_L $EE_TOPO_FILE $BEND_ENERGY 0 $SHUFFLE $FIRST_DIR/contact_map.dat $HP_STRENGTH"
	$FIRST_EXEC $SEED $FIRST_DIR $DENSITY $FIRST_TIME $FIRST_INTERVAL $START_LENGTH $START_L $EE_TOPO_FILE $BEND_ENERGY 0 $SHUFFLE $FIRST_DIR/contact_map.dat $HP_STRENGTH || exit $?
# 	$EXEC $SEED $FIRST_DIR $DENSITY $FIRST_TIME $FIRST_INTERVAL $START_LENGTH $START_L $EE_TOPO_FILE $BEND_ENERGY 0 0 $FIRST_DIR/contact_map.dat $HP_STRENGTH_2 || exit $?
# 	$EXEC $SEED $FIRST_DIR $DENSITY $FIRST_TIME $FIRST_INTERVAL $START_LENGTH $START_L $EE_TOPO_FILE $BEND_ENERGY 0 0 $FIRST_DIR/contact_map.dat $HP_STRENGTH_3 || exit $?
# 	exit 0
	./check_similarity.sh $DIR $FIRST_DIR > "$FIRST_DIR/similarity.dat"
	$BIN_DIR/contact_map $FIRST_DIR/`get_last_tfile $FIRST_DIR` "$FIRST_DIR/contact_map_new.dat"
fi

DBL_STEP=1
CUR_LENGTH=$START_LENGTH
CUR_L=$START_L

while [ $CUR_LENGTH -lt $LENGTH ]; do
	SECOND_DIR="${BASE_DEST_DIR}_b${DBL_STEP}"
	mkdir -p ${SECOND_DIR}
	
	CUR_TFILE=`get_last_tfile $SECOND_DIR`
	let "CUR_LENGTH*=8"
	let "CUR_L*=2"
	
	LAST_TFILE=$FIRST_DIR/`get_last_tfile $FIRST_DIR`
	LAST_SAVFILE=$FIRST_DIR/`get_last_savfile $FIRST_DIR`
	
	if [ $CUR_TFILE == "NOT_FOUND" ]; then
		$BIN_DIR/denspol_scaleup $LAST_TFILE $LAST_SAVFILE "$SECOND_DIR/t=0_dev=0.res" "$SECOND_DIR/sav_t0_dev=0.res" ./denspol/straight_topo.dat || exit $?
		cp $FIRST_DIR/contact_map.dat $SECOND_DIR
		$EXEC $SEED $SECOND_DIR $DENSITY $TIME $INTERVAL $CUR_LENGTH $CUR_L $EE_TOPO_FILE $BEND_ENERGY 1 0 $SECOND_DIR/contact_map.dat $HP_STRENGTH
		$BIN_DIR/contact_map $SECOND_DIR/`get_last_tfile $SECOND_DIR` $SECOND_DIR/contact_map_new.dat
	fi
	
	SIM_FILE="$SECOND_DIR/similarity_hp${HP_STRENGTH}.dat"
	
	if `needs_update $SECOND_DIR/$(get_last_tfile $SECOND_DIR) $SIM_FILE`; then
		./check_similarity.sh $DIR $SECOND_DIR > $SIM_FILE
	fi
	
	FIRST_DIR=$SECOND_DIR
	let "DBL_STEP++"
done

# echo "SHUFFLE=$SHUFFLE"
# exit 192

# echo "$BIN_DIR/denspol_ring_hp $SEED $FIRST_DIR $DENSITY $TIME $INTERVAL $START_LENGTH $START_L $EE_TOPO_FILE $BEND_ENERGY 0 $SHUFFLE $FIRST_DIR/contact_map.dat $HP_STRENGTH"
# exit 0


# ./update.sh

# for i in `seq 1 144`; do 
# ./bin/mobility $FIRST_DIR/ptl/pol\=${i}_dev\=0.res | awk '{sum+=$(2); n++;}END{print sum/n}'; done

# while [ $ISTEP -lt $NSTEP ]; do
# 	
# 	
# 	let "ISTEP++"
# done

# echo $LENGTH $CUR_LENGTH $NSTEP $START_LENGTH

#plot "pc_test.dat" u ($1*100+$2):($3*100+$4):5 w circles