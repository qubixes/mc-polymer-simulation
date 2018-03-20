#!/bin/bash

. ./scripts/util.sh

BIN_DIR=./bin

function abs_path {
	X=`cd $(dirname $1); pwd -P`
	echo $X
}

TOPO=""
NPC_SAMPLES=""
EE_TOPO_FILE=./denspol/ee_topo_comp.dat

BOUNDARY_COND="periodic"
LATTICE_SHAPE=""
HP_STRENGTH=0.35

SEED=12938173
TIME=1e5
INTERVAL=1e3
DENSITY=7.2
START_L=6
POL_TYPE="original"

BEND_ENERGY=0.3
SHUFFLE=1

DIR="$1"
shift

while (( "$#" )); do
	case $1 in 
		-t|--time)
			if [ $# -lt 2 ]; then
				echo "Need amount of time after -t/--time option."
				exit 1
			fi
			TIME=$2
			shift;;
		-s|--seed)
			if [ $# -lt 2 ]; then
				echo "Need seed after -s/--seed option"
				exit 1
			fi
			SEED=$2
			shift;;
		-i|--interval)
			if [ $# -lt 2 ]; then
				echo "Need interval after -i/--interval option"
				exit 1
			fi
			INTERVAL=$2
			shift;;
		--phantom|--topo)
			TOPO="--topo" ;;
		-n|--samples)
			if [ $# -lt 2 ]; then
				echo "Need number of samples after -n/--samples option"
				exit 1
			fi
			NPC_SAMPLES=$2
			shift;;
		-e|--hpstrength)
			if [ $# -lt 2 ]; then
				echo "Need attraction strength after -e/--hpstrength option"
				exit 1
			fi
			HP_STRENGTH=$2
			shift;;
		--hardsphere)
			LATTICE_SHAPE="--lattice-sphere"
			;;
		--linear)
			FORCE_LINEAR="1"
			POL_TYPE="lin"
			;;
		--ring)
			FORCE_RING="1"
			POL_TYPE="ring"
			;;
		*)
			echo "Error: unknown option $1"
			exit 1 ;;
	esac
	shift
done

HP_STRENGTH_2=`echo "$HP_STRENGTH+0.05" | bc -l`
HP_STRENGTH_3=`echo "$HP_STRENGTH+0.10" | bc -l`

FIRST_TIME=`echo "print $TIME*10" | gnuplot 2> /dev/stdout`
FIRST_INTERVAL=`echo "print $INTERVAL*10" | gnuplot 2> /dev/stdout`

EXEC=$BIN_DIR/denspol_hp


if [ ! -f "$DIR/simulation_settings.txt" ]; then
	echo "Directory $DIR is not a valid directory (missing simulation_settings.txt)"
	exit 193;
fi

#make all install &> /dev/null || { echo "Error in compilation."; exit 192; }

SCRIPT_DIR=`pwd`
cd $DIR/..
BASE_DIR=`pwd`
cd $SCRIPT_DIR
LENGTH=`get_attr 'Length' $DIR`
L=`get_attr 'Latsize' $DIR`
DENSITY=`get_attr 'Density' $DIR`
NPOL=`get_attr 'Npol' $DIR`


L=$START_L;


BASE_DEST_DIR="$BASE_DIR/reconstruction/rec_hp${HP_STRENGTH}_${POL_TYPE}_s${SEED}`echo $TOPO | sed 's/--/_/'`"

if [ "$NPC_SAMPLES" != "" ]; then
	BASE_DEST_DIR="${BASE_DEST_DIR}_ns${NPC_SAMPLES}"
fi

FIRST_DIR=$BASE_DEST_DIR
CONTACT_FILE="$FIRST_DIR/contact_map.dat"
T_CONTACT_FILE="$FIRST_DIR/temp_contact_map.dat"
mkdir -p $FIRST_DIR
SRC_FILE=$DIR/`get_last_tfile $DIR`
$BIN_DIR/contact_map $SRC_FILE $T_CONTACT_FILE $SEED $NPC_SAMPLES || exit $?

if [ "$FORCE_LINEAR" == "1" ]; then
	cat $T_CONTACT_FILE | sed 's/ring/lin/g' > $CONTACT_FILE
elif [ "$FORCE_RING" == "1" ]; then
	cat $T_CONTACT_FILE | sed 's/lin/ring/g' > $CONTACT_FILE
else
	cp $T_CONTACT_FILE $CONTACT_FILE
fi

rm $T_CONTACT_FILE

if [ "$LATTICE_SHAPE" == "--lattice-sphere" ]; then
	BOUNDARY_COND="static"
	EXTRA_SCALEUP_OPTS="sphere $CONTACT_FILE $DENSITY"
fi

LAST_TFILE=`get_last_tfile $FIRST_DIR`
LAST_SAVFILE=`get_last_savfile $FIRST_DIR`

COMMON_OPTS="-x denspol -s $SEED --double 1 --hp $HP_STRENGTH -b $BEND_ENERGY -d $DENSITY $LATTICE_SHAPE"

if [ $LAST_TFILE == "NOT_FOUND" -o $LAST_TFILE == "t=0_dev=0.res" ]; then
	./do_run.sh $COMMON_OPTS -t $FIRST_TIME -i $FIRST_INTERVAL -m $FIRST_DIR/contact_map.dat --outdir $FIRST_DIR -g $L || exit $?
	
	./check_similarity.sh $DIR $FIRST_DIR $BOUNDARY_COND > "$FIRST_DIR/similarity.dat"
	$BIN_DIR/contact_map $FIRST_DIR/`get_last_tfile $FIRST_DIR` "$FIRST_DIR/contact_map_new.dat"
fi

DBL_STEP=1
MAX_L=20


while [ $L -lt $MAX_L ]; do
	SECOND_DIR="${BASE_DEST_DIR}_b${DBL_STEP}"
	mkdir -p ${SECOND_DIR}
	
	CUR_TFILE=`get_last_tfile $SECOND_DIR`
	let "L*=2"
	
	LAST_TFILE=$FIRST_DIR/`get_last_tfile $FIRST_DIR`
	LAST_SAVFILE=$FIRST_DIR/`get_last_savfile $FIRST_DIR`
	
	if [ $CUR_TFILE == "NOT_FOUND" ]; then
		$BIN_DIR/denspol_scaleup $LAST_TFILE $LAST_SAVFILE "$SECOND_DIR/t=0_dev=0.res" "$SECOND_DIR/sav_t0_dev=0.res" ./denspol/straight_topo.dat $EXTRA_SCALEUP_OPTS 
		cp $FIRST_DIR/contact_map.dat $SECOND_DIR
		
		./do_run.sh $COMMON_OPTS -t $TIME -i $INTERVAL -m $SECOND_DIR/contact_map.dat --outdir $SECOND_DIR -g $L || exit $?
		
		$BIN_DIR/contact_map $SECOND_DIR/`get_last_tfile $SECOND_DIR` $SECOND_DIR/contact_map_new.dat
	fi
	
	SIM_FILE="$SECOND_DIR/similarity_hp${HP_STRENGTH}.dat"
	
	if `needs_update $SECOND_DIR/$(get_last_tfile $SECOND_DIR) $SIM_FILE`; then
		./check_similarity.sh $DIR $SECOND_DIR $BOUNDARY_COND > $SIM_FILE
	fi
	
	FIRST_DIR=$SECOND_DIR
	let "DBL_STEP++"
done
