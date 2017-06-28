#!/bin/bash

. $HOME/.bashrc
# shopt -s extglob

USAGE="usage: do_run.sh [OPTIONS]\nAvailable Options:\n\n-s/--seed   Seed\n-n/--nmono  Polymer length\n-t/--time   Simulation time\n-h/--help   Print usage of do_run.sh\n-o/--opencl Use openCL\n-c/--cpu    Use CPU\n-u/--cuda   Use CUDA\n"


function abs_path {
	X=`cd $(dirname $1); pwd -P`
	echo $X
}

function last_modification {
	if [ `uname` == "Darwin" ]; then
		stat -f "%m" $1
	else
		date -r $1 +%s
	fi
}

function build_args {
	OPTS=""
	for j in `seq $#`; do
		let "jm=$j-1"
		OPTS="$OPTS -${OPTNAME[jm]} $1"; shift
	done
}

function is_num { 
	NUM=`echo $1 | sed 's/[0-9]//g;'`
	if [ $NUM ]; then
		return 0
	else
		return 1
	fi
}

function is_val_func { 
	NUM=`echo $1 | sed 's/[0-9L\+\*\^]//g'`
	if [ $NUM ]; then
		return 0
	else
		return 1
	fi
}

function remove_schar {
	echo $1 | sed 's/[\+\*\^\/]//g'
}

CUR_DIR=`abs_path $0`
DATA_DIR="$CUR_DIR/data"
BIN_DIR="$CUR_DIR/bin"
SCRIPT_DIR="$CUR_DIR/scripts"

. $SCRIPT_DIR/util.sh

START_TIME=`date`
SEED="12396143"
TIME="1e4"
SIM_TYPE="ring"
DENSITY="1.2"
INTERVAL="1e2"
FAST_EQ="0"
NMONO="345"
MODEL="sl_equal"
L="96"
DBL_STEP="0"
SHORT="0"
B_EXEC="gpupol3"
BEND_ENERGY="0.0"

while (( "$#" )); do
	case $1 in 
		-x|--exec)
			if [ $# -lt 2 ]; then 
				echo "Need exec after -x/--exec option: $1"
				exit 1
			fi
			B_EXEC=$2
			shift;;
		-r|--dir)
			if [ $# -lt 2 ]; then
				echo "need directory after -r/--dir option"
				exit 2
			fi
			echo "Set alternative dir: $2"
			DATA_DIR=$2
			shift ;;
		-b|--bend)
			if [ $# -lt 2 ]; then
				echo "Need an energy after -b/--bend option: $1"
				exit 1
			fi
			echo "Set bending energy: $2"
			BEND_ENERGY=$2
			shift;;
		-u|--double)
			if [ $# -lt 2 ]; then
				echo "Need number after double option: $1"
				exit 1
			fi
			echo "Set number of doubling steps: $2"
			DBL_STEP=$2
			shift;;
		-g|--grid)
			if [ $# -lt 2 ]; then
				echo "Need number after grid option: $1"
				exit 1
			fi
			echo "Set initial lattice size: $2"
			L=$2
			shift;;
		-s|--seed)
			if [ $# -lt 2 ]; then
				echo "Need number after seed option: $1"
				exit 1
			fi
			echo "Set alternate seed: $2"
			SEED=$2
			shift ;;
		-t|--time)
			if [ $# -lt 2 ]; then
				echo "need number after -t/--time option"
				exit 2
			fi
			echo "Set simulation time: $2"
			TIME=$2
			shift ;;
		-n|--nmono) 
			if [ $# -lt 2 ]; then
				echo "Need number after -n/--nmono option"
				exit 1
			fi
			echo "Set number of monomers: $2"
			NMONO=$2
			shift ;;
		-i|--interval)
			if [ $# -lt 2 ]; then
				echo "Need number after -i/--interval option."
				exit 3
			fi
			echo "Set interval: $2"
			INTERVAL=$2
			shift;;
		-h|--help)
			printf "$USAGE"
			exit 0 ;;
		-r|--ring)
			echo "Doing simulation with only ring polymers"
			SIM_TYPE="ring" ;;
		-d|--density)
			if [ $# -lt 2 ]; then
				echo "Need number after -d/--density option."
				exit 3
			fi
			echo "Set density: $2"
			DENSITY=$2
			shift;;
		-l|--linear)
			echo "Doing simulation, including linear polymers"
			SIM_TYPE="lin" ;;
		-m|--model)
			if [ $# -lt 2 ]; then
				echo "Need model after -m/--model option."
				exit 3
			fi
			echo "Using \"$2\" simulation model"
			MODEL=$2
			shift;;
		-f|--fast-eq)
			if [ $# -lt 2 ]; then
				echo "Need interval after -f/--fast-eq option. It must be smaller than the given interval."
				exit 3
			fi
			echo "Using fast equilibration, using redistribution of stored length."
			FAST_EQ=$2
			shift ;;
		--short)
			echo "Doing short time simulations"
			SHORT="1" ;;
		*)
			echo "Unknown option: $1"
			exit 3
	esac
	shift
done

if is_num $NMONO ; then
	echo "$NMONO is not a valid number"
	exit 192
fi

if is_num $SEED ; then
	echo "$SEED is not a valid number"
	exit 192
fi

if is_num $L ; then
	echo "$L is not a valid number"
	exit 192
fi

if [ $B_EXEC == "efpol" ]; then 
	EXEC="$BIN_DIR/efpol_$SIM_TYPE"
	DIR="$DATA_DIR/${SIM_TYPE}_efpol_l${NMONO}_g${L}_b${DBL_STEP}_s${SEED}_d${DENSITY}_t${TIME}"
	BASE_DIR=$DIR
elif [ $B_EXEC == "denspol" ]; then
	EXEC="$BIN_DIR/denspol"
	BASE_DIR="$DATA_DIR/ring_denspol_l${NMONO}_g${L}_s${SEED}_d${DENSITY}_t${TIME}_b${BEND_ENERGY}"
elif [ $B_EXEC == "gpupol2" -o $B_EXEC == "gpupol3" ]; then
	EXEC="$BIN_DIR/${B_EXEC}_cuda_$SIM_TYPE"
	
	if [ $FAST_EQ != "0" ]; then
	BASE_DIR="$DATA_DIR/${SIM_TYPE}_${B_EXEC}_l${NMONO}_g${L}_s${SEED}_d${DENSITY}_f${FAST_EQ}"
	else 
		BASE_DIR="$DATA_DIR/${SIM_TYPE}_${B_EXEC}_l${NMONO}_g${L}_s${SEED}_d${DENSITY}"
	fi
# 	DIR="$BASE_DIR/long"
# 	EXEC_LINE="$EXEC $NMONO $TIME $SEED $DIR $DENSITY 0 $INTERVAL $L $L $L"
else 
	echo "Not a valid algorithm: $B_EXEC"
	exit 192
fi

if [ $SHORT == "1" -a $B_EXEC != "efpol" ]; then
	if [ ! -d $BASE_DIR/long ]; then 
		echo "Error: running short simulation without the long [$BASE_DIR]"; exit 192; 
	fi
	
	DIR=$BASE_DIR/short_$INTERVAL
	if [ -d $DIR ]; then
		echo "Directory $DIR is already there, remove it first."
		exit 192;
	fi
	mkdir -p $DIR || { echo "Error creating directory $DIR"; exit 192; } 
	IN_FILE=`get_last_tfile $BASE_DIR/long/` || { echo "Error finding file in $BASE_DIR/long"; exit $?; }
	cp "$BASE_DIR/long/$IN_FILE" $DIR/t=0_dev=0.res || exit $?
	echo "Using file:  $IN_FILE"
elif [ $DBL_STEP -gt 0 -a $B_EXEC != "efpol" ]; then
	if [ ! -d $BASE_DIR/long ]; then 
		echo "Error: running short simulation without the long [$BASE_DIR]"; exit 192; 
	fi
	
	DBL_DIRS=(`echo $BASE_DIR/double_*`)

	
	for i in `seq $DBL_STEP`; do
		let "L=L*2"
	done
	
	DIR=$BASE_DIR/double_${DBL_STEP}
	
	if [ $DBL_STEP -eq 1 ]; then
		DBL_STEP=1
		SRC_DIR=$BASE_DIR/long
	else
		let "DBL_PREV=$DBL_STEP-1"
		SRC_DIR=$BASE_DIR/double_$DBL_PREV
	fi
	
	if [ ! -f "$DIR/t\=0_dev\=0.res" ]; then
		SRC_FILE=$SRC_DIR/`get_last_tfile $SRC_DIR` || { echo "File to scale not found: $SRC_DIR/`get_last_tfile $SRC_DIR`"; exit 192; }
		if [ -f ${SRC_FILE} ]; then
			mkdir -p $DIR || { echo "Error creating directory $DIR"; exit 192; } 
			$BIN_DIR/scaleup $SRC_FILE $DIR/t\=0_dev\=0.res || { echo "No file found in ${SRC_DIR}"; exit 192; }
		fi
	fi
else
	DIR=$BASE_DIR/long
fi


if [ $B_EXEC == "efpol" ]; then
	EXEC_LINE="$EXEC $SEED $DIR $DENSITY $TIME $INTERVAL $NMONO $DBL_STEP $L $MODEL"
elif [ $B_EXEC == "denspol" ]; then
	EXEC_LINE="$EXEC $SEED $DIR $DENSITY $TIME $INTERVAL $NMONO $L $CUR_DIR/denspol/ee_topo.dat $BEND_ENERGY"
elif [ $B_EXEC == "gpupol3" ]; then
	EXEC_LINE="$EXEC $NMONO $TIME $SEED $DIR $DENSITY $FAST_EQ $INTERVAL $L $L $L $SHORT $DBL_STEP"
fi

# DESTDIR=$DIR


echo "Seed           : $SEED"
echo "Dest dir       : $DIR"
echo "Max time       : $TIME"
echo "Polymer type   : $SIM_TYPE"
echo "Density        : $DENSITY"
echo "Polymer size   : $NMONO"
echo "Lattice size   : $L"
echo "Double count   : $DBL_STEP"
echo "Interval       : $INTERVAL"
echo "Command        : $EXEC_LINE"

# exit 0
$EXEC_LINE || { echo "Error during polymer simulation execution"; exit 192; }

cat > $BASE_DIR/simulation_settings.txt << EOFCAT
`grep 'RELEASE=' $CUR_DIR/Makefile | sed 's/=/ = /'`
Start_seed = $SEED
Length = $NMONO
Time = $TIME
Polytype = $SIM_TYPE
Start_time = $START_TIME
End_time = `date`
Density = $DENSITY
Latsize = $L
Start_polysize = $NMONO
Double_step = $DBL_STEP
Interval = $INTERVAL
Equilibrated = 0
EOFCAT

# DIRS=$DESTDIR/*/

# for DIR in ${DIRS[*]}; do
# 	cp $DESTDIR/simulation_settings.txt $DIR
# done
	