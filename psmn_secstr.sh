#!/bin/bash
. scripts/util.sh

START_N=$1
SUB_DIR=$2
NPROC=$3
QUEUE=$4
NRUN=$5
BASE_DIR=$6

# BASE_DIR=./data/
DIR=$BASE_DIR/ring_gpupol3_l${START_N}_g128_s12396143_d1.2/$SUB_DIR/
CUR_DIR=`pwd`

if [ ! -d $DIR ]; then
	echo "Directory $DIR does not exist"
	exit 102
fi

OUT_DIR=$DIR/sec_struct
mkdir -p $OUT_DIR

NPOL=`get_attr 'Npol' $DIR`
LAST_T=`get_last_t $DIR`

IPOL=0

BATCH_DIR="./batch/N${START_N}_${SUB_DIR}"
mkdir -p $BATCH_DIR
while [ $IPOL -lt $NPOL ]; do
	BATCH_FILE="$BATCH_DIR/batch_${START_N}_${SUB_DIR}_${IPOL}.sh"
	
	cat << EOFC > $BATCH_FILE
#!/bin/bash
#$ -S /bin/bash
#$ -N secstr_${START_N}_${SUB_DIR}_${IPOL}
#$ -q $QUEUE
#$ -V
#$ -pe openmp2 $NPROC
#$ -cwd

cd $CUR_DIR

IPROC=0
CURPOL=$IPOL
IRUN=0
while [ \$IRUN -lt $NRUN -a \$CURPOL -lt $NPOL ]; do
	while [ \$IPROC -lt $NPROC -a \$CURPOL -lt $NPOL ]; do
		./bin/secstr_wrap $DIR ./bin/ $LAST_T \$CURPOL &
		let "IPROC++"
		let "CURPOL++"
	done
	let "IRUN++"
	wait
done

EOFC
	chmod 755 $BATCH_FILE
	let "IPOL+= $NPROC*$NRUN"
done

# 
# source /usr/share/modules/init/bash
# module use /applis/PSMN/Modules
# module load Base/psmn
# 
# 
# echo "qsub $BATCH_FILE"
# chmod 755 $BATCH_FILE
# cat $BATCH_FILE

# #$ -pe openmp2 2
