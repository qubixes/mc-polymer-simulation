#!/bin/bash

QUEUE=$1; shift

cat << EOFCAT

#!/bin/bash

CUR_DIR=/home/rschram/conring
REC_DIR=/home/rschram/conring/data/pascal_data/long/

#$ -S /bin/bash
#$ -N denspol_ring
#$ -q $QUEUE
#$ -V
#$ -pe mpi8_debian 8
#$ -cwd


# source /usr/share/modules/init/bash
# module use /applis/PSMN/Modules
# module load Base/psmn\\\\\


#SEEDS=("191823" "1298401" "12083701" "1947613" "61928" "16284712" "32948710" "2794812")
SEEDS="191823 1298401 12083701 1947613 61928 16284712 32948710 2794812"
COMMON="-t 1e5 -i 1e3 --hardsphere"
EXTRA_OPTS="$*"

cd \$CUR_DIR

for SEED in \${SEEDS}; do 
	./reconstruct.sh \$REC_DIR \${COMMON} \${EXTRA_OPTS} -s \$SEED &
done

wait

EOFCAT