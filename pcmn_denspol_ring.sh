#!/bin/bash

CUR_DIR=/home/rschram/conring
DEST_DIR=/scratch/e5-267-gfs-scratch/schram/denspol_data/

#$ -S /bin/bash
#$ -N denspol_ring
#$ -q E5-2667v4deb128nl
#$ -V
#$ -pe openmp2 1
#$ -cwd


source /usr/share/modules/init/bash
module use /applis/PSMN/Modules
module load Base/psmn


DENSITY=7.2
BEND=0.3

COMMON="-x denspol -d $DENSITY -b $BEND --ring --dir $DEST_DIR"

cd $CUR_DIR
./do_run.sh --nmono 30  -g 15 --time 1e5 --interval 1e2 $COMMON &
./do_run.sh --nmono 50  -g 15 --time 2e5 --interval 1e2 $COMMON &
./do_run.sh --nmono 70  -g 15 --time 4e5 --interval 1e3 $COMMON & 
./do_run.sh --nmono 100 -g 15 --time 2e6 --interval 1e3 $COMMON &
./do_run.sh --nmono 150 -g 15 --time 5e6 --interval 1e4 $COMMON &
./do_run.sh --nmono 200 -g 15 --time 1e7 --interval 1e4 $COMMON &
./do_run.sh --nmono 300 -g 20 --time 3e7 --interval 1e5 $COMMON &
./do_run.sh --nmono 400 -g 20 --time 6e7 --interval 1e5 $COMMON &

wait

