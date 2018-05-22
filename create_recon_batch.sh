#!/bin/bash

QUEUE=$1

if [ "$QUEUE" == "" ]; then
	echo "Need queue name to create batch"
	exit 192
fi

NSAMPLES=("2e3" "5e3" "1e4" "2e4" "5e4")
# POL_TYPES=("--linear" "--ring")
TOPO_TYPES=(" " "--topo")
HP_STRENGTHS=("0.35" "0.45" "0.55" "0.65")

for HP_STRENGTH in ${HP_STRENGTHS[*]}; do
	for NS in ${NSAMPLES[*]}; do
			let "NTOP_M_ONE=${#TOPO_TYPES[*]}-1"
			for ITOP in `seq 0 $NTOP_M_ONE`; do
				BATCH_FILE=`echo "./batch/batch_ns${NS}_${HP_STRENGTH}$(echo ${TOPO_TYPES[ITOP]}| sed 's/ //g')" | sed 's/--/_/g'`.sh
				./psmn_reconstr.sh $QUEUE --nsamples $NS ${TOPO_TYPES[ITOP]} --hpstrength ${HP_STRENGTH} > $BATCH_FILE
				echo $BATCH_FILE
			done
	done
done