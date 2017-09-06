#!/bin/bash

L=$1

DIR=../data/ring_gpupol3_l${L}_g128_s12396143_d1.2/long/
LENGTH=(`grep 'len' $DIR/t=0_dev=0.res | head -1`)
LENGTH=${LENGTH[1]}
OUT_DIR=$DIR/sec_struct
mkdir -p $OUT_DIR
sed 's/RING_LENGTH_SUB/'$LENGTH'/' ./input_template.dat > $OUT_DIR/input.dat

./secstr_wrap $DIR ./ 10000000 2