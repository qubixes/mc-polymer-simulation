#!/bin/bash

L="128"
EXEC="gpupol3"
D="1.2"
DIR="/local/schram/gpupol_data"

EX_ARG="-g $L -x $EXEC -d $D -r $DIR"

./do_run.sh --time 1e5 --interval 1e2 --nmono 30   -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e2 --interval 1   --nmono 30   -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 5e5 --interval 1e3 --nmono 60   -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e3 --interval 10  --nmono 60   -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 1e6 --interval 1e3 --nmono 100  -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e3 --interval 10  --nmono 100  -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 3e6 --interval 1e4 --nmono 150  -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e4 --interval 1e2 --nmono 150  -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 5e6 --interval 1e4 --nmono 200  -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e4 --interval 1e2 --nmono 200  -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 3e7 --interval 1e5 --nmono 300  -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e5 --interval 1e3 --nmono 300  -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 7e7 --interval 1e5 --nmono 500  -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e5 --interval 1e3 --nmono 500  -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 2e8 --interval 1e5 --nmono 700  -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e5 --interval 1e3 --nmono 700  -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 3e8 --interval 1e6 --nmono 1000 -g $L -x $EXEC -d $D -r $DIR
./do_run.sh --time 3e6 --interval 1e4 --nmono 1000 -g $L -x $EXEC -d $D -r $DIR --short

./do_run.sh --time 5e8 --interval 1e6 --nmono 1500 $EX_ARG
./do_run.sh --time 3e6 --interval 1e4 --nmono 1500 $EX_ARG --short

./do_run.sh --time 1e9 --interval 1e6 --nmono 2000 $EX_ARG
./do_run.sh --time 3e6 --interval 1e4 --nmono 2000 $EX_ARG --short

./do_run.sh --time 3e9 --interval 1e7 --nmono 3000 $EX_ARG 
./do_run.sh --time 3e7 --interval 1e5 --nmono 3000 $EX_ARG --short

# wait