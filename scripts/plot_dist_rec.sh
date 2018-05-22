#!/bin/bash

. util.sh

REC_DIR=$1
PAS_DIR=$REC_DIR/../../long/

REC_FILE=$REC_DIR/`get_last_tfile $REC_DIR`

BIN_FILE=$PAS_DIR/bins*
CONFIG_FILE=$PAS_DIR/config*

OUT_FILE="$REC_DIR/rec_comparison.dat"
PAS_HIST_FILE="$REC_DIR/exp_hist.dat"
SIM_HIST_FILE="$REC_DIR/sim_hist.dat"

if [ ! -f $OUT_FILE ]; then
	../bin/detailed_comp_pascal $BIN_FILE $CONFIG_FILE $REC_FILE $OUT_FILE $PAS_HIST_FILE $SIM_HIST_FILE || { echo "../bin/detailed_comp_pascal $BIN_FILE $CONFIG_FILE $REC_FILE $OUT_FILE $PAS_HIST_FILE $SIM_HIST_FILE"; exit 192; }
fi

gnuplot --persist << EOF_GNU

plot "$OUT_FILE" u 1:2:3 w errorbars
EOF_GNU

gnuplot --persist << EOF_GNU
plot "$OUT_FILE" u 4:5:6 w errorbars, x, 0.95*x
EOF_GNU


PRE_FAC=0.07
gnuplot --persist << EOF_GNU
plot "$PAS_HIST_FILE" u 1:(\$2/\$1**0.5) w l, "$SIM_HIST_FILE" u (\$1/$PRE_FAC):(\$2/(\$1/$PRE_FAC)**0.5) w l
EOF_GNU

# gnuplot --persist << EOF_GNU
# plot "$SIM_HIST_FILE" u 1:2 w l
# EOF_GNU

