#!/bin/bash

DIR="../data"

gnuplot -persist << EOFGNU
set log x
set log y

plot [][0.04:0.1] "$DIR/rgyr_d1.0.dat" u (0.9*\$1):(1.011*\$2/\$1) , "$DIR/rgyr_d1.5.dat" u (\$1*1.27):(1.26*\$2/\$1) , "$DIR/../gpupol_data/rgyr_ring.dat" u (\$1*1.9):(0.827*\$2/\$1), 0.666*x**(-1./3.), "$DIR/rgyr_d2.0.dat" u (\$1*1.0):(1.56*\$2/\$1), "$DIR/rgyr_d1.1.dat" u (1.75*\$1):(0.92*\$2/\$1), "$DIR/rgyr_d1.6.dat" u (1.75*\$1):(1.18*\$2/\$1), "$DIR/rgyr_d0.9.dat" u (2.5*\$1):(0.74*\$2/\$1), "$DIR/rgyr_d1.2.dat" u (2.3*\$1):(0.89*\$2/\$1)

EOFGNU