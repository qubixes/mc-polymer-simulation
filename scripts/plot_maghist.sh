#!/bin/bash


. ./util.sh
DIRS=(`get_dirs $* -e`) || { echo ${DIRS[*]}; exit $?; }

PLOT="plot [:7][1e-8:] "
PLOT2="plot "

for DIR in ${DIRS[*]}; do
	MAGDIP=(`head -1 $DIR/magdip.dat`)
	MAGDIP=${MAGDIP[1]}
	MAG1=(`head -1 $DIR/magdip_hist.dat`)
	MAG2=(`head -2 $DIR/magdip_hist.dat | tail -1`)
	
	DBIN=`echo "$MAG2-$MAG1" | bc -l`
	N=`get_attr 'Length' $DIR`
	
# 	PLOT="$PLOT '$DIR/magdip_hist.dat' u (\$1/$MAGDIP):($MAGDIP*\$3/$DBIN/(\$1/$MAGDIP)**2/(4*pi)):($MAGDIP*\$4/$DBIN/(\$1/$MAGDIP)**2/(4*pi)) w error title 'N=$N', "
	PLOT="$PLOT '$DIR/magdip_hist.dat' u (\$1/$MAGDIP):(func_rem(func($MAGDIP, \$1, \$3, $DBIN), func($MAGDIP, \$1, \$4, $DBIN))):(func($MAGDIP, \$1, \$4, $DBIN)) w error title 'N=$N', "

	PLOT2="$PLOT2 '$DIR/magdip_hist.dat' u ((\$1/$MAGDIP)**0.5):(\$3/$DBIN*2*sqrt(\$1*$MAGDIP)) title 'N=$N', "
done


PLOT=${PLOT:0:${#PLOT}-2}
PLOT2=${PLOT2:0:${#PLOT2}-2}


gnuplot << EOF
set xlabel 'r_m^2/<r_m^2>'
set ylabel 'P_{bin} <r_m^2>^3 / (r_m^4*4*{/Symbol p})'
set log y

func(rm_avg, rm, p, dbin) = (rm_avg*p/dbin/(rm/rm_avg)**2/(4*pi));
func_rem(x,sigma) = (sigma>x*0.3)?(1/0):x;

$PLOT, 0.45*exp(-3*x)
EOF

gnuplot << EOF
set log y
sigma=0.25
mu=1
set key bottom left
$PLOT2, (2*pi*sigma**2)**-0.5*exp(-(x-mu)**2/(2*sigma**2)) title 'gaussian(mu,sigma)'
EOF
# 
# gnuplot <<EOF
# set term 