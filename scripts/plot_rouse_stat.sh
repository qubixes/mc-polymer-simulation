#!/bin/bash
OFILE="rouse_stat_trans"
# ODIR=`./paper_dir.sh`
./self_modes.sh $*

TYPE="ring"
SELF_FILE="self_modes.dat"
TRANS_FILE="trans_modes.dat"
DIV=(5./3.)

PVAL=(`head -1 $SELF_FILE`)

let "NP=${#PVAL[*]}-1"

if [ $TYPE == "lin" ]; then
# 	PLOT="plot [2:5][0.3:0.9]"
	PLOT="plot [2:1e4][0.3:2]"
	B=0.7095350
else
	PLOT="plot [0.1:50][0.02:0.2]"
	B=0.7120154
# 	PLOT="plot [2:5][0.1:0.5]"
# 	PLOT="plot [][] "
fi

# NP=2
IP=2
while [ $IP -lt $NP ]; do
	let "RCOL=IP+3"
#	let "NCOL=1"
#	let "NE_COL=2"
#	let "RE_COL=3"
	P=${PVAL[IP]}
	
	if [ $P -lt 3 ]; then
		TITLE="title 'p=$P'"
# 		LINE="w l"
		LINE=""
	else
		TITLE="notitle"
		LINE=""
	fi
	PLOT="${PLOT} \"$SELF_FILE\" u (\$1/(\$2*$P)):((\$3*\$$RCOL*($P/\$1)**($DIV))) w l $TITLE $LINE,"
	let "IP=IP+1"
done
if [ $TYPE == "lin" ]; then
	PLOT="$PLOT 0.08*x**0.067, 0.097*x**0.022, gauss(x)"
else
	PLOT="$PLOT  f(x), g(x), h(x), fg(x), fgh(x)"
fi

# PLOT="$PLOT (f(x)**(-a)+g(x)**(-a))**(-1.0/a) title 'fit', f(x), g(x)"
gnuplot -persist << EOF
set term aqua enhanced font 'Helvetica, 20'
set key top left
set log x
set log y
a1=-10.0
a2=7.0
nu=0.588


b=$B
gauss(x) = 2**0.5*0.25*b**2/sin(pi/x)**2/x**2;

# f(x) = 0.175*x**(-1.8)
# g(x) = 0.0126*x**(0.26)
# h(x) = 0.19*x**(-1./3.)
# f(x) = 0.015*x**0.45
# g(x) = 0.037*x**0.23333
f(x) = 0.145*x**0.45
g(x) = 0.12*x**0.23333
h(x) = 0.17
fg(x) = (f(x)**a1+g(x)**a1)**(1/a1)
fgh(x) = (fg(x)**-a2+h(x)**-a2)**(-1/a2)

k(x) = 0.0335
m(x) = 0.18*x**(-1./3.)
km(x) = (k(x)**-a1+m(x)**-a1)**(-1/a1)

# set ytics 0.01,1.1*g(0)/0.01, 1.1*g(0)
set ylabel '|X_l|^2 / l^2'
set xlabel 'l'

$PLOT
EOF
