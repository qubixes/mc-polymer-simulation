#!/bin/bash
OFILE="rouse_stat_trans"
# ODIR=`./paper_dir.sh`
./self_modes.sh $*

TYPE="ring"
SELF_FILE="self_modes.tmp"
TRANS_FILE="trans_modes.tmp"
DIVS=(0 5./3.)
DIV_TITLES=("0" "5/3")
BOUNDS=("" "0.04:0.2")

PVAL=(`head -1 $SELF_FILE`)

let "NP=${#PVAL[*]}-1"


IDIV=0
for DIV in ${DIVS[*]}; do

if [ ${DIV_TITLES[IDIV]} == "0" ]; then
	DIV_TITLE_TOT=""
	EX_FIT=", (x>5 && x<30)?(0.27*x**(5./3.)):(1/0) lw 2 dt 2 lc rgb \"black\" notitle"
	EX_COMMAND="set label \"l_z^{5/3}\" at 10,40"
else
	DIV_TITLE_TOT="/ l_z^{${DIV_TITLES[IDIV]}}"
	EX_COMMAND=""
	EX_FIT=""
fi


if [ $TYPE == "lin" ]; then
# 	PLOT="plot [2:5][0.3:0.9]"
	PLOT="plot [2:1e4][0.3:2]"
	B=0.7095350
else
# 	PLOT="plot [0.1:200][0.02:0.2]"
	B=0.7120154
# 	PLOT="plot [2:5][0.1:0.5]"
	PLOT="plot [0.2:200] "
# 	echo "bounds: $IDIV ${BOUNDS[IDIV]}"
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
	PLOT="${PLOT} \"$SELF_FILE\" u (\$1/(\$2*$P)):((\$3*\$2**(-5./3.)*\$$RCOL*(\$2*$P/\$1)**($DIV))) w l $TITLE $LINE,"
	let "IP=IP+1"
done
if [ $TYPE == "lin" ]; then
	PLOT="$PLOT 0.08*x**0.067, 0.097*x**0.022, gauss(x)"
else
	PLOT="$PLOT fghi(x)*(x)**(5./3.-$DIV) title \"fit\""
fi

# echo $PLOT
# PLOT="$PLOT (f(x)**(-a)+g(x)**(-a))**(-1.0/a) title 'fit', f(x), g(x)"
gnuplot -persist << EOF
set term aqua enhanced dashed font 'Times Roman, 24'
set key bottom right
set log x
set log y
set yrange [${BOUNDS[IDIV]}]
set xrange [0.01:10]
a1=-10.0
a2=7.0
a3=2.34
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
i(x) = 0.000055*x**-(5./3.)

fg(x) = (f(x)**a1+g(x)**a1)**(1/a1)
fgh(x) = (fg(x)**-a2+h(x)**-a2)**(-1/a2)
fghi(x) = (fgh(x)**a3+i(x)**a3)**(1/a3)
k(x) = 0.0335
m(x) = 0.18*x**(-1./3.)
km(x) = (k(x)**-a1+m(x)**-a1)**(-1/a1)

# set ytics 0.01,1.1*g(0)/0.01, 1.1*g(0)
set ylabel '|X_{l_z}|^2 $DIV_TITLE_TOT * c'
set xlabel 'l_z'

$EX_COMMAND
$PLOT $EX_FIT
# print fghi(0.03)
# plot fgh(x), i(x), fghi(x)
EOF
let "IDIV++"
done
