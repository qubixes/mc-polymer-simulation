#!/bin/bash

. ./util.sh
DIRS=(`get_dirs $* -e`) || { echo ${DIRS[*]}; exit $?; }
TYPE="ring"

OFILE1="rouse_dyn_$TYPE"
OFILE2="rouse_dyn_lin_$TYPE"

PLOT="plot [][]"
PLOT_LIN="plot [:1e1][]"

J=0
SET_LABEL=0
USE_PCORR=1;
PFAC=2.15;

if [ $USE_PCORR == "1" ]; then
	PLOT="plot [12:2e8][]"
else
	PLOT="plot [][]"
fi
for DIR in ${DIRS[*]}; do
	FILE="$DIR/rouse_dyn.dat"
	LENGTH=`get_attr 'Length' "$DIR"`
	TFILE="./rtemp_${LENGTH}.txt"
	
	VAR=(`awk '{
		if( $1 !~ /#/ ){
			for(i=2; i<=NF; i++){
				if(zero[i]>0 && zero[i]<12){
					var[i] += $(i)*$(i);
					zero[i]++;
				}
				else if($(i) <= 0 && !zero[i]){
					zero[i]=1;
					var[i] += $(i)*$(i);
				}
			}
		}
		}
		END{
		for(i=0; i<=NF; i++){
			if(!zero[i])
				printf "0 ";
			else
				print sqrt(var[i]/zero[i]);
		}
		}
		' $FILE`)
# 	echo $DIR
# 	echo ${VAR[*]}
# 	exit 0
	awk -v str="${VAR[*]}" 'BEGIN{
	split(str, var, " ");
	}
	{
		if( $1 !~ /#/ ){
			printf "%i ", $1;
			for(i=2; i<=NF; i++){
				if(zero[i])
					printf "0 ";
				else if($(i) <= var[i] ){
					zero[i]=1;
					printf "0 ";
				}
				else
					printf "%le ", $(i);
			}
			printf "\n";
		}
		else
			print $0;
		}
# 	END{ for(i=0; i<=NF; i++) printf "%s ", var[i]
# 		printf "\n"
# 		print str
# 		}
		' $FILE > $TFILE
# 	cat $TFILE
# 	exit 0
	
	#tau ~ p^beta
	PNORM=(`head -2 $TFILE | tail -1`)
	PVAL=(`head -1 $TFILE`)
	let "PMAX=${#PVAL[*]}-2"
	for IP in `seq 1 $PMAX`; do
		P=${PVAL[IP+1]}
		let "COL=IP+2"
		if [ $SET_LABEL -eq 1 ]; then
			TITLE="title '\$p=${P}\$'"
		else
			TITLE="notitle"
		fi
		PLOT="$PLOT '$TFILE'  u (\$1*($USE_PCORR?$P**$PFAC:1)):(log(-log(\$$COL/${PNORM[COL-1]})) + beta*log($P/(1.0*${LENGTH})) + delta*log(\$1)/log(1+$P) + $USE_PCORR*0.85*$PFAC*log($P)) w l $TITLE lt $P, "
		PLOT_LIN="$PLOT_LIN '$TFILE' u (\$1*10**-6):( log(\$$COL/${PNORM[COL-1]})*($P/(1.0*${LENGTH}))**beta) w l $TITLE lc $IP, "
	done
	SET_LABEL=0
	let "J=J+1"
done
# exit 0
# PLOT="$PLOT x**0.7*25, x*25"
PLOT=${PLOT:0:${#PLOT}-2}
PLOT_LIN=${PLOT_LIN:0:${#PLOT_LIN}-2}
# echo $PLOT
gnuplot -persist << EOF
set term aqua enhanced font 'Helvetica, 22'
set log x
# set log y
set key top left
k=4.0
Ne=100
f(t,l) = (t/l**2)
g(t,l) = (t/l**2.5+Ne**2*(l**-2-l**-2.5))
time_tild(t,l) = ((t/l**2)**-k+(t/l**2.5+Ne**2*(l**-2-l**-2.5))**-k)**(-1/k)
set xlabel 'log(t)'
set ylabel 'log(-log(<X_l(t)*X_l(t+dt)>/<X_l(t)*X_l(t)>))+{/Symbol b} log (l)'
set grid
beta=-1.88
delta=-0.0
# plot [Ne**2/100:Ne**2*100] time_tild(x, 500), f(x, 500), g(x, 500)
h(x) = exp(5.02)*x**0.61
i(x) = exp(1.32)*x**0.85
hi(x) = (h(x)**-4+i(x)**-4)**(-1./4.)
$PLOT, log(h(x)) lw 1.5 notitle, log(i(x)) lw 1.5 notitle
EOF

# gnuplot << EOFGNU
# set term epslatex color standalone colortext 14
# beta=-1.88
# delta=-0.0
# set key at 8.5,19
# set ylabel '\$\\log\\left[-\\log\\left(\\frac{\\left<\\mathbf{X}_\\ell(t)\\cdot \\mathbf{X}_\\ell(t+dt)\\right>}{\\left<\\mathbf{X}_\\ell(t)\\cdot \\mathbf{X}_\\ell(t)\\right>}\\right)\\right]+ \\beta \\log \\ell\$'
# set xlabel '\$\\log t\$'
# # set key spacing 1.5
# set out "$OFILE1.tex"
# f(x) = x**0
# g(x) = x**-1
# h(x) = x**1.0
# $PLOT, 0.65*x+4.3 lw 2 notitle, 0.85*x+1.2 notitle
# EOFGNU
exit 0





# 
# gnuplot -persist << EOF
# set term aqua enhanced font "Helvetica, 22"
# beta=-1.88
# delta=-0.0
# set xlabel 't (*10^6)'
# set ylabel 'log(<X_l(t)*X_l(t+dt)>/<X_l(t)*X_l(t)>)*l^{/Symbol b}'
# set grid
# $PLOT_LIN, -exp(4.3)*(x*10**6)**0.65 lw 2 lc rgb 'black' title "-74 t^{0.65}"
# EOF
# 
# gnuplot << EOFGNU
# set term epslatex color standalone colortext 14
# beta=-1.88
# delta=-0.0
# # set key at 7.5,19
# set key off
# set ylabel '\$\\log\\left(\\frac{\\left<\\mathbf{X}_\\ell(t)\\cdot \\mathbf{X}_\\ell(t+dt)\\right>}{\\left<\\mathbf{X}_\\ell(t)\\cdot \\mathbf{X}_\\ell(t)\\right>}\\right)\\ell^\\beta\$'
# set xlabel '\$t (\\times 10^6)\$'
# # set key spacing 1.5
# set out "$OFILE2.tex"
# $PLOT_LIN, -exp(4.3)*(x*10**6)**0.65 lw 2 lc rgb 'black' title "\$-74 t^{0.65}\$"
# EOFGNU

# latex $OFILE1.tex
# dvips -E -o full_$OFILE1.eps $OFILE1.dvi
# mv full_$OFILE1.eps $OFILE1.eps
# 
# latex $OFILE2.tex
# dvips -E -o full_$OFILE2.eps $OFILE2.dvi
# mv full_$OFILE2.eps $OFILE2.eps
# 
exit 0
make {$OFILE1,$OFILE2}.eps
cp {$OFILE1,$OFILE2}.eps $ODIR

rm ./rtemp_*