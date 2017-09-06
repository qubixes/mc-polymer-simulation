#!/bin/bash
. ./util.sh

DIRS=(`get_dirs $* -e`) || { echo "${DIRS[*]}"; exit $?; }

NDIR=${#DIRS[*]}

MODE_FILE="self_modes.tmp"
TRANS_FILE="self_modes_trans.tmp"

rm -f $MODE_FILE $TRANS_FILE

PLINE=(`head -1 "${DIRS[${#DIRS[*]}-1]}/rouse_stat.dat"`)
# echo $PLINE "${DIRS[${#DIRS[*]}-1]}/rouse_stat.dat"
let "NP_MAX=${#PLINE[*]}-1"
echo ${PLINE[*]} >> $MODE_FILE
J=0
for DIR in ${DIRS[*]}; do
	FILE="$DIR/rouse_stat.dat"
	LENGTH=`get_attr 'Length' "$DIR"`
	PVAL=(`head -1 $FILE`)
	let "NP=${#PVAL[*]}-1"
	
	I=0
	DENSITY=`get_attr 'Density' $DIR`
	EXEC=`get_attr 'Executable' $DIR`
	LINE=(`grep '^'"$EXEC ${DENSITY:0:3}" "ne_list.dat"`)
	NE_FAC=${LINE[2]}
	R_FAC=${LINE[3]}
# 	echo ${LINE[1]} ${LINE[*]} $DENSITY
	printf "%i %f %f " $LENGTH $NE_FAC $R_FAC
# 	echo $FILE
	while read -r LINE ; do
		if [ $I -eq "0" ]; then
			let "I=I+1"
			continue;
		fi
		
		X=($LINE)
		
		printf "%lf " ${X[I-1]}
		
		let "I=I+1"
	done < "$FILE"
	let "NNAN=NP_MAX-NP"
	if [ $NNAN -gt 0 ]; then
		for K in `seq $NNAN`; do 
			printf "%s " "nan"
		done
	fi
	echo ""
	let "J=J+1"
done >> $MODE_FILE


awk ' {
if(iLine==0){
	for(i=2; i<=NF; i++){
		p[NP++]= $(i);
	}
}
else{
	len[iLine-1] = $(1);
	for(i=2; i<=NF; i++){
		mat[iLine-1, i-2] = $(i);
	}
}
iLine++
}
END{
	NL=iLine-1;
	printf "# "
	for(i=0; i<NL; i++) printf "%i ", len[i]
	printf "\n"
	for(iP=1; iP<NP; iP++){
		printf "%lf ", 1/p[iP]
		for(iL=0; iL<NL; iL++){
			printf "%f ", mat[iL, iP];
		}
		printf "\n"
	}
	}
' $MODE_FILE > $TRANS_FILE

