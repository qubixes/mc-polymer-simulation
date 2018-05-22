#!/bin/bash

. util.sh

FIN_POV=final.pov
TEMP_POV=temp.pov
HEAD_POV=head.pov
PNG_FILE=pol3d.png

if [ $# -gt 1 ]; then
	PNG_FILE=$2
fi

if [ $# -gt 2 ]; then
	POLID=$3
fi


function last_modification {
	if [ `uname` == "Darwin" ]; then
		stat -f "%m" $1
	else
		date -r $1 +%s
	fi
}

if [ $# -lt 1 ]; then
	echo "Need one argument: either a directory or file to plot"
	exit 193;
fi

if [ -d $1 ]; then
	FILE_TUV=$1/`get_last_tfile $1`
else
	FILE_TUV=$1
fi

../bin/tuv2pov $FILE_TUV $TEMP_POV $POLID || { echo "../bin/tuv2pov $FILE_TUV $TEMP_POV"; exit $?; }


#Build title


cat $HEAD_POV $TEMP_POV > $FIN_POV || { echo "Error opening files $HEAD_POV, $TEMP_POV"; exit 192; }
povray $FIN_POV -V -P -D +W600 +H600 +A0.3 +AM1 +Q2 +O$PNG_FILE &> /dev/null || exit 192

gnuplot << EOFGNU
# set title "$TITLE" offset 0,-2
# set label "$TITLE" at 0,0 front
unset border
unset xtics
unset ytics
set size ratio 1
set term aqua size 600 600 
plot "$PNG_FILE" binary filetype=png w rgbimage notitle
EOFGNU