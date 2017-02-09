#!/bin/bash

FILE=$1
FILE_OUT=`echo $FILE | sed 's/pc.dat/pc_point.dat/g'`

awk 'BEGIN{l=0}{
for (i=1; i<=NF; i++){
	print l,i-1,$(i);
}
l++;
}' $FILE > $FILE_OUT

# gnuplot << EOF_GNU
# plot "$FILE_OUT" u 1:2:(log(\$3)) w image
# EOF_GNU