#!/bin/bash

FILE=$1

awk '

function abs (x){
 if(x>=0) return x;
 else return -x;
}
{
	if(abs($(2)-$(1))>max)
		max = abs($(2)-$(1));
	avg_pc[abs($(2)-$(1))] += $(3);
	n[abs($(2)-$(1))]++;
}END{
	for(i=0; i<=max; i++){
		print i, avg_pc[i]/n[i];
	}
}' $FILE
