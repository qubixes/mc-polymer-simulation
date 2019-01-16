#!/bin/bash

FILE=$1
POLID=$2

TFILE="xyz.tmp"

../bin/tuv2xyz $FILE $TFILE $POLID

gnuplot <<EOFGNU

set term x11
splot "$TFILE" u 1:2:3 w l

pause mouse keypress
EOFGNU