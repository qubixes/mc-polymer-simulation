#!/bin/bash

# . util.sh

DIRS=../data/ring_gpupol*/short*/

for DIR in $DIRS; do
# 	if [ "`grep 'Equilibrated = 0' "$DIR/simulation_settings.txt"`" != "" ]; then
# 		echo "$DIR"
# 		cat "$DIR/simulation_settings.txt" | sed 's/Equilibrated = 0/Equilibrated = 1/g' > temp
# 		cp temp "$DIR/simulation_settings.txt"
# 	fi
done
rm -f temp