#!/bin/bash

HOSTS=(`cat resources.dat`)

THRES="25"
MAIN=`hostname`


declare -i USE_THRES=0
if [ $# -gt 0 ]; then
	USE_THRES=$1
fi

for HOST in ${HOSTS[*]}; do
	ssh -t -o ConnectTimeout=3 -o ServerAliveInterval=3 -o ServerAliveCountmax=1 $HOST << EOFSSH 2> /dev/null  | grep "$HOST"

CPU_UTIL=\`top -b -n2 -p 1 | fgrep "Cpu(s)" | tail -1 | awk -F'id,' -v prefix="\$prefix" '{ split(\$1, vs, ","); v=vs[length(vs)]; sub("%", "", v); printf "%s%.1f%%\n", prefix, 100 - v }'\`

echo "$HOST": "\$CPU_UTIL", \`grep  processor /proc/cpuinfo | wc -l \` cores

# if [ \${EMPTY} -eq 1  -a ${USE_THRES} -ne 1 ]; then
# 	echo \`hostname\`
# fi
# 
# if [  ${USE_THRES} -eq 1 ]; then
# 	echo "$HOST" \$CPU_UTIL%, , user=\${USER[1]}
# fi
EOFSSH
done
