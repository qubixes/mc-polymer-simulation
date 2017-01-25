#!/bin/bash

#SBATCH -N 1
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH --constraint=haswell

module load cuda



NT_MAX=48
cd ~/gpupol_cpu/bin

for i in `seq $NT_MAX`; do
	./bench --cpu $i >> ../bench.txt
done

#./bench --opencl 1 >> ../bench.txt
#./bench --cuda 1 >> ../bench.txt
