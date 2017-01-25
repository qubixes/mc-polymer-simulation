#!/bin/bash

NT_MAX=16
cd ./bin

for i in `seq $NT_MAX`; do
	./bench --cpu $i >> ../bench.txt
done

./bench --opencl 1 >> ../bench.txt
./bench --cuda 1 >> ../bench.txt
