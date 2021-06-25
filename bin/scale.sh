#!/bin/bash

dp=''  # path to data folder
graph=''  # graph file
reads=()  # list of read files
blocks=('136' '68' '34' '16' '8' '4' '2' '1')

echo "" > hga.log

# Varied number of GPUs
for reads in "${reads[@]}"; do
    for ((t = 8; t >= 1; t /= 2));do
        echo "> Run Read ${reads} Devices ${t}" >> hga.log
        /usr/bin/time ./hga -g ${dp}/${graph} -r ${dp}/${reads} -x 1 -b 68 -t 128 -d ${t} >> hga.log
    done
done

# Varied number of threads per block
for ((t = 128; t >= 1; t /= 2));do
    for reads in "${reads[@]}"; do
        echo "> Run Read ${reads} Threads ${t}" >> hga.log
        /usr/bin/time ./hga -g ${dp}/${graph} -r ${dp}/${reads} -x 1 -b 68 -t ${t} -d 1 >> hga.log
    done
done

# Varied number of thread blocks
for t in "${blocks[@]}";do
    for reads in "${reads[@]}"; do
        echo "> Run Read ${reads} Blocks ${t}" >> hga.log
        /usr/bin/time ./hga -g ${dp}/${graph} -r ${dp}/${reads} -x 1 -b ${t} -t 128 -d 1 >> hga.log
    done
done
