#!/bin/bash

#set -x


# Arguments: start, step, end
START=$2
STEP=$3
END=$4

export OMP_NUM_THREADS=6

# Vary grid size at fixed arg
ARG=$1
#for N in $(seq $START $STEP $END); do
#    N=$N ./main_taylor_green $ARG | awk -v n=$N '/L2-norm/{print n, $NF}'
#done


# Vary arg at fixed grid size
N=$1
for ARG in $(seq $START $STEP $END); do
    N=$N ./main_taylor_green $ARG | awk -v a=$ARG '/L2-norm/{print a, $NF}'
done

