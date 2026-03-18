#!/bin/bash
# Minimalist MLUPS benchmark

echo "N | 1 Thread | 6 Threads"
echo "--------------------------"

for N in $(seq 10 5 150); do
    # Run Single-Threaded
    MLUPS_1=$(OMP_NUM_THREADS=1 N=$N ./main_taylor_green 0.2 | awk '/MLUPS/{print $NF}')
    
    # Run Multi-Threaded (using 4 threads as an example)
    MLUPS_2=$(OMP_NUM_THREADS=6 N=$N ./main_taylor_green 0.2 | awk '/MLUPS/{print $NF}')
    
    echo "$N | $MLUPS_1 | $MLUPS_2"
done

