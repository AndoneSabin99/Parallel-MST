#!/bin/bash
#PBS -l select=1:ncpus=4:mem=2gb
# set max execution time
#PBS -l walltime=0:03:00
# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load mpich-3.2
mpirun.actual -np 4 /home/sabin.andone/Parallel-MST/mst /home/sabin.andone/Parallel-MST/graph.txt 8