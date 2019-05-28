#!/bin/bash
#$ -cwd
#$ -j y
#$ -N Tokunaga
#$ -o /data/Geog-c2s2/toku
#$ -pe smp 1
#$ -l node_type=sm
#$ -l h_vmem=256G
#$ -l h_rt=12:0:0
#$ -l highmem

./data/home/faw513/LSDTopoTools2/src/lsdtt-drivers/strahler-hpc.out
