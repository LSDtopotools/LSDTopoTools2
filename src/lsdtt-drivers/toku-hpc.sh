#!/bin/bash
#$ -cwd
#$ -j y
#$ -N Tokunaga
#$ -o /data/home/faw513/Alex-Toku/LSDTopoTools2/src/lsdtt-drivers
#$ -pe smp 1
#$ -l h_vmem=16G
#$ -l h_rt=0:20:0


cd /data/home/faw513/Alex-Toku/LSDTopoTools2/src/lsdtt-drivers

./strahler-hpc.out /data/home/faw513/Alex-Toku/LSDTopoTools2/src/lsdtt-drivers/ /data/home/faw513/Alex-Toku/LSDTopoTools2/src/lsdtt-drivers/ OCR
