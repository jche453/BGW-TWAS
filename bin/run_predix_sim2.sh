#!/bin/sh

he2=$2
causal_prop=$1
n_causal=$3
n_thread=$4

# i=$SGE_TASK_ID

module load R

Rscript /mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/run_predixcan_simuProp.R ${he2} ${causal_prop} ${n_causal} ${n_thread}

exit
