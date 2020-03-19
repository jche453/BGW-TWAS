#!/bin/sh

he2=$2
causal_prop=$1
n_thread=$3

# i=$SGE_TASK_ID

module load R

Rscript /mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/run_predixcan_simu.R ${he2} ${causal_prop} ${n_thread}

exit
