#!/usr/bin/bash

wkdir=$1
line=$2
N=$3
pv=$4
burnin=$5
Nmcmc=$6

# summary stat directories
LDdir=$7
Scoredir=$8

# Input files
hfile=$9
start_pos=${10}
end_pos=${11}
target_chr=${12}
rand_seed=$RANDOM


# echo Curerent working directory
# pwd

#echo Run Estep with $line $N $pv $hfile $start_pos $end_pos $target_chr

/home/jchen/BayesianTWAS/bin/Estep_mcmc -inputSS \
-score ${Scoredir}/${line}.score.txt.gz \
-LDcorr ${LDdir}/${line}.LDcorr.txt.gz \
-target_chr ${target_chr} -start_pos ${start_pos} -end_pos ${end_pos} \
-hfile ${hfile} -max_iter 100 -convergence 1e-20 \
-n ${N} -pv ${pv} -maf 0.01 -r2 0.0001 -smin 0 -smax 100 -win 100 \
-o ${line} -w ${burnin} -s ${Nmcmc} -bvsrm -initype 3 -seed ${rand_seed} \
> ${wkdir}/OUT/${line}.output.txt

exit
