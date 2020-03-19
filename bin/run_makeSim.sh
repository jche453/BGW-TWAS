#!/bin/sh

#  run_makeGW.sh
#
# create makefiles and run bfgwas for a list of genes in parallel
# provide the submission script the list of genes in a text file


cis_prop=$1
he2=$2
n_causal=$3
task=$SGE_TASK_ID

pv=0.99


bfGWAS_SS_dir="/mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS"
# first test on a small number of blocks, use alternative filehead file

# mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_LDs
LDdir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments/Simulation/LDs
# Score statistics directory
#mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_scores
Score_dir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments/Simulation/Score_Stat2
LDwindow=1

cd ${Score_dir}

echo SNPS_${n_causal}_cis_${cis_prop}_He2_${he2}_trainScore_${task} > SNPS_${n_causal}_cis_${cis_prop}_He2_${he2}_trainScore_${task}_filehead.txt

filehead=SNPS_${n_causal}_cis_${cis_prop}_He2_${he2}_trainScore_${task}_filehead.txt


# get train, get test sets
# match IDs to gene_exp_dat
N=499


EMdir=/mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation

# Hyper parameter file (a seperate file from hypval.current under the working directory)
hfile=${EMdir}/hypval.txt

# Specify computation specs
em=5; # EM steps
burnin=10000; # Burn-in iterations in MCMC
Nmcmc=10000; # MCMC iteration number
study=SNPS_${n_causal}_cis_${cis_prop}_He2_${he2}_trainScore_${task}
mkdir /mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/Results/${study}

wkdir=/mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/Results/${study}

mkfile=${wkdir}/${study}.mk

${bfGWAS_SS_dir}/bin/gen_mkf.pl \
--EMdir ${EMdir} --hyp ${hfile} \
-n ${N} --pv ${pv} \
-w ${wkdir} --geno sumstat --gd ${geno_dir} \
-f ${filehead} -l local --pheno /mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/Output/nSNPs_${n_causal}/cis_${cis_prop}/He2_${he2}_trainExpr_${task}.txt \
--targ 19 --start 1040101 --end 1065571 \
--LDdir ${LDdir} --Scoredir ${Score_dir} \
-j BFGWAS_${study} --em ${em} -b ${burnin} -N ${Nmcmc} \
--mf ${mkfile}

cd ${wkdir}
echo run simulation study
make -f ${mkfile} clean


echo Run make with $mkfile

make -k -C ${wkdir} -f ${mkfile} -j 1 > ${wkdir}/make.output 2> ${wkdir}/make.err

exit

#!/bin/sh

#  backup.sh
#
#
#  Created by Justin on 5/13/19.
#
