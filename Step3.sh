#!/bin/sh

################################################################
################################################################
# Step 3: EM-MCMC training
# use Make with Perl script to create file
################################################################
################################################################

gene=$1
geneFile=$2
geno_dir=$3
Scripts_dir=$4
Res_dir=$5
LDdir=$6
N=$7
num_cores=$8

# Set up Working directory
mkdir ${Res_dir}/${gene}_GREX
wkdir=${Res_dir}/${gene}_GREX
cd ${wkdir}
#rm core.*

Score_dir=${Res_dir}/${gene}_scores

# get gene info to provide to the makefile arguments

target_chr=$(grep ${gene} ${geneFile} | awk 'FS {print $1}'); echo $target_chr
start_pos=$(grep ${gene} ${geneFile} | awk 'FS {print $2}'); echo $start_pos
end_pos=$(grep ${gene} ${geneFile} | awk 'FS {print $3}'); echo $end_pos

pheno=${Score_dir}/${gene}_exp_dat.txt
pv=$(awk '{delta=$2; sum+=$2; ssq+=(delta - avg/NR)^2} END {print ssq/(NR-1)}' ${pheno})


filehead=${Score_dir}/${gene}_signif_segments.txt
#filehead=${geno_dir}/ROSMAP_GWAS_Seg_filehead2.txt
# Hyper parameter file (a seperate file from hypval.current under the working directory)
hfile=${Scripts_dir}/hypval.txt

# Specify computation specs
em=5; # EM steps
burnin=10000; # Burn-in iterations in MCMC
Nmcmc=10000; # MCMC iteration number
mkfile=${wkdir}/${gene}_BFGWAS.mk

${Scripts_dir}/bin/gen_mkf.pl \
--EMdir ${Scripts_dir}/bin --hyp ${hfile} \
-n ${N} --pv ${pv} \
-w ${wkdir} --geno sumstat \
-f ${filehead} -l local --pheno ${pheno} \
--targ ${target_chr} --start ${start_pos} --end ${end_pos} \
--LDdir ${LDdir} --Scoredir ${Score_dir} \
-j BFGWAS_${gene} --em ${em} -b ${burnin} -N ${Nmcmc} \
--mf ${mkfile}


######### Submit the job for running the makefile (run_Estep.sh and run_Mstep.sh are called by make)
#j=14 # Number of cores to run the job in parallel


make -f ${mkfile} clean


echo Run make with $mkfile and j=$j parallel jobs

make -k -C ${wkdir} -f ${mkfile} -j ${num_cores} > ${wkdir}/make.output 2> ${wkdir}/make.err

exit
