#!/bin/sh

#  get_score_stat.sh
#
#
#  Created by Justin on 5/10/19.
#
cis_prop=$1
he2=$2
n_causal=$3

task=$SGE_TASK_ID

module load tabix

#gene=$(head -n $SGE_TASK_ID $geneList | tail -n1)

#echo Gene: ${gene}
echo task ${cis} ${he2} ${task}

bfGWAS_SS_dir="/home/jluningham/Projects/bfGWAS_SS"
# first test on a small number of blocks, use alternative filehead file

# mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_LDs
LDdir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments/Simulation/LDs
# Score statistics directory
#mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/Simulation/Score_Stat2
Score_dir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments/Simulation/Score_Stat2
LDwindow=1

cd ${Score_dir}
rm *.core
echo Output/nSNPs_{n_causal}/cis_${cis_prop}/He2_${he2}_trainExpr_${task}.txt

#echo LDwindow is $LDwindow

### With input genotype file in dosage format
${bfGWAS_SS_dir}/bin/Estep_mcmc -g /mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/data/simulation_dosages_all.geno -p /mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/Simulation/Output/nSNPs_${n_causal}/cis_${cis_prop}/He2_${he2}_trainExpr_${task}.txt -o SNPS_${n_causal}_cis_${cis_prop}_He2_${he2}_trainScore_${task} -LDwindow ${LDwindow} -GTfield DS -saveSS -zipSS

mv ${Score_dir}/output/SNPS_${n_causal}_cis_${cis_prop}_He2_${he2}_trainScore_${task}.score.txt.gz ${Score_dir}


#echo task ${task} finished!



