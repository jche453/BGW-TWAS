#!/bin/sh

################################################################
################################################################
# Step 1: obtain summary statistics (aka Score Statistics)
# Runs single-variant GWAS on the training sample genotypes and available expression data.
# first extracts gene location information for target gene from geneFile.
################################################################
################################################################

gene=$1
geneFile=$2
geno_dir=$3
Scripts_dir=$4
Res_dir=$5
LDdir=$6
Genome_Seg_File=$7
num_segments=$8
num_cores=$9

echo ${gene}
# Score statistics directory
mkdir ${Res_dir}/${gene}_scores
Score_dir=${Res_dir}/${gene}_scores
cd ${Score_dir}
# geneFile columns are Gene Name, Chr, Pos, start, end, expr_data
# the following creates a phenotype file for target gene that includes subject IDs and expression scores only.
head -1 ${geneFile} | awk '{$1=$2=$3=$4=$5=""; print substr($0,6)}' | tr ' ' '\n' > temp_ID.txt
grep ${gene} ${geneFile} | awk '{$1=$2=$3=$4=$5=""; print substr($0,6)}' | tr ' ' '\n' > exp_temp.txt
paste temp_ID.txt exp_temp.txt > ${gene}_exp_dat.txt
pv=$(awk '{delta=$2; sum+=$2; ssq+=(delta - avg/NR)^2} END {print ssq/(NR-1)}' ${gene}_exp_dat.txt)
echo expression pheno variance = $pv
paste ${gene} ${pv} >> ${Res_dir}/geneExp_var.txt
rm temp_ID.txt exp_temp.txt

pheno=${gene}_exp_dat.txt

seq 1 ${num_segments} | xargs -I{block} -n 1 -P ${num_cores} sh ${Scripts_dir}/get_score_stat.sh ${pheno} ${geno_dir} ${Score_dir} ${Scripts_dir} ${LDdir} ${Genome_Seg_File} ${block}

echo Step 1 complete!

exit
