#!/bin/sh

#  get_score_stat.sh
#  
#
#  Created by Justin on 5/10/19.
#

pheno=$1
geno_dir=$2
Score_dir=$3
Scripts_dir=$4
LDdir=$5
filehead=$6
block=$7

# Score statistics directory
LDwindow=1000000
cd ${Score_dir}

line=$(head -n $block $filehead | tail -n1)

if [ -f ${LDdir}/${line}.LDcorr.txt.gz ] ; then
#echo ${LDdir}/${line}.LDcorr.txt.gz exists!
LDwindow=1
fi

#echo LDwindow is $LDwindow

### With input genotype file in dosage format
${Scripts_dir}/bin/Estep_mcmc -vcf ${geno_dir}/${line}.vcf.gz -p ${pheno} -maf 0 -o ${line} -LDwindow ${LDwindow} -GTfield DS -saveSS -zipSS

mv ${Score_dir}/output/${line}.score.txt.gz ${Score_dir}

echo task ${task} finished!



