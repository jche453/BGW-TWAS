#!/bin/sh

#geneList=$1

#gene=$(head -n $SGE_TASK_ID $geneList | tail -n1)
gene=${1}
Res_dir=${2}
geno_dir=${3}
pred_geno_filenames=$4
pheno_file=${5}
GT=${6}
tab_mod=${7}

wkdir=${Res_dir}/${gene}_GREX
cd ${wkdir}
#module load tabix/0.2.6
module load ${tab_mod}


## make header row for vcf file from one of the existing genotype files
head -n 1 ${geno_dir}/${pred_geno_filenames} | zcat | head -n 300 | grep CHROM > ${wkdir}/${gene}_grex_signal.vcf

## Loop through SNP ID and use tabix
cat ${wkdir}/Eoutput/grexsnps.txt | while read line ; do
chr=$(echo $line | awk -F[:_/] '{print $1}' )
start=$(echo $line | awk -F[:_/] '{print $2}' )
ref=$(echo $line | awk -F[:_/] '{print $3}' )
alt=$(echo $line | awk -F[:_/] '{print $4}' )
echo ${chr}:${start}-${start}
file_temp=$(grep CHR${chr} $pred_geno_filenames | awk '{print }' )
tabix {geno_dir}/${file_temp}.vcf.gz ${chr}:${start}-${start} | awk -v ref=${ref} '$4==ref{print }' >> ${wkdir}/${gene}_grex_signal.vcf
done

## Generate genotype file (dosages) with cesdsum_meta_signal.vcf
${Scripts_dir}/bin/Estep_mcmc -vcf ${wkdir}/${gene}_grex_signal.vcf -p ${pheno_file} -o ${gene}_grex_genotypes -GTfield ${GT} -saveGeno -maf 0.01


## genotype file will be saved under ./output/${gene}_grex_genotypes.geno

## Remove the "#" from the header row

sed -i 's/#//g' ./output/${gene}_grex_genotypes.geno
mv ./output/${gene}_grex_genotypes.geno ${wkdir}

exit

