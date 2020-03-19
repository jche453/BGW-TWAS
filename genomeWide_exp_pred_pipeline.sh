## Workflow for genome-wide bfGWAS for each of 14772 genes.

## first step: set up score statistics for each gene
##geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/random_selected_genes.txt" ##used a subset of genes to first apply tool, Oct/Nov 2018
geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/pilot_gene_list.txt" ##used a subset of genes to first apply tool, Oct 2018

wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"

cd ${wkdir}

gene=${LINE}
echo ${gene}
bfGWAS_SS_dir="/home/jluningham/Projects/bfGWAS_SS"
geno_dir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments

module load R
Rscript make_exp_files2.R

pheno=/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/${gene}_exp_dat.txt
Rscript split_test_train.R ${pheno} ${gene}
pheno2=/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/${gene}_exp_dat_train.txt

done < ${geneDir}/geneList_tmp.txt

while read LINE ; do

gene=${LINE}
echo ${gene}
bfGWAS_SS_dir="/home/jluningham/Projects/bfGWAS_SS"
geno_dir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments

pheno=/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/${gene}_exp_dat_train.txt

filehead=${geno_dir}/ROSMAP_GWAS_Seg_filehead2.txt
#filehead=${geno_dir}/ROSMAP_GWAS_Seg_filehead_short.txt

######### First generate LD files
# LD coefficient directory
# mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_LDs
LD_dir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments/rosmap_LDs
# Score statistics directory
mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_scores
Score_dir=/mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_scores
LDwindow=1000000
study=${gene}

#Submit Array jobs with GetRefLD.sh
qsub -j y -wd ${Score_dir} -N GetRefLD_${study} -t 1-1703 -tc 30 /home/jluningham/Projects/bfGWAS_SS/bin/GetRefLD.sh ${geno_dir} ${pheno} ${filehead} ${LD_dir} ${Score_dir} ${LDwindow}

done < ${geneDir}/geneList_tmp.txt

find . -name \GetRefLD_*.o* -type f -delete
find . -type f \
\( -name 'GetRefLD_*.o*' \) -delete

#commands to run bfgwas_ss_vb on all genes
# initially, I created a list of only genes that had non-zero PrediXcan R2
# I then applied the same pipeline to a second list of the remaining genes
# this is why there are two gene lists used

geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/all_predix_genes.txt"
wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"
j=4

## predixcan gene list
qsub -q b.q -j y -pe smp ${j} -M jluning@emory.edu -m n -t 1-6007 -tc 24 -wd ${wkdir} -N predixGene ${wkdir}/run_makeGW.sh ${geneList} ${j}

## remaining gene list
geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/remaining_genes_list.txt"

# moved working directory away from /home/ to avoid using too much storage
wkdir="/mnt/YangFSS/data/ROSMAP_Expr_BVSRM"
j=4

qsub -q b.q -j y -pe smp ${j} -M jluning@emory.edu -m n -t 1-8292 -tc 24 -wd ${wkdir} -N predixGene ${wkdir}/run_makeGW_compare.sh ${geneList} ${j}


##################################
# after training the prediction model, obtain a list of non-zero PIP SNPs
# (this script was also used to save eQTL summary information and cis/trans proportions)
# make sure directories and gene lists are correct within the R script.
# extract genotypes for the SNPs with non-zero PIP


geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/all_predix_genes.txt"
wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"
cd ${wkdir}
qsub -j y -wd /home/jluningham/Projects/BFGWAS/ROSMAP/Scripts -N mayo_geno2 -pe smp 1 -t 1-6007 /home/jluningham/Projects/BFGWAS/ROSMAP/Scripts/Extract_SNP_VCF_MAYO2.sh ${geneList}

geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/all_predix_genes.txt"

geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/remaining_genes_list.txt"
wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"
cd ${wkdir}
qsub -j y -wd /home/jluningham/Projects/BFGWAS/ROSMAP/Scripts -N ros_geno3 -pe smp 1 -t 11-8292 /home/jluningham/Projects/BFGWAS/ROSMAP/Scripts/Extract_SNP_VCF.sh ${geneList}

## compute GReX
wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"
cd ${wkdir}

qsub -q b.q -M jluning@emory.edu -j y -cwd -N getGrex1 -pe smp 8 Rscript --vanilla compute_grex_predix.R ${NSLOTS}
qsub -q b.q -M jluning@emory.edu -j y -cwd -N getGrex2 -pe smp 8 Rscript --vanilla compute_grex_remain.R ${NSLOTS}
qsub -q b.q -M jluning@emory.edu -j y -cwd -N getGrex3 -pe smp 8 Rscript --vanilla compute_R2_499_mayo.R ${NSLOTS}
qsub -q b.q -M jluning@emory.edu -j y -cwd -N getGrex4 -pe smp 8 Rscript --vanilla compute_R2_499_mayo2.R ${NSLOTS}


## conduct the TWAS

qsub -q b.q -M jluning@emory.edu -j y -cwd -N TWAS -pe smp 8 Rscript --vanilla twas.R ${NSLOTS}
## within this script, it is possible to flag predix gene list or the remaining gene list, or obtain results for both
## computes meta-analyzed effect sizes, saves p-values

## combine TWAS results, compare p-values, obtain manhattan plot

"/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts/twas_result.R"

################################################################

######################################################################################################
######################################################################################################
######################################################################################################
## older pipelines and shell scripts, when first testing VB vs MCMC, cis-only vs genome wide
wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"
j=1
qsub -q b.q -j y -pe smp ${j} -M jluning@emory.edu -m e -t 1-2500 -tc 1 -wd ${wkdir} -N cis_only ${wkdir}/run_makeGW_cis.sh ${geneList} ${j}

geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/all_predix_genes.txt"
#geneList="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/pilot_gene_list.txt"
#geneDir="/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/"
wkdir="/home/jluningham/Projects/BFGWAS/ROSMAP/Scripts"
j=4
nodes="node02.local|node03.local|node04.local|node05.local|node06.local|node07.local|node08.local|node09.local"

## normally 400, but start with 5 test genes
qsub -q b.q -j y -pe smp ${j} -M jluning@emory.edu -m e -t 1600-2100 -tc 30 -l h=${nodes} -wd ${wkdir} -N predix_genes ${wkdir}/run_makeGW_vb.sh ${geneList} ${j}
