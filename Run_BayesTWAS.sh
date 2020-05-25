
################################################################
################################################################
# Step 0: Set up directores and files.
# need a name/identifier for each gene.
# need a file that has gene information with following columns:
# CHR, START POS, END POS, GENE Name, Gene ID
# from column 6 on is the subject ID and expression scores for each gene.
# need directory with genotype, and a Results directory to send output.
# need a file that contains a list of filenames for the genotype data.
# For computation speed, recommend segmenting into many small chunks.
# need scripts directory with perl script, shell files, R script.
################################################################
################################################################


geneFile=/Example/Data/Gene_Exp_combination.txt
gene=ABCA7
Scripts_dir=/BayesianTWAS #current directory
Res_dir=/Example/ROSMAP_Expr_BVSRM
geno_dir=/Example/Data/Geno #we were not about to share the original genotype
LDdir=/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs
Genome_Seg_File=/Example/Data/ROSMAP_GWAS_Seg_filehead.txt


################################################################
################################################################
# Step 1: obtain summary statistics (aka Score Statistics)
# Runs single-variant GWAS on the training sample genotypes and available expression data.
# first extracts gene location information for target gene from geneFile.
################################################################
################################################################

num_segments=2
num_cores=4

${Scripts_dir}/Step1.sh ${gene} ${geneFile} ${geno_dir} ${Scripts_dir} ${Res_dir} ${LDdir} ${Genome_Seg_File} ${num_segments} ${num_cores}

################################################################
################################################################
# Step 2: Prune blocks
# Select a limited number of genome segments to consider in the EM training model
# cis blocks are always included (when available), then blocks are filtered
# by p-value threshold, ranked by smallest p-value, and set at a max number of blocks.
################################################################
################################################################

p_thresh=0.00000005
max_blocks=100

${Scripts_dir}/Step2.sh ${gene} ${geneFile} ${Res_dir} ${p_thresh} ${max_blocks}


filehead=${Res_dir}/${gene}_scores/${gene}_signif_segments.txt


################################################################
################################################################
# Step 3: EM-MCMC training
# use Make with Perl script to create file
################################################################
################################################################


N=499

${Scripts_dir}/Step3.sh ${gene} ${geneFile} ${geno_dir} ${Scripts_dir} ${Res_dir} ${LDdir} ${N} ${num_cores}

################################################################
################################################################
# Step 4: Extract genotypes for prediction sample that
# align with eQTL with non-zero effect size from the bayesian training result
################################################################
################################################################


pred_geno_dir=/mnt/YangFSS/data/AMP-AD/Mayo/Genotype/Impute2_1KG
pred_geno_filenames=/mnt/YangFSS/data/AMP-AD/Mayo/Genotype/Impute2_1KG/file_names.txt
pred_pheno_file=/mnt/YangFSS/data/AMP-AD/Mayo/Phenotype/MayoPhenoAD.txt
genotype_format=GT #or DS
tabix_mod=tabix/0.2.6

${Scripts_dir}/Step4.sh ${gene} ${Res_dir} ${pred_geno_dir} ${pred_geno_filenames} ${pred_pheno_file} ${genotype_format} ${tabix_mod}

################################################################
################################################################
# Step 5: Obtain predicted GREX
################################################################
################################################################

${Scripts_dir}/Step5.sh ${gene} ${Scripts_dir} ${Res_dir}


# optional script for removing files
