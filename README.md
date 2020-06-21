# Bayesian Genome-wide (BGW) TWAS Software Usage
### Justin Luningham
### 4/04/2020

This manual describes the functionality and usage of software for a novel Bayesian Transcriptome-wide Study. This methodology is equipped to include potential eQTL outside of the $cis$ region for a given gene when training the model to predict genetically regulated gene expression (GReX). 

The method is implemented through several command-line argument described in this manual. Note that the genome-wide Bayesian method used for prediction is further described in the supplement for Yang et al. (https://github.com/yjingj/bfGWAS/blob/master/bfGWAS_Manual.pdf). 

# Required Files

The following information is required, stored in various text files with pre-defined formats: a gene expression file, genotype files for the training sample, a list of filenames for training genotype files, genotype files for the prediction sample, a list of filenames for prediction genotype files, and phenotype files for the prediction sample.

## Gene Expression File 

A file containing a list of genes and gene expression scores for the training sample. Each gene should be listed in a row. The first five columns need to contain the chromosome, starting position, ending position, Gene ID, and Gene Name. Columns from 6 to $N_{train}$ should contain the subject IDs for each training sample, with gene expression scores for each gene. An example of the first row and 6 columns of such a file is below: 

| CHROM |	GeneStart |	GeneEnd |	    TargetID     | GeneName |	 ROS20275399 |
| ----- | --------- | ------- | ---------------- | -------- | ------------ |
|  19	  |  1040101	| 1065571 |	 ENSG00000064687 |	 ABCA7  | 0.6707739044 |


## Genotype Files - Training Sample

Genotype files for both training and prediction sets are required for analysis. The genotype files should be a VCF text file zipped by `gzip` (file names ending with .vcf.gz), e.g., `chr1_seg1.vcf.gz`
(http://samtools.github.io/hts-specs/VCFv4.1.pdf).

Note that the training software can handle `plink` PED, bim, and fam files or genotype text file in a specific format to carry out summary statistics for GReX training. However, the process currently requires VCF files for the prediction dataset, and the scripts as constructed utilize VCF for training. 

## Genotype Files - Prediction Sample

The genotype files for the prediction sample should be stored by chromosome, as is common for GWAS samples. The files should be VCF files, zipped with `gzip`, and the file name should include `CHR${chr}` for each chromosome. 


## Genotype file lists 

For the training procedure, users can split the genotype information into many segments based on LD information, each segment containing approximately 3,000-5,000 SNPs. To facilitate computation, users should provide a text file that is a list names for the genome segment files. These should contain a common name (e.g., the name of study or sample), plus the chromosome, start position, and end position of the genome block - all seperated by underscore `_`. E.g.:

|            #filename             |
| -------------------------------- | 
| Rosmap_GWAS_19_610729_2098396    | 
| Rosmap_GWAS_15_38530777_40384132 |


The prediction genotype file list should contain the names of the files for each chromosome. This is useful in the event that files include the sample or study names. E.g.: 

|                #filename                |
| --------------------------------------- | 
| MayoLOADGWAS_CHR1_filtered.dose.geno.gz |
| MayoLOADGWAS_CHR2_filtered.dose.geno.gz |
| MayoLOADGWAS_CHR3_filtered.dose.geno.gz |
| MayoLOADGWAS_CHR4_filtered.dose.geno.gz |
                                

Note that the file extension suffix `.vcf.gz` should not be listed with the file names. 

## Phenotype File 

A phenotype file is required for the prediction dataset. This file contains only two columns: the first column is a list of IDs that match the IDs in the prediction VCF files, and the second column is a  quantitative phenotype score. Column names should not be included. 

For creating the $\hat{GReX}$ values, the actual phenotypes scores are not used, but this file is needed to match phenotype IDs with genotype IDs. Therefore, the second column does not have to be a meaningful variable and can be a placeholder of some kind. 

## Directory Paths

The entire process requires the storage of many results files for each gene (storage is optionally temporary). As the tool is only feasible in a high-performance computing environment, users should supply the following paths: 

- `training_genotype` directory: a directory where all genotype segment files are stored for the training sample
- `prediction_genotype` directory: a directory where all genotype files are stored for the prediction sample (by chromosome)
- `results_directory`: a location where summary stats and training results can be stored 
- `LD_directory`: a directory to store LD information corresponding to genome blocks. This step only needs to be completed once, and the LD information among the genotypes can then be used for all subsequent genes 
- `Scripts_directory`: file path to the "Scripts" folder downloaded from this site, containing the software and scripts needed to implement the Bayesian $\hat{GReX}$ method 


## Other necessary info: 

- Total number of segments in the training genotype files
- Sample size of training sample
- $p$-value threshold for including genome segments in training model 
- maximum number of segments to include in training model 
- Genotype format (`GT` for ${0,1,2}$ Genotypes or `DS` for $[0-2]$ Dosages) for training and for prediction
- Name of `tabix` module available in your computing environment
- Number of available cores for parallel computing 

## Example

```
geneFile=/Example/Data/Gene_Exp_combination.txt
gene=ABCA7
Scripts_dir=/BayesianTWAS #current directory
Res_dir=/Example/ROSMAP_Expr_BVSRM
geno_dir=/Example/Data/Geno #we were not about to share the original genotype
LDdir=/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs
Genome_Seg_File=/Example/Data/ROSMAP_GWAS_Seg_filehead.txt

```

# Usage 

A simple command executing a shell script for each step is detailed below, along with corresponding necessary information. 


## Step 1: Obtain Summary Statistics 

The first shell script will run single-variant TWAS and save summary statistics in the necessary format. This shell script requires 9 arguments: 1) `${gene}` 2) `${geneFile}` 3) `${geno_dir}` 4) `${Scripts_dir}` 5) `${Res_dir}` 6) `${LDdir}` 7) `${Genome_Seg_File}` 8) `${number_segments}` 9) `${number_cores}`


```
################################################################
################################################################
### Step 1: obtain summary statistics (aka Score Statistics)
### Runs single-variant GWAS on the training sample genotypes and available expression data.
### first extracts gene location information for target gene from geneFile.
################################################################
################################################################

num_segments=2
num_cores=4

${Scripts_dir}/Step1.sh ${gene} ${geneFile} ${geno_dir} ${Scripts_dir} ${Res_dir} 
${LDdir} ${Genome_Seg_File} ${num_segments} ${num_cores}

```


This shell script will create a directory called `${gene}_scores` within the `${Res_dir}`, where results for each genome segment will be stored. 

**Note**: for the first gene, this script will also compute LD information and save the information per genome segment in the `${LDdir}`. This step is time consuming, but only needs to be completed one time. 

## Step 2: Prune genome segments 

Step 2 reduces the number of genome segments considered for the Bayesian training model for $GReX$. Unique arguments required for this shell script are the $p$-value threshold for inclusion and the maximum number of segments to include. The arguments are 1) `${gene}` 2) `${geneFile}` 3) `${Res_dir}` 4) `${p_thresh}` 5) `${max_blocks}`

```
################################################################
################################################################
### Step 2: Prune blocks
### Select a limited number of genome segments to consider in the EM training model
### cis blocks are always included (when available), then blocks are filtered
# by p-value threshold, ranked by smallest p-value, and set at a max number of blocks.
################################################################
################################################################

p_thresh=0.00000005
max_blocks=100

${Scripts_dir}/Step2.sh ${gene} ${geneFile} ${Res_dir} ${p_thresh} ${max_blocks}


filehead=${Res_dir}/${gene}_scores/${gene}_signif_segments.txt

```

This shell script creates the file `${gene}_signif_segments.txt`, which needs to be passed to the commands for Step 3. 

## Step 3: Execute Bayesian Training model on selected genome segments 

Step 3 will use summary statistics to carry out Bayesian variable selection regression with the EM-MCMC algorithm. This process is limited to the pre-selected genome segments from Step 2, rather than the entire genome. The sample size is the only unique variable that needs to be created for this script. The arguments passed to the script are 1) `${gene}` 2) `${geneFile}` 3) `${geno_dir}` 4) `${Scripts_dir}` 5) `${Res_dir}` 6) `${LDdir}` 7) `${N}` 8) `${num_cores}`

```

################################################################
################################################################
### Step 3: EM-MCMC training
### use Make with Perl script to create file
################################################################
################################################################


N=499

${Scripts_dir}/Step3.sh ${gene} ${geneFile} ${geno_dir} ${Scripts_dir} ${Res_dir} 
${LDdir} ${N} ${num_cores}
```

The final iteration of the EM-MCMC algorithm will result in a file that lists all variants found to have  non-zero effect sizes (specifically, a non-zero posterior probability of inclusion, $PP_i$). The chromosome, position, rsID, reference and alternative alleles, $maf$, $cis$ or $trans$ status relative to the gene, $PP_i$, and effect size $w_i$ will all be saved in the file `${Res_dir}/${gene}_TWAS/Eoutput/grexparam.txt`.

**note**: within this script, there are several optional arguments that control parameters of the EM-MCMC algorithm, such as the number of EM iterations and the number of MCMC draws. These are set to 5 EM iterations with 10,000 burn-in and 10,000 additional MCMC iterations. If desired, these can be modified by the user within the script. 

## Step 4: Extract Genotypes from the Prediction Sample 

In step 4, the variants found to have non-zero effect sizes are matched to the genotype information available in the prediction sample and extracted into a text file. The required arguments for this step are 1) `${gene}` 2) `${Res_dir}` 3) `${pred_geno_dir}` 4) `${pred_geno_filenames}` 5) `${pred_pheno_file}` 6) `${genotype_format}` 7) `${tabix_mod}`


```
################################################################
################################################################
### Step 4: Extract genotypes for prediction sample that
### align with eQTL with non-zero effect size from the bayesian training result
################################################################
################################################################

pred_geno_dir=/mnt/YangFSS/data/AMP-AD/Mayo/Genotype/Impute2_1KG
pred_geno_filenames=/Example/pred_geno_list.txt
pred_pheno_file=/mnt/YangFSS/data/AMP-AD/Mayo/Phenotype/MayoPhenoAD.txt
genotype_format=GT #or DS
tabix_mod=tabix/0.2.6

${Scripts_dir}/Step4.sh ${gene} ${Res_dir} ${pred_geno_dir} ${pred_geno_filenames} 
${pred_pheno_file} ${genotype_format} ${tabix_mod}

```

This step will produce a file with genotype information for all of the prediction samples for those matching variants that were found to have non-zero effect size in the training step. The file is called `${gene}_grex_genotypes.geno`.

## Step 5

Step 5 uses the prediction sample genotypes and effect sizes from the training model to compute $\hat{GReX}$. The arguments required are 1) `${gene}` 2) `${Scripts_dir}` 3) `${Res_dir}`

```
################################################################
################################################################
### Step 5: Obtain predicted GREX
################################################################
################################################################

${Scripts_dir}/Step5.sh ${gene} ${Scripts_dir} ${Res_dir}
```

This script provides the two primary results of ultimate interest from the method. First, for each individual in the prediction sample, $\hat{GReX}$ for person $j$ is calculated as 

$$\hat{GReX}_j = \sum_{i=1}^p x_{ij} \hat{PP}_i \hat{w}_i$$ 

where $x_{ij}$ is the genotype for person $j$ on variant $i$, $\hat{PP}_i$ is the estimated posterior probability for variant $i$ from the training model, and $\hat{w}_i$ is the estimated regression effect size for variant $i$ from the training model. These are stored for each individual in a file called `${Res_dir}/${gene}_GREX/${gene},_predicted_GREX.txt`. 

Secondly, the sums of the posterior probabilities are also provided. These represent the expected number of eQTL for a given gene. Users are provided the total expected eQTL as well as the expected *cis*- and *trans*-eQTL. These are provided in the file `${Res_dir}/${gene}_GREX/${gene}_PPis.txt`.

# Additional Details 

The current Bayesian TWAS is an extension of existing Bayesian Functional GWAS. Instead of requiring a functional annotation file, the current method automatically assigns $cis$ or $trans$ annotation based on the location information of the outcome gene. Numerous arguments can be used to modify the EM-MCMC algorithm in Step 3, but should be done with caution. These arguments are detailed in Yang et al. 2017 (https://github.com/yjingj/bfGWAS/blob/master/bfGWAS_Manual.pdf). 


