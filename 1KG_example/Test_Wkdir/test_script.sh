
# Script for testing bfGWAS and bfGWAS_SS

################ setup directories
bfGWAS_SS_dir="/home/jyang/GIT/bfGWAS_SS"
data_dir="/home/jyang/GIT/bfGWAS_SS/1KG_example/ExData"

######### Example script for running BFGWAS_SS 
wkdir="/home/jyang/GIT/bfGWAS_SS/1KG_example/Test_Wkdir"
# makdir -p $wkdir
cd $wkdir

######### Run BFGWAS_SS with a single genome block (08/21/2018)
filehead=CFH_REGION_1KG

## Run BFGWAS_SS with individual-level data
${bfGWAS_SS_dir}/bin/Estep_mcmc \
-vcf ${data_dir}/vcfs/${filehead}.vcf.gz \
-p ${data_dir}/phenoAMD_1KG.txt \
-a ${data_dir}/annos/Anno_CFH_REGION_1KG.gz \
-fcode ${data_dir}/AnnoCode6.txt \
-hfile ${wkdir}/hypval.current \
-GTfield GT \
-bvsrm -maf 0.005 -win 100 -smin 0 -smax 10 -w 10000 -s 10000 \
-o ${filehead}_indv -initype 3 -seed 2017 \
-LDwindow 1000000  -saveSS -zipSS  > ${filehead}_indv.output.txt &

### Run BFGWAS_SS with summary level file
# Need the MAF of analyzed SNPs in the summary data
${bfGWAS_SS_dir}/bin/Estep_mcmc \
-score ${wkdir}/output/${filehead}_indv.score.txt.gz \
-LDcorr ${wkdir}/output/${filehead}_indv.LDcorr.txt.gz \
-a ${data_dir}/annos/Anno_CFH_REGION_1KG.gz \
-fcode ${data_dir}/AnnoCode6.txt \
-n 2504 -pv 0.25 -inputSS \
-hfile ${wkdir}/hypval.current \
-bvsrm -maf 0.005 -win 100 -smin 0 -smax 10 -w 10000  -s 10000  \
-o ${filehead}_ss -initype 3 -seed 2017 > ${filehead}_ss.output.txt  &

### example R code for analyzing the results
see "/1KG_example/AnalyzeResults/Analysis.r" for details of loading data and make plots


############## Run BFGWAS_SS with multiple genome block
##### To be revised (08/21/2018)
##### perl script is used to generate makefile #############

${bfGWAS_dir}/bin/gen_mkf.pl \
-w ${wkdir} \
--Estep ${bfGWAS_dir}/bin/Estep_mcmc \
--ad ${data_dir}/annos \
--geno vcf --ac ${data_dir}/AnnoCode6.txt  \
--gd ${data_dir}/vcfs \
--pheno ${data_dir}/phenoAMD_1KG.txt \
--hyp ${data_dir}/InitValues6.txt \
-f ${data_dir}/fileheads_4region.txt \
--rs ${bfGWAS_dir}/bin/Mstep.r \
-G GT --maf 0.005 --smax 10 -b 10000 -N 10000 --NL 10000 \
--em 5 -j testjob -l slurm --part nomosix,main,assembly --mem 3000 \
--mf ${wkdir}/Test_bfGWAS.mk

## Run Makefile with 4 parallel jobs
make -k -C ${wkdir} -f ${wkdir}/Test_bfGWAS.mk -j 10 > ${wkdir}/make.output 2> ${wkdir}/make.err  &

## Clean all jobs when you need rerun everything
# make -f ${wkdir}/Test_bfGWAS.mk clean


#################
############## Output under ${wkdir} ###################

######### /OUT/ : saves all screen outputs

######### /output/ : saves all MCMC output from E-step / will be overridden by the next E-step
# filehead.log.txt contains log file for MCMC 
# filehead.paramtemp contains estimates for each variant with columns "ID", "chr", "bp", "ref", "alt", "maf", "func", "beta", "pi", "Zscore", "SE_beta", "LRT", "pval_LRT", "rank";

# "ID" : variant ID
# "chr" : chromosome number 
# "bp" : base pair position
# "REF" : reference allel
# "ALT" : alternative allel
# "maf" : MAF of the variant
# "func" : annotation code used in --ac FuncAnno4.txt
# "beta" : effect size estimate
# "pi" : causal probability for each variant
# "Zscore" : Zscore by single likelihood ratio test
# "SE_beta" : from single likelihood ratio test
# "LRT" : test statistic by single likelihood ratio test
# "pval_LRT" : pvalue by single likelihood ratio test
# "rank" : rank within block based on p-values, 0:top significant variant by pvalue

# function LoadEMdata() in bin/R_funcs.r can be used to read this paramtemp file into R, requiring library "data.table" and "ggplot2"

# filehead.hyptemp contains estimates required for M-step 
# filehead.mcmc contains all included variants (id:chr:pos:ref:alt, seperated by ";") in the M-step, one row per MCMC iteration (can be used for calculating regional posterior inclusion probabilities)


######### /slurm_err/ : saves all error file for slurm jobs


###### under folder /Eoutput/ : main results ########
# let i denote the ith EM iteration
# log${i}.txt contains all log files for all blocks
# hyptemp${i}.txt contains all hyptemp files for all blocks
# paramtemp${i}.txt contains all paramtemp files for all blocks
# EM_result.txt contains all hyper parameter estimates with columns "EM_iter_#", "PVE/heritability", "likelihood", every 4 of the following columns denotes the "causal probability" "causal probability SE" "effect size variance" "effect size variance SE" for group 0, 1, 2, ...
# R function LoadEMhyp() in bin/R_funcs.r can be used to read EM_result.txt file
# R function CItable() in bin/R_funcs.r can be used to convert one row of EM_result.txt to a table of annotations


############ example R code for analysis
see "/1KG_example/AnalyzeResults/Analysis.r" for details of loading data and make plots


## Run BFGWAS individual-level data
bfGWAS_dir="/home/jyang/GIT/bfGWAS"

${bfGWAS_dir}/bin/Estep_mcmc \
-vcf ${data_dir}/vcfs/${filehead}.vcf.gz \
-p ${data_dir}/phenoAMD_1KG.txt \
-a ${bfGWAS_dir}/1KG_example/ExData/annos/Anno_CFH_REGION_1KG.gz \
-fcode ${bfGWAS_dir}/1KG_example/ExData/AnnoCode6.txt \
-hfile ${wkdir}/hypval.current \
-GTfield GT -rmin 1 -rmax 1 -rv 1 \
-bvsrm -maf 0.005 -win 100 -smin 0 -smax 10 -w 10000 -s 10000 \
-o ${filehead}_bfgwas -initype 3 -seed 2017 \
> ${filehead}_bfgwas.output.txt &