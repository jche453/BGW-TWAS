# rm(list=ls(all=TRUE))

####### Source Util Functions and set data directories
source("/home/jyang/GIT/bfGWAS_SS/bin/R_funcs.r")

DataDir = "/home/jyang/GIT/bfGWAS_SS/1KG_example/ExData/" # example data directory
OutDir = "/home/jyang/GIT/bfGWAS_SS/1KG_example/AnalyzeResults/" # directory to save plots
ResultDir = "/home/jyang/GIT/bfGWAS_SS/1KG_example/Test_Wkdir" # result directory

setwd(ResultDir)

######## Compare results
paramdata_bfgwas = LoadEMdata_bfgwas(filename="./output/CFH_REGION_1KG_bfgwas.paramtemp")
# head(paramdata_bfgwas)

paramdata_indv = LoadEMdata(filename="./output/CFH_REGION_1KG_indv.paramtemp")
# head(paramdata_indv)

paramdata_ss = LoadEMdata(filename="./output/CFH_REGION_1KG_ss.paramtemp")
# head(paramdata_ss)

##
qplot(paramdata_ss$beta, paramdata_indv$beta) 
ggsave("comp_beta_ss.pdf")

qplot(paramdata_bfgwas$beta, paramdata_indv[paramdata_bfgwas$ID, ]$beta) 
ggsave("comp_beta_bfgwas.pdf")

##
sum(paramdata_bfgwas$pi)
sum(paramdata_indv$pi)
sum(paramdata_ss$pi)

##
paramdata_bfgwas[paramdata_bfgwas$pi > 0.1, ]
# paramdata_indv[paramdata_bfgwas$ID[paramdata_bfgwas$pi > 0.1], ]
paramdata_indv[paramdata_indv$pi > 0.1, ]
paramdata_ss[paramdata_ss$pi > 0.1, ]

## True causal variants 
VCF_Data = fread(paste(DataDir, "vcfs/causalSNP.vcf", sep = ""), sep="\t", header=TRUE)
Geno = VCF_Data[, 10:2513] # genotype of true causal SNPs
Geno[Geno == "0|0"] = 0
Geno[Geno == "1|1"] = 2
Geno[Geno == "1|0" | Geno == "0|1"] = 1

Geno = apply(Geno, 1, as.numeric)
colnames(Geno) = VCF_Data$ID

maf_df = data.frame(VCF_Data[, 1:5],  maf = apply(Geno, 2, calcMAF))
maf_df = maf_df[maf_df$maf > 0.005, ]
maf_df[1:3, ]

# Variants that were excluded from the analysis will show NA values
paramdata_bfgwas[maf_df$ID[1:3], ]
paramdata_indv[maf_df$ID[1:3], ]
paramdata_ss[maf_df$ID[1:3], ]

## Load Phenotype
pheno = read.table("/home/jyang/GIT/bfGWAS/1KG_example/ExData/phenoAMD_1KG.txt", header = FALSE)
y = pheno$V2
names(y) = pheno$V1

x = Geno[, 1]
names(x) = colnames(VCF_Data)[-(1:9)]

x1 = x - mean(x)
y1 = y - mean(y)

sum(y1 * x1) / sum(x1 * x1)
sqrt( var(y1 )  * sum(x1 * x1))

summary(lm(y[names(x)] ~ x))

################################################
######## Load and analyze GWAS results by bfGWAS
paramdata = LoadEMdata(filename="./Eoutput/paramtemp5.txt")
head(paramdata)

#### Variants with association probabilities > 0.1068, potential causals
print(paramdata[paramdata$pi>0.1068, ])

#### Variants with SVT pvalue < 5e-8, by likelihood ratio tests
paramdata[paramdata$pval_LRT<5e-8, ]

#### True variants used to simulate this dataset
VCF_Data = fread(paste(DataDir, "vcfs/causalSNP.vcf", sep = ""), sep="\t", header=TRUE)
Causal_SNP_Result = paramdata[VCF_Data$ID, ]
# variants that were excluded from the analysis will show NA values
print(Causal_SNP_Result[, c(1:3, 6:10, 13), with = FALSE])

###### REMARK: you might want to check if the estimated variants with high association probability are in high LD with the true causal ones

###### Although we are not able to identify all causal variants, we can estimate the probability that there exist at least one causal variant per genome block
system("cat \`ls output/** | grep log\` | grep Region_PIP | cut -d\" \" -f4 > Eoutput/Region_PP.txt")
Region_PP = scan("Eoutput/Region_PP.txt", what = double())
sum(Region_PP > 0.95)

######## Load results of hyper parameters ####### 
test_hyp <- LoadEMhyp(filename = "./Eoutput/EM_result.txt")

# the labels should be consistant to 0,1,2... defined in the annotation code file
group_labels <- as.factor(c("Coding", "UTR", "Promoter", "DHS", "Intronic", "Intergenic"))

test_CI_table <- CItable(test_hyp[6, ], n_type = 6, alpha = 0.95, 
		funcgroup = group_labels)

# Set the values of these categories without associations at NAs, change prior_pp value accordingly
prior_pp = 1e-6
test_CI_table [test_CI_table$pi == prior_pp, 1:6] <- NA 
 
# plot causal proportion estimates, requiring library "ggplot2"
PlotCI_groupPP(hyp_table=test_CI_table, pdfname=paste(OutDir, "test_groupPP.pdf", sep =""), size = 18, tit = "")

# plot effect size estimates
PlotCI_groupESvar(hyp_table=test_CI_table, pdfname=paste(OutDir, "test_ESvar.pdf", sep =""), size = 18, tit = "")


######## Compare hyper estimates to the genome-wide averages
n_group = 6
pp_cols = (1:n_group) * 4
group_pp = as.vector( test_hyp[6, (pp_cols)])
no_asso_groups = (group_pp == prior_pp)
group_pp[no_asso_groups] = 0
pp_se <- as.vector( test_hyp[6, (pp_cols + 1)])
pp_se[no_asso_groups] = 0

# as.vector(table(paramdata$func))

n_vec = c(9221, 203, 0, 1354, 0, 5627) # number of variants per annotation, 
# should be consistant to 0,1,2... defined in the annotation code file

ncausal_vec = n_vec * group_pp # number of associations
group_sigma2 <- as.vector( test_hyp[6, (pp_cols + 2)])
group_sigma2[no_asso_groups] = 0
sigma2_se <- as.vector( test_hyp[6, (pp_cols + 3)])
sigma2_se[no_asso_groups] = 0

### Construct comparison results
comp_group_sigma2 <-comp_group_pp <- data.frame(log_ratio_lcl = rep(NA, n_group), 
							log_ratio=rep(NA, n_group), 
							log_ratio_ucl = rep(NA, n_group))

for(i in 1:n_group){
		comp_group_pp[i, ] = comp_groupEst(est = group_pp, est_se = pp_se, n_vec = n_vec, i = i, conf = 0.95)
		comp_group_sigma2[i, ] = comp_groupEst(est = group_sigma2, est_se = sigma2_se, n_vec = ncausal_vec, i = i, conf = 0.95)
}

exp(comp_group_pp[, 2])
exp(comp_group_sigma2[, 2])


### Plot comparison results
PlotRatio(comp_dat = data.frame(exp(comp_group_pp), funcgroup=group_labels), 
	tit ="Regulatory Annotation", 
	pdfname = paste(OutDir, "comp_test_pp.pdf", sep =""), 
	ymode = 1, size = 28, wid = 8) 


PlotRatio(comp_dat = data.frame(exp(comp_group_sigma2), funcgroup=group_labels), 
	tit ="Regulatory Annotation", 
	pdfname = paste(OutDir, "comp_test_sigma2.pdf", sep =""), 
	ymode = 0, size = 28, wid = 8) 

######## END ################










