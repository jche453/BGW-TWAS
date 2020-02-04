
setwd("/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files")

library(data.table)
library(dplyr)
library(parallel)

gene_list <- readLines("remaining_genes_list.txt") ##list of gene names
gene_exp <- fread("Gene_Exp_combination.txt", header=T)  ##gene expression data
target_dir<-"/mnt/YangFSS/data/ROSMAP_Expr_BVSRM/exp_dat"

get_exp_all<-function(index){tryCatch({
gene<-gene_list[index]

dat <- gene_exp %>% filter(GeneName == paste0(gene))

dat2 <- data.frame(exp.dat = t(dat[,-c(1:5)]))
dat2$ID <- rownames(dat2)

dat2 <- dat2 %>% select(ID, exp.dat)

fwrite(dat2, file = paste0(target_dir, "/",gene,"_exp_dat.txt"), quote=F, col.names=F, row.names=F, sep="\t")

pv <- var(dat2$exp.dat)
dat3 <- cbind(dat[,1:5], pv=pv)
write.table(x=dat3, file = paste0(target_dir, "/",gene,"_info.txt"), quote=F, col.names=F, row.names=F )
},error=function(e) NULL)}

mclapply(1:length(gene_list), function(d) get_exp_all(index=d),mc.cores=1)

