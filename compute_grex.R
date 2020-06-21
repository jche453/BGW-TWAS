options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)

gene<-args[1]
wkdir<-args[2]


#print(genes)
###### Analyze BFGWAS results

library(data.table)
library(dplyr)

setwd(paste0(wkdir))

geno_pred<- fread(paste0(gene, "_grex_genotypes.geno"), header=T)

gwas_res <- fread(paste0("Eoutput/grexparam.txt" ), header=F)

names(gwas_res)<- c("CHROM", "POS", "ID", "REF", "ALT", "maf", "location", "pi", "beta")

tmp_test<-(merge(gwas_res, geno_pred))

X_geno <- data.matrix(t(tmp_test[,10:length(tmp_test)]))


beta<-tmp_test$beta
pi <- tmp_test$pi
weight<-beta*pi

GREX_hat <- X %*% weight


GREX_hat<-data.frame(GREX_hat=GREX_hat)
GREX_hat$ID <- row.names(GREX_hat)

write.table(GREX_hat, file=paste0(gene,"_predicted_GREX.txt"))

cis_PPi<-gwas_res %>% filter(location==0) %>% summarize(sum(pi))
trans_PPi<-gwas_res %>% filter(location==1) %>% summarize(sum(pi))
total_PPi<-cis_PPi+trans_PPi

PPi<-data.frame(cbind(Gene=paste0(gene), total_PPi=total_PPi, cis_PPi=cis_PPi, trans_PPi=trans_PPi))

write.table(PPi, file=paste0(gene, "_PPis.txt", col.names=T, row.names=F, quote=F))

