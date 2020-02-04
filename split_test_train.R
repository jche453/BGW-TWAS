## separate the training IDs

options(stringsAsFactors = FALSE)
args <- as.numeric(commandArgs(trailingOnly = TRUE))
print(args)

train_IDs <- readLines("/home/jyang/Collaborations/WingoMAP/RNAseq/data/sample_iid_train.txt")
geneList <- readLines("/mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/all_predix_genes.txt")

i<-args[1]
gene<-geneList[i]
print(gene)

gene_exp <- read.table(paste0("/mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/",gene,"_exp_dat.txt"), header=F)
colnames(gene_exp) <- c("ID", "exp")

w<- which(gene_exp$ID %in% train_IDs)
train_file <- gene_exp[w,]
test_file<- gene_exp[-w,]

write.table(train_file, file=paste0("/mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/",gene,"_exp_dat_train.txt"), quote=F, col.names=F, row.names=F, sep="\t")
write.table(test_file, file=paste0("/mnt/icebreaker/data2/home/jluningham/Projects/BFGWAS/ROSMAP/expression_files/",gene,"_exp_dat_test.txt"), quote=F, col.names=F, row.names=F, sep="\t")








