#!/bin/sh

# create makefiles and run bfgwas for a list of genes in parallel
# provide the submission script the list of genes in a text file

### first step: set up score statistics for each gene
#!/bin/sh

cd /mnt/YangFSS/data/jchen/BayesianTWAS/Test

module load tabix
#directory of the tool
bfGWAS_SS_dir=/home/jchen/BayesianTWAS

# Set up a LD and Score statistics directory
# mkdir /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD
LDdir=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD
# mkdir /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat
Score_dir=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat
#Set LDwindow = 1 if LD files have already been generated, if not, set LDwindow = 1000000
#LDwindow=1
LDwindow=1000000

# -g With input genotype file in dosage format, calculate LD and score statistics
# Use -vcf if genotype file is in vcf format
# -p specify phenotype file
# -o output name
${bfGWAS_SS_dir}/bin/Estep_mcmc -g /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Data/simulation_dosages_all.geno \
                                -p /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Data/cis_0.5_He2_0.2_trainExpr_11.txt \
                                -o cis_0.5_He2_0.2_trainExpr_11 -LDwindow ${LDwindow} -GTfield DS -saveSS -zipSS

if [ ! -f ${LDdir}/cis_0.5_He2_0.2_trainExpr_11.LDcorr.txt.gz ] ; then
  echo Copy log.txt and LDcorr.txt.gz back to $LDdir
  cp -f  output/cis_0.5_He2_0.2_trainExpr_11.LDcorr.txt.gz ${LDdir}/
  cp -f  output/cis_0.5_He2_0.2_trainExpr_11.log.txt ${LDdir}/
else
  echo Copy log.txt back to $Score_dir
  cp -f  output/cis_0.5_He2_0.2_trainExpr_11.log.txt ${Score_dir}/
fi

echo Copy score.txt back $Score_dir
cp -f  output/cis_0.5_He2_0.2_trainExpr_11.score.txt.gz ${Score_dir}/

rm -rf output

###Second Step: run genome-wide bfGWAS model
#!/bin/sh
bfGWAS_SS_dir=/home/jchen/BayesianTWAS
#variance of phenotype
pv=0.996475

# mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_LDs
LDdir=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD
# Score statistics directory
#mkdir /mnt/YangFSS/data/ROSMAP_GWAS_Segments/${gene}_scores
Score_dir=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat
LDwindow=1

filehead=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Data/filehead.txt

# get train, get test sets
# match IDs to gene_exp_dat
N=499

#directory where the pl and sh scripts were located
EMdir=/home/jchen/BayesianTWAS/bin

# Hyper parameter file (a seperate file from hypval.current under the working directory)
hfile=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Data/hypval.txt

# Specify computation specs
em=5; # EM steps
burnin=10000; # Burn-in iterations in MCMC
Nmcmc=10000; # MCMC iteration number
study=test

wkdir=/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Results

mkfile=${wkdir}/${study}.mk

${bfGWAS_SS_dir}/bin/gen_mkf.pl \
--EMdir ${EMdir} --hyp ${hfile} \
-n ${N} --pv ${pv} \
-w ${wkdir} --geno sumstat \
-f ${filehead} -l local --pheno /mnt/YangFSS/data/jchen/BayesianTWAS/Test/cis_0.5_He2_0.2_trainExpr_11.txt \
--targ 19 --start 1040101 --end 1065571 \
--LDdir ${LDdir} --Scoredir ${Score_dir} \
-j BFGWAS_${study} --em ${em} -b ${burnin} -N ${Nmcmc} \
--mf ${mkfile}

cd ${wkdir}
echo run simulation study
make -f ${mkfile} clean


echo Run make with $mkfile

make -k -C ${wkdir} -f ${mkfile} -j 1 > ${wkdir}/make.output 2> ${wkdir}/make.err

exit
