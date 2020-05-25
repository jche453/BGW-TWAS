.DELETE_ON_ERROR:

all: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/pre_em.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.1.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param1.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R1.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.2.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param2.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R2.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.3.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param3.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R3.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.4.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param4.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R4.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.5.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param5.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R5.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/pre_em.OK: 
	rm -f -r /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT
	mkdir -p /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT
	cp -f /mnt/YangFSS/data/jchen/BayesianTWAS/hypval.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current
	> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt
	> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/pre_em.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.0.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/pre_em.OK
	/mnt/YangFSS/data/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX Rosmap_GWAS_19_610729_2098396 499 1.02401 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_scores /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current 1040101 1065571 19 1e6 
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.0.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param0.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.0.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep paramtemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/paramtemp0.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp0.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/log0.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param0.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R0.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/pre_em.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp0.txt 0 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R0.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.1.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R0.OK
	/mnt/YangFSS/data/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX Rosmap_GWAS_19_610729_2098396 499 1.02401 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_scores /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.1.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param1.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.1.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/paramtemp1.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp1.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/log1.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param1.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R1.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param1.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp1.txt 1 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R1.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.2.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R1.OK
	/mnt/YangFSS/data/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX Rosmap_GWAS_19_610729_2098396 499 1.02401 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_scores /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.2.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param2.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.2.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/paramtemp2.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp2.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/log2.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param2.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R2.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param2.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp2.txt 2 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R2.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.3.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R2.OK
	/mnt/YangFSS/data/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX Rosmap_GWAS_19_610729_2098396 499 1.02401 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_scores /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.3.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param3.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.3.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/paramtemp3.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp3.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/log3.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param3.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R3.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param3.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp3.txt 3 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R3.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.4.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R3.OK
	/mnt/YangFSS/data/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX Rosmap_GWAS_19_610729_2098396 499 1.02401 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_scores /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.4.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param4.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.4.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/paramtemp4.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp4.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/log4.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param4.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R4.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param4.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp4.txt 4 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R4.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.5.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R4.OK
	/mnt/YangFSS/data/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX Rosmap_GWAS_19_610729_2098396 499 1.02401 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/Data/ROSMAP_GWAS_Segments/rosmap_LDs /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_scores /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.5.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param5.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/Rosmap_GWAS_19_610729_2098396.5.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/paramtemp5.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp5.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/log5.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param5.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R5.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/cp_param5.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/hyptemp5.txt 5 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/R5.OK

clean_err: 
	-rm -rf /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/slurm_err/*.err
clean: 
	-rm -rf /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/*.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/Eoutput/*.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/OUT/*.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Example/ROSMAP_Expr_BVSRM/ABCA7_GREX/slurm_err/*.err