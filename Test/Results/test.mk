.DELETE_ON_ERROR:

all: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/pre_em.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.1.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param1.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R1.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.2.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param2.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R2.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.3.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param3.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R3.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.4.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param4.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R4.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.5.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param5.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R5.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/pre_em.OK: 
	rm -f -r /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT
	mkdir -p /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT
	cp -f /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current
	> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt
	> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/pre_em.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.0.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/pre_em.OK
	/home/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Test cis_0.5_He2_0.2_trainExpr_11 499 0.996475 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current 1040101 1065571 19 1e6 
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.0.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param0.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.0.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep paramtemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/paramtemp0.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp0.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/log0.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param0.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/R0.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param0.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/pre_em.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp0.txt 0 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R0.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.1.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R0.OK
	/home/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Test cis_0.5_He2_0.2_trainExpr_11 499 0.996475 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.1.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param1.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.1.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/paramtemp1.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp1.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/log1.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param1.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/R1.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param1.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp1.txt 1 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R1.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.2.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R1.OK
	/home/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Test cis_0.5_He2_0.2_trainExpr_11 499 0.996475 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.2.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param2.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.2.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/paramtemp2.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp2.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/log2.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param2.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/R2.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param2.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp2.txt 2 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R2.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.3.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R2.OK
	/home/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Test cis_0.5_He2_0.2_trainExpr_11 499 0.996475 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.3.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param3.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.3.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/paramtemp3.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp3.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/log3.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param3.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/R3.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param3.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp3.txt 3 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R3.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.4.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R3.OK
	/home/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Test cis_0.5_He2_0.2_trainExpr_11 499 0.996475 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.4.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param4.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.4.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/paramtemp4.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp4.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/log4.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param4.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/R4.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param4.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp4.txt 4 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R4.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.5.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R4.OK
	/home/jchen/BayesianTWAS/bin/run_Estep.sh /mnt/YangFSS/data/jchen/BayesianTWAS/Test cis_0.5_He2_0.2_trainExpr_11 499 0.996475 10000 10000 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/LD /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Score_stat /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current 1040101 1065571 19 1e6
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.5.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param5.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/cis_0.5_He2_0.2_trainExpr_11.5.OK 
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep paramtemp | sort ` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/paramtemp5.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep hyptemp | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp5.txt
	cat `ls -d -1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/output/** | grep log | sort` > /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/log5.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param5.OK

/mnt/YangFSS/data/jchen/BayesianTWAS/Test/R5.OK: /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/cp_param5.OK
	Rscript --no-save --no-restore --verbose /mnt/icebreaker/data2/home/jluningham/Projects/bfGWAS_SS/bin/Mstep.R /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/hyptemp5.txt 5 1e-6 0.1 /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/EM_result.txt /mnt/YangFSS/data/jchen/BayesianTWAS/Test/hypval.current MCMC >> /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Rout.txt
	touch /mnt/YangFSS/data/jchen/BayesianTWAS/Test/R5.OK

clean_err: 
	-rm -rf /mnt/YangFSS/data/jchen/BayesianTWAS/Test/slurm_err/*.err
clean: 
	-rm -rf /mnt/YangFSS/data/jchen/BayesianTWAS/Test/*.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/Eoutput/*.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/OUT/*.OK /mnt/YangFSS/data/jchen/BayesianTWAS/Test/slurm_err/*.err
