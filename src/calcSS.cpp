/*
	Bayesian Functional GWAS with Summary Statistics --- MCMC (bfGWAS_SS:MCMC)
    Copyright (C) 2017  Jingjing Yang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "calcSS.h"


void CALCSS::CopyFromParam (PARAM &cPar) 
{

    zipSS=cPar.zipSS;

    UnCompBufferSize = cPar.UnCompBufferSize;
    CompBuffSizeVec = cPar.CompBuffSizeVec;
    Compress_Flag = cPar.Compress_Flag;

    file_out=cPar.file_out;
        
    ni_total=cPar.ni_total;
    ns_total=cPar.ns_total;
    ni_test=cPar.ni_test;
    ns_test=cPar.ns_test;
    n_type = cPar.n_type;

    LDwindow=cPar.LDwindow;
    
    indicator_idv=cPar.indicator_idv;   
    indicator_snp=cPar.indicator_snp;

    SNPmean = cPar.SNPmean; 
    pheno_mean = cPar.pheno_mean;
    pheno_var = cPar.pheno_var;

    snpInfo=cPar.snpInfo;
    
    return;
}

//calculat LD correlation matrix
void CALCSS::GetSS(uchar **X, gsl_vector *y, vector< vector<double> > &LD, vector<double> &beta,vector<double> &beta_sd){

    cout << "\nStart calculating summary statistics ... \n";

    Gvec.assign(n_type, 0.0); // set 

    // y is centered by cPar.CopyPheno() 
    double yty;
    gsl_blas_ddot(y, y, &yty); // yty is the sum square of total SST
    pheno_var = yty / ((double)(ni_test-1)) ;
    cout << "pheno_var = " << pheno_var << endl;

    // standardize phenotype y
    // gsl_vector_scale(y, 1.0 / sqrt(pheno_var)); 
    gsl_blas_ddot(y, y, &yty); // calculate yty for one more time after standardization
    
    //cout << "create UcharTable ...\n";
    CreateUcharTable(UcharTable);

    // Create a vector of "SNPPOS" structs snp_pos (snpInfo will be cleared)
    CreateSnpPosVec(snp_pos, snpInfo, ns_total, indicator_snp); //ordered by position here
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp); // order snp_pos by chr/bp

    // define used variables 
    gsl_vector *xvec_i = gsl_vector_alloc(ni_test);
    gsl_vector *xvec_j = gsl_vector_alloc(ni_test);
    gsl_vector *xbeta_i = gsl_vector_alloc(ni_test);

    double xtx_ij, xty, xtx_i, beta_i, beta_sd_i;

    LD.clear();
    beta.clear();

    for (size_t i=0; i<ns_test; ++i) {

        LD.push_back(vector<double>());

        //calculate xtx_i
        getGTgslVec(X, xvec_i, snp_pos[i].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
        gsl_blas_ddot(xvec_i, xvec_i, &xtx_i);
        LD[i].push_back(xtx_i);

        for (size_t j=0; j<n_type; j++) {
            if (snp_pos[i].indicator_func[j]) {
                Gvec[j] += (xtx_i / (double)ni_test);
                //Gvec[j] += (xtx );
                continue;
            }
        }

        //calculate effect-size
        gsl_blas_ddot(xvec_i, y, &xty);
        beta_i = xty / xtx_i;
        beta.push_back(beta_i);

        gsl_vector_memcpy(xbeta_i, xvec_i);
        gsl_vector_scale(xbeta_i, -beta_i);
        gsl_vector_add(xbeta_i, y);
        gsl_blas_ddot(xbeta_i, xbeta_i, &beta_sd_i);
        beta_sd_i = sqrt(beta_sd_i / ((double)ni_test * xtx_i));
        beta_sd.push_back(beta_sd_i);

        if(i < (ns_test-1) ){
            //calculate xtx_ij 
            for(size_t j=(i+1); j < ns_test; ++j){
                if( (snp_pos[j].chr == snp_pos[i].chr) && (snp_pos[j].bp <= snp_pos[i].bp + LDwindow) ){
                    getGTgslVec(X, xvec_j, snp_pos[j].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
                    gsl_blas_ddot(xvec_i, xvec_j, &xtx_ij);
                    LD[i].push_back(xtx_ij);
                }else{break;}
            }
        }

    }

    gsl_vector_free(xvec_i);
    gsl_vector_free(xvec_j);
    gsl_vector_free(xbeta_i);

    return;
}


void CALCSS::WriteSS(const vector< vector<double> > &LD, const vector<double> &beta, const vector<double> &beta_sd)
{
    cout << "\nStart writing summary statistics ... \n";
    String fout = file_out.c_str();

    String LD_file_str = "./output/" + fout;
    String beta_file_str = "./output/" + fout;

    IFILE LDout=NULL;
    IFILE Beta_out=NULL;

    if(zipSS){
        LD_file_str +=".LD.gz";
        LDout = ifopen(LD_file_str, "w", InputFile::BGZF);

        beta_file_str += ".beta.gz";
        Beta_out = ifopen(beta_file_str, "w", InputFile::BGZF);

        if(LDout == NULL || Beta_out == NULL){
            perror("Fail to open LD or beta file!!! \n");
        }
    }else{
        LD_file_str +=".LD";
        LDout = ifopen(LD_file_str, "w", InputFile::UNCOMPRESSED);

        beta_file_str += ".beta";
        Beta_out = ifopen(beta_file_str, "w", InputFile::UNCOMPRESSED);

        if(LDout == NULL || Beta_out == NULL){
            perror("Fail to open LD or beta file!!! \n");
        }
    }

    ifprintf(Beta_out, "#ID\tCHR\tPOS\tREF\tALT\tbeta\tbeta_sd\tMAF\n");
    ifprintf(LDout, "#ID\tCHR\tPOS\tREF\tALT\tPOSvec\tLDvec\n");
    
    //Write file for LD matrix by the order of chr/bp
    for(size_t i=0; i<ni_test; i++){
        // write beta
        ifprintf(Beta_out, "%s\t%s\t%ld\t%s\t%s\t%g\t%g\t%g\n", snp_pos[i].rs.c_str(), snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), beta[i], beta_sd[i], snp_pos[i].maf);

        // write LD
        ifprintf(LDout, "%s\t%s\t%ld\t%s\t%s\t", snp_pos[i].rs.c_str(), snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str());

        for(size_t j=0; j<LD[i].size(); j++){
            ifprintf(LDout, "%ld,", snp_pos[i+j].bp);
        }
        ifprintf(LDout, "\t");

        for(size_t j=0; j<LD[i].size(); j++){
            ifprintf(LDout, "%g,", LD[i][j]);
        }
        ifprintf(LDout, "\n");

    }

    ifclose(LDout);
    ifclose(Beta_out);

    // tabix zipped files
    String cmd;
    int sys_status=1;

    if(zipSS){
        printf("\nTabixing .LD.gz files ... \n");
        cmd = String("tabix -c \"#\" -s 2 -b 3 -e 3 -f ") + LD_file_str;
        sys_status = system(cmd.c_str());
        if ( sys_status == 0 ) {
            printf( "\nLD output %s has been tabixed\n", LD_file_str.c_str() );
        }
        else {
            printf("\nUnable to tabix %s\n", LD_file_str.c_str());
        }

        printf("\nTabixing .beta.gz files ... \n");
        cmd = String("tabix -c \"#\" -s 2 -b 3 -e 3 -f ") + beta_file_str;
        sys_status = system(cmd.c_str());
        if ( sys_status == 0 ) {
            printf( "\nSS output %s has been tabixed\n", beta_file_str.c_str() );
        }
        else {
            printf("\nUnable to tabix %s\n", beta_file_str.c_str());
        }
    }

    return;
}


void getXy(const vector< vector<double> > &LD, const vector<double> &beta, vector <double> &Xty)
{
    Xty.clear();
    double xty_i;
    for(size_t i=0; i<beta.size(); i++){
        xty_i = beta[i] * LD[i][0];
        Xty.push_back(xty_i);
    }
    return;
}

void getPval(const vector<double> &beta, const vector<double> &beta_sd, vector <double> &pval, vector<pair<size_t, double> > &pos_ChisqTest)
{
    pval.clear();
    pos_ChisqTest.clear();
    double pval_i, chisq_i;

    for(size_t i=0; i<beta.size(); i++){
        chisq_i = pow(beta[i] / beta_sd[i], 2);
        pos_ChisqTest.push_back( make_pair(i, chisq_i) );

        pval_i = gsl_cdf_chisq_Q (chisq_i, 1);
        pval.push_back(pval_i);
    }
    return;
}

double getXtX(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j ){

    double xtx_ij = 0.0;

    if(pos_i == pos_j){
        xtx_ij = LD[pos_i][0];
    }
    else if(pos_i < pos_j)
    {
        if( (pos_j - pos_i) < LD[pos_i].size()  ) xtx_ij = LD[pos_i][pos_j - pos_i];        
    }else{
        if( (pos_i - pos_j) < LD[pos_j].size() ) xtx_ij = LD[pos_j][pos_i - pos_j];
    }

    return xtx_ij;
}


double getR2(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j ){

    double xtx_i, xtx_j, xtx_ij=0.0, r2_ij = 0.0;

    xtx_ij = getXtX(LD, pos_i, pos_j);

    if(xtx_ij > 0) { 
        xtx_i = LD[pos_i][0];
        xtx_j = LD[pos_j][0];
        r2_ij = pow(xtx_ij, 2) / ( xtx_i * xtx_j ) ; 
    }

    return r2_ij;
}





double CalcResVar(const gsl_matrix *XtX_cond, const gsl_vector * Xty_cond, const gsl_vector * beta_cond, const double &yty)
{
    double rtr, xtyb, bxtxb;
    gsl_vector *XtXb = gsl_vector_alloc(beta_cond->size);
    gsl_vector_set_zero (XtXb);

    gsl_blas_ddot(Xty_cond, beta_cond, &xtyb);

    gsl_blas_dgemv(CblasNoTrans, 1.0, XtX_cond, beta_cond, 0.0, XtXb);
    gsl_blas_ddot(XtXb, beta_cond, &bxtxb);

    rtr = yty - 2 * xtyb + bxtxb;

    gsl_vector_free(XtXb);

    if(rtr <= 0)
        perror("Nonpositive residual variance!\n");

    return rtr;
}














