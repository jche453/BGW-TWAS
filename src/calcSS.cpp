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


//calculat summary statistics of score statistics and LD matrix
void CALCSS::GetSS(uchar **X, gsl_vector *y, vector< vector<double> > &LD, vector<double> &beta, vector<double> &beta_SE, vector<double> &U_STAT, vector<double> &SQRT_V_STAT, vector<double> &pval, vector<pair<size_t, double> > &pos_ChisqTest){

    cout << "\nStart calculating summary statistics ... \n";

    // Center y is centered by cPar.CopyPheno()
    gsl_blas_ddot(y, y, &pheno_var); 
    pheno_var /= ((double)(ni_test-1)) ;
    cout << "pheno_var = " << pheno_var << "\n";

    //cout << "create UcharTable ...\n";
    CreateUcharTable(UcharTable);

    // Create a vector of "SNPPOS" structs snp_pos for analyzed SNPs (ns_test)
    // (snpInfo will be cleared)
    CreateSnpPosVec(snp_pos, snpInfo, ns_total, indicator_snp); // The same order (pos) as in genotype file
    stable_sort(snp_pos.begin(), snp_pos.end(), comp_snp); // order snp_pos by chr/bp

    // define used variables 
    gsl_vector *xvec_i = gsl_vector_alloc(ni_test);
    gsl_vector *xvec_j = gsl_vector_alloc(ni_test);
    gsl_vector *xbeta_i = gsl_vector_alloc(ni_test);

    double xtx_ij, xty, xtx_i, beta_i, v2, chisq_i, beta_SE_i;
    beta.clear();
    xtx_vec.clear();
    LD.clear(); 
    pval.clear();
    pos_ChisqTest.clear();
    beta_SE.clear();

    cout << "calculate xtx, beta, score statistics by the order of chr/bp ... \n";
    for (size_t i=0; i<ns_test; ++i) {
        //calculate xtx_i
        getGTgslVec(X, xvec_i, snp_pos[i].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
        gsl_blas_ddot(xvec_i, xvec_i, &xtx_i);
        xtx_vec.push_back( xtx_i );

        //calculate effect-size
        gsl_blas_ddot(xvec_i, y, &xty);
        if(xtx_i > 0) beta_i = xty / xtx_i;
        else beta_i = 0.0;
        beta.push_back(beta_i); // effect size
        U_STAT.push_back(xty); // score statistic 
        v2 = pheno_var * xtx_i ;
        SQRT_V_STAT.push_back( sqrt(v2) ); // score statistic standard deviation
        chisq_i = xty * xty / v2; // chisq test statistic
        pval.push_back( gsl_cdf_chisq_Q (chisq_i, 1.0) ); // pvalue needed for BVSRM
        pos_ChisqTest.push_back( make_pair(i, chisq_i) ) ; // pos_ChisqTest needed for BVSRM

        gsl_vector_memcpy(xbeta_i, xvec_i);
        gsl_vector_scale(xbeta_i, -beta_i);
        gsl_vector_add(xbeta_i, y);
        gsl_blas_ddot(xbeta_i, xbeta_i, &beta_SE_i); // effect-size deviation
        if(xtx_i > 0) beta_SE_i = sqrt( beta_SE_i / ((double)ni_test * xtx_i) );
        else beta_SE_i = 0.0; 
        beta_SE.push_back(beta_SE_i);
        
        // saving X'X to LD
        LD.push_back(vector<double>()); // save r2
        LD[i].push_back(xtx_i / ((double)ni_test)); // varianct of x_i

        if(i < (ns_test-1) ){
            //calculate xtx_ij 
            for(size_t j=(i+1); j < ns_test; ++j){
                if( (snp_pos[j].chr == snp_pos[i].chr) && (snp_pos[j].bp <= snp_pos[i].bp + LDwindow) )
                {
                    getGTgslVec(X, xvec_j, snp_pos[j].pos, ni_test, ns_test, SNPmean, CompBuffSizeVec, UnCompBufferSize, Compress_Flag, UcharTable);
                    gsl_blas_ddot(xvec_i, xvec_j, &xtx_ij);
                    LD[i].push_back(xtx_ij / ((double)ni_test)); // covariance between x_i and x_j
                }
                else{break;}
            }
        }

    }

    gsl_vector_free(xvec_i);
    gsl_vector_free(xvec_j);
    gsl_vector_free(xbeta_i);

    return;
}


void CALCSS::WriteSS(const vector< vector<double> > &LD, const vector<double> &beta, const vector<double> &beta_SE, const vector<double> &U_STAT, const vector<double> &SQRT_V_STAT, const vector<double> &pval)
{
    cout << "\nStart writing summary statistics ... \n";
    String fout = file_out.c_str();

    // output files matches RareMetalWorker outputs
    String cov_file_str = "./output/" + fout;
    String score_file_str = "./output/" + fout;

    IFILE cov_out=NULL;
    IFILE score_out=NULL;

    if(zipSS){
        cov_file_str +=".cov.txt.gz";
        cov_out = ifopen(cov_file_str, "w", InputFile::BGZF);

        score_file_str += ".score.txt.gz";
        score_out = ifopen(score_file_str, "w", InputFile::BGZF);

        if(cov_out == NULL || score_out == NULL){
            perror("Fail to open LD or beta file!!! \n");
        }
    }else{
        cov_file_str +=".cov.txt";
        cov_out = ifopen(cov_file_str, "w", InputFile::UNCOMPRESSED);

        score_file_str += ".score.txt";
        score_out = ifopen(score_file_str, "w", InputFile::UNCOMPRESSED);

        if(cov_out == NULL || score_out == NULL){
            perror("Fail to open LD or beta file!!! \n");
        }
    }

    // write an extra column saving xtx with centered genotypes
    ifprintf(score_out, "#CHR\tPOS\tID\tREF\tALT\t N_INFORMATIVE\tFOUNDER_AF\tALL_AF\tINFORMATIVE_ALT_AC\tCALL_RATE\tHWE_PVALUE\tN_REF\tN_HET\tN_ALT\tU_STAT\tSQRT_V_STAT\tALT_EFFSIZE\tBETA_SE\tPVALUE\n");
    // assuming variants have unique CHR:POS 
    ifprintf(cov_out, "#CHR\tCURRENT_POS\tCURRENT_ID\tCURRENT_REF\tCURRENT_ALT\tMARKERS_IN_WINDOW\tCOV_MATRICES\n");
    
    double alt_ac;

    //Write files by the order of chr/bp
    for(size_t i=0; i<ns_test; i++){

        alt_ac = 2 * (double)ni_test * snp_pos[i].maf;
        // write score statistics
        ifprintf(score_out, "%s\t%ld\t%s\t%s\t%s\t%u\t%g\t%s\t%g\t%s\t%s\t%s\t%s\t%s\t%g\t%g\t%g\t%g\t%g\n", snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str(), ni_test, snp_pos[i].maf, "NA", alt_ac, "NA", "NA", "NA", "NA", "NA", U_STAT[i], SQRT_V_STAT[i], beta[i], beta_SE[i], pval[i]);

        // write banded covariance matrix: chr pos ref alt
        ifprintf(cov_out, "%s\t%ld\t%s\t%s\t%s\t", snp_pos[i].chr.c_str(), snp_pos[i].bp, snp_pos[i].rs.c_str(), snp_pos[i].a_major.c_str(), snp_pos[i].a_minor.c_str());

        for(size_t j=0; j<LD[i].size(); j++){
            ifprintf(cov_out, "%ld,", snp_pos[i+j].bp);
        }
        ifprintf(cov_out, "\t");

        for(size_t j=0; j<LD[i].size(); j++){
            ifprintf(cov_out, "%g,", LD[i][j]);
        }
        ifprintf(cov_out, "\n");

    }

    ifclose(cov_out);
    ifclose(score_out);

    // tabix zipped files
    String cmd;
    int sys_status=1;

    if(zipSS){
        printf("\nTabixing .cov.txt.gz files ... \n");
        cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + cov_file_str;
        sys_status = system(cmd.c_str());
        if ( sys_status == 0 ) {
            printf( "\nCOV output %s has been tabixed\n", cov_file_str.c_str() );
        }
        else {
            printf("\nUnable to tabix %s\n", cov_file_str.c_str());
        }

        printf("\nTabixing .score.txt.gz files ... \n");
        cmd = String("tabix -c \"#\" -s 1 -b 2 -e 2 -f ") + score_file_str;
        sys_status = system(cmd.c_str());
        if ( sys_status == 0 ) {
            printf( "\nScore statistics output %s has been tabixed\n", score_file_str.c_str() );
        }
        else {
            printf("\nUnable to tabix %s\n", score_file_str.c_str());
        }
    }

    return;
}

void Convert_LD(vector< vector<double> > &LD, vector<double> &xtx, const size_t &ns_test, const size_t &ni_test, const vector<SNPPOS> &snp_pos, const bool &refLD){
// convert cov matrix to LD r2 matrix, save n*diagonal to xtx vector
    xtx.clear();
    double r2, v2;
    vector<double> xtx_var;
    xtx_var.clear();

    cout << "Convert_LD ns_test = " << ns_test << "; ni_test = " << ni_test << endl;

    if(refLD){
        //use maf to calculate xtx
        cout << "Using MAF from score.txt \n";
        cout << "Print out LD from cov.txt:\n";
        for(size_t i=0; i<snp_pos.size(); ++i){
            //if(i < 10)
                //cout << snp_pos[i].key  << " maf = " << snp_pos[i].maf << "; " ;
            v2 = 2.0 * snp_pos[i].maf * (1.0 - snp_pos[i].maf);
            xtx_var.push_back(v2)  ;
            xtx.push_back( (double)ni_test * v2 ); // 

            for(size_t j = 0; j < LD[i].size(); j++){
                cout << LD[i][j] << ",";
            }
            cout << endl;
        }
    }
    else{
        //use xtx from cov matrix
        cout << "Using xtx/n from cov.txt \n";
        for(size_t i=0; i<ns_test; ++i){
            xtx_var.push_back(LD[i][0])  ;
            xtx.push_back( (double)ni_test * LD[i][0] ); // 
        }
    }
    cout << "set xtx vector success ! "<< endl;

    //cout << "LD size is " << LD.size() << endl;
    //cout << "LD[0] size is " << LD[0].size() << endl;

    for(size_t i=0; i<snp_pos.size(); i++){
        LD[i][0] = 1.0;
        //cout << xtx_var[i] << ",";

        if(i < (snp_pos.size() - 1)){
            if(xtx_var[i] == 0){
                for(size_t j=1; j< (LD[i].size()-1) ; j++)
                { LD[i][j] = 0.0 ; }
            }
            else{
                for(size_t j=1; j< (LD[i].size()-1) ; j++){
                    if(xtx_var[i+j] == 0.0)
                    { LD[i][j] = 0.0 ; }
                    else{
                        r2 = LD[i][j] / sqrt( xtx_var[i] * xtx_var[i+j] ) ;
                        //if(r2 < 1e-4) r2 = 0.0; // set r2 to 0 if r2 < 1e-4
                        LD[i][j] = r2 ;
                    }
                    //cout << LD[i][j] << "," ;             
                }
            }
        }
        //cout << endl;
    }

    return;
} 


void getXty(const vector<double> &beta, const vector<double> &xtx, vector <double> &Xty)
{
	//n is the sample size
    cout << "Calculate Xty ... \n";
    Xty.clear();
    for(size_t i=0; i<beta.size(); i++){
        Xty.push_back( beta[i] * xtx[i] );
    }
    return;
}

void getPval(const vector<double> &beta, const vector<double> &beta_sd, vector <double> &pval, vector<pair<size_t, double> > &pos_ChisqTest)
{
    cout << "Calculate pval ... \n";
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

double getXtX(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j, const vector<double> &xtx){

    double xtx_ij = 0.0;

    if(pos_i == pos_j){
        xtx_ij = xtx[pos_i];
    }
    else 
    {
        if( (pos_j - pos_i) > 0 && (pos_j - pos_i) < LD[pos_i].size()  ) 
            {
                xtx_ij = LD[pos_i][pos_j - pos_i] * sqrt(xtx[pos_i] * xtx[pos_j]);   
            }     
        else if( (pos_i - pos_j) > 0 && (pos_i - pos_j) < LD[pos_j].size() ) 
            {
                xtx_ij = LD[pos_j][pos_i - pos_j] * sqrt(xtx[pos_i] * xtx[pos_j]);
            }

    }

    return xtx_ij;
}






double getR2(const vector< vector<double> > &LD, const size_t &pos_i, const size_t &pos_j ){

    double r2_ij = 0.0;

    if(pos_i == pos_j) {
        r2_ij = 1.0;
    }
    else {
        if( (pos_j - pos_i) > 0 && (pos_j - pos_i) < LD[pos_i].size()  ) 
            {
                r2_ij = LD[pos_i][pos_j - pos_i] ;   
            }     
        else if( (pos_i - pos_j) > 0 && (pos_i - pos_j) < LD[pos_j].size() ) 
            {
                r2_ij = LD[pos_j][pos_i - pos_j]  ;
            }
    }

    r2_ij = r2_ij * r2_ij;

    return r2_ij;
}




double CalcResVar(const gsl_vector * Xty_cond, const gsl_vector * beta_cond, const double &yty)
{
    double rtr, xtyb;

    gsl_blas_ddot(Xty_cond, beta_cond, &xtyb);

    //cout << "Regression R2 in calcResVar = " << xtyb / yty << endl;

    rtr = yty - xtyb ;

    if(rtr <= 0){
        cout << "Regression R2 in calcResVar = " << xtyb / yty << endl;
        perror("Nonpositive residual variance!\n");  
    }

    return rtr;
}


void CalcBeta(const gsl_matrix *XtX_cond, const gsl_vector * Xty_cond, gsl_vector * beta_cond)
{
    size_t s_size = Xty_cond->size;

    gsl_matrix *XtXinv = gsl_matrix_alloc(s_size, s_size);

    gsl_matrix_memcpy(XtXinv, XtX_cond);

    LapackSolve(XtXinv, Xty_cond, beta_cond);

    gsl_matrix_free(XtXinv);

    return ;
}












