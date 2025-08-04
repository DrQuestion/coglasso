#include "math.h"
#include <vector>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

void coglasso_sub(Eigen::MatrixXd& S, Eigen::MatrixXd& W, Eigen::MatrixXd& T, int pX, int p, double ialpha, double ic, Eigen::MatrixXd& Lambda_star, int& df, int& converged, bool scr, unsigned int index);

//[[Rcpp::export]]
List co_glasso(Eigen::Map<Eigen::MatrixXd> S, int pX, Eigen::Map<Eigen::MatrixXd> hpars, Eigen::Map<Eigen::MatrixXd> T_guess, bool scr, bool verbose, bool cov_output)
{
  // Caller function, initiates what must be initiated and calls the main algorithm for 
  // every combination of hyperparameters.
    unsigned int nhpars = hpars.rows();
    int nfeatures = S.rows();
    auto S_diag = S.diagonal().array();
    int& p = nfeatures;
    
    // Anti explosion warm start initialization strategy
    int size_c = 0;
    int nexploded = 0;
    int last_c;
    IntegerVector idx_c(nhpars);

    bool zero_sol = true;

    // define result list: loglik
    NumericVector loglik(nhpars, -nfeatures), density(nhpars);
    IntegerVector df(nhpars), convergence(nhpars, 0);
    List result;
    result["loglik"] = loglik;
    result["density"] = density;
    result["df"] = df;
    // temporary returned info to debug explosions
    result["convergence"] = convergence;

    std::vector<Eigen::MatrixXd*> tmp_icov_p, tmp_cov_p, tmp_path_p;

    // HYPERPARAMETERS FOR LOOP - iterates through the combinations of hyperparameters from 
    // the least to the most stringent, and for any of those executes the coglasso algorithm.
    for (int i = nhpars - 1; i >= 0; i--) {

        // Shut down warm start - Warm start initialization is suppressed for the moment with
        // the following line
        zero_sol = true;

        // Warm start initialization
        int converged = 0;

        // TODO: try to generalize out of pX and PM
        // TODO: try reasoning with "sorted" matrix ncombinations x 3 for (alpha,lambda_w,lambda_b) 
        // rather than iterating through parameters --> done, not sure of sorting strategy
        
        double alpha = hpars(i, 0);
        double lambda_w = hpars(i, 1);
        double lambda_b = hpars(i, 2);
        double c = hpars(i, 3);

        // Matrix of lambda_w and lambda_b multiplied by alpha is built here
        MatrixXd Lambda_star_ij(p, p);
        for (int row_i = 0; row_i < p; row_i++) {
            for (int col_i = 0; col_i < p; col_i++) {
                if ((row_i < pX && col_i < pX) || (row_i >= pX && col_i >= pX)) {
                    Lambda_star_ij(row_i, col_i) = alpha * lambda_w;
                    //Lambda_ij(row_i, col_i) = lambda_w;
                }
                else {
                    Lambda_star_ij(row_i, col_i) = alpha * lambda_b;
                    //Lambda_ij(row_i, col_i) = lambda_b;
                }
            }
        }
        
        // pre-screening by z
        // z vector of row (hence col) indexes for which more than one element 
        // is bigger than the ith lambda
        // q is the number of rows for which more than one element is bigger
        // than the ith lambda
        vector<int> z;
        for (int row_i = 0; row_i < p; row_i++) {
            int break_flag = 0;
            for (int col_i = 0; col_i < p; col_i++) {
                if (break_flag > 1) break;
                // should I use alpha*S (or maybe just lambda) since I'll be subtracting from alpha*S?
                if (alpha * S(row_i, col_i) > Lambda_star_ij(row_i, col_i) or alpha * S(row_i, col_i) < -Lambda_star_ij(row_i, col_i)) break_flag++;
                // if (S(row_i, col_i) > Lambda_star_ij(row_i, col_i) or S(row_i, col_i) < -Lambda_star_ij(row_i, col_i)) break_flag++;
            }
            if (break_flag > 1) z.push_back(row_i);
        }
        int q = z.size();
        MatrixXd sub_S(p, p), sub_W(p, p), sub_T(p, p);
        int sub_df = 0;

        // Rcout << "\n" << q << "\n";

        // code for z screening is commented (and development paused) until it is not clear how (and why) it works
        // #pragma omp parallel for
        // create the submatrix of S that corresponds to the indexes for which
        // more than one element is bigger than the ith lambda, and W and T of same dimensions
        // if zero_sol, sub_W is equal to sub_S and sub_T identity, otherwhise the covariance
        // and the inverse covariance found with the previous ilambda as a "warm" starting 
        // point for sub_W and sub_T
        // int sub_pT; would have been useful if screening with z
        // for (int ii = 0; ii < q; ii++) {
        //     for (int jj = 0; jj < q; jj++) {
        //         sub_S(ii, jj) = S(z[ii], z[jj]);
        //         // find new sub_pT if (z[ii])
        //         if (zero_sol) {
        //             sub_W(ii, jj) = S(z[ii], z[jj]);
        //             sub_T(ii, jj) = ii == jj ? 1 : 0;
        //         }
        //         else {
        //             sub_W(ii, jj) = (*(tmp_cov_p.back()))(z[ii], z[jj]);
        //             sub_T(ii, jj) = (*(tmp_icov_p.back()))(z[ii], z[jj]);
        //         }
        //     }
        // }

        // Without anti explosion
        for (int ii = 0; ii < p; ii++) {
            for (int jj = 0; jj < p; jj++) {
                sub_S(ii, jj) = S(ii, jj);
                if (zero_sol) {
                    // maybe initialize W to alpha * S?
                    sub_W(ii, jj) = S(ii, jj);
                    sub_T(ii, jj) = ii == jj ? 1 : T_guess(ii, jj);
                  }
                else {
                    sub_W(ii, jj) = (*(tmp_cov_p.back()))(ii, jj);
                    sub_T(ii, jj) = (*(tmp_icov_p.back()))(ii, jj);
                }
            }
        }

        // Anti explosion warm start
        // for (int ii = 0; ii < p; ii++) {
        //     for (int jj = 0; jj < p; jj++) {
        //         sub_S(ii, jj) = S(ii, jj);
        //         if (zero_sol) {
        //             if (size_c > 0) {
        //                 last_c = idx_c[size_c - 1];
        //                 sub_W(ii, jj) = (*(tmp_cov_p[last_c]))(ii, jj);
        //                 sub_T(ii, jj) = (*(tmp_icov_p[last_c]))(ii, jj);
        //             }
        //             else {
        //                 sub_W(ii, jj) = S(ii, jj);
        //                 sub_T(ii, jj) = ii == jj ? 1 : 0;
        //             }
        //         }
        //         else {
        //             if (size_c > 0) {
        //                 last_c = idx_c[size_c - 1];
        //                 sub_W(ii, jj) = (*(tmp_cov_p[last_c]))(ii, jj);
        //                 sub_T(ii, jj) = (*(tmp_icov_p[last_c]))(ii, jj);
        //             }
        //             else {
        //                 sub_W(ii, jj) = S(ii, jj);
        //                 sub_T(ii, jj) = ii == jj ? 1 : 0;
        //             }
        //             
        //         }
        //     }
        // }

        if (q > 0)
        {
            if (verbose) {
                if (scr)
                    Rcout << "Conducting the collaborative graphical lasso (coglasso) wtih lossy screening....in progress: " << floor(100 * (1. - (1. * i / nhpars))) << "%\r";
                if (!scr)
                    Rcout << "Conducting the collaborative graphical lasso (coglasso)....in progress: " << floor(100 * (1. - (1. * i / nhpars))) << "%\r";
            }

            coglasso_sub(sub_S, sub_W, sub_T, pX, p, alpha, c, Lambda_star_ij, sub_df, converged, scr, i);
            zero_sol = false;
        }

        if (q == 0) { 
            zero_sol = true; 
            nexploded -= 1;
            // Rcout << "q was zero for index=" << i << "\n\n";
        }

        // Rcout << converged << "\n";

        if (converged == 1) {
            idx_c[size_c] = nhpars - 1 - i;
            size_c += 1;
            convergence[i] = 1;
        }
        else {
            // Rcout << "\n" << alpha << "\t" << lambda_w << "\t" << lambda_b << "\n\n";
            nexploded += 1;
        }

        

        // update result list
        tmp_path_p.push_back(new Eigen::MatrixXd(p, p));
        tmp_icov_p.push_back(new Eigen::MatrixXd(p, p));
        tmp_cov_p.push_back(new Eigen::MatrixXd(p, p));

        Eigen::MatrixXd* tmp_icov, * tmp_cov, * tmp_path;
        tmp_icov = tmp_icov_p.back();
        tmp_cov = tmp_cov_p.back();
        tmp_path = tmp_path_p.back();
        tmp_icov->setZero();
        tmp_icov->diagonal() = 1 / (S_diag + double(lambda_w));
        tmp_cov->setZero();
        tmp_cov->diagonal() = S_diag + double(lambda_w);
        tmp_path->setZero();

        if (!zero_sol)
        {
            //#pragma omp parallel for
            for (int ii = 0; ii < p; ii++) {
                for (int jj = 0; jj < p; jj++) {
                    (*tmp_icov)(ii, jj) = sub_T(ii, jj);
                    (*tmp_cov)(ii, jj) = sub_W(ii, jj);
                    (*tmp_path)(ii, jj) = sub_T(ii, jj) == 0 ? 0 : 1;
                    (*tmp_path)(ii, jj) = ii == jj ? 0 : (*tmp_path)(ii, jj);
                }
            }
            
            // No screening over z is done, next chunk commented
            // for (int ii = 0; ii < q; ii++) {
            //     for (int jj = 0; jj < q; jj++) {
            //         (*tmp_icov)(z[ii], z[jj]) = sub_T(ii, jj);
            //         (*tmp_cov)(z[ii], z[jj]) = sub_W(ii, jj);
            //         (*tmp_path)(z[ii], z[jj]) = sub_T(ii, jj) == 0 ? 0 : 1;
            //         (*tmp_path)(z[ii], z[jj]) = ii == jj ? 0 : (*tmp_path)(z[ii], z[jj]);
            //     }
            // }
            density[i] = 1.0 * sub_df / p / (p - 1);
            df[i] = sub_df / 2;
            loglik[i] = log(sub_T.determinant()) - (sub_T * sub_S).diagonal().sum();
        }

    }
    
    if (verbose) {
      if (scr)
        Rcout << "Conducting the collaborative graphical lasso (coglasso) wtih lossy screening....done            \n";
      if (!scr)
        Rcout << "Conducting the collaborative graphical lasso (coglasso)....done            \n";
    }
    
    List path, icov, cov;
    for (unsigned int i = 0; i < nhpars; i++) {
        path.push_back(*(tmp_path_p[nhpars - 1 - i]));
        icov.push_back(*(tmp_icov_p[nhpars - 1 - i]));
        if (cov_output) cov.push_back(*(tmp_cov_p[nhpars - 1 - i]));
    }

    // Rcout << nexploded << "\t";

    result["path"] = path;
    result["icov"] = icov;
    if (cov_output) result["cov"] = cov;
    
    // temporary returned info to debug explosions
    result["nexploded"] = nexploded;
    return result;
}

// Main function, run for each combination of hyperparameters. This algorithm uses the coglasso update 
// rule to optimize the matrixes T and W based on the current combination of hyperparameters lambda_w,
// lambda_b and alpha. For more info on what ere these obects, please read the inputs to the algorithm,
// that are described in the following lines.
// 
// S is the covariance matrix computed from data;
// W is the estimated covariance matrix. When warm start initialization is possible, the last available W 
//      matrix will be used, otherwise W is initialized equal to S (then inside the main algorithm we would 
//      add lambda_w to the diagonal);
// T is the matrix of coefficients (betas). Its diagonal-normalized form corresponds to the partial 
//      correlation matrix. When warm start initialization is possible, the last available T matrix will be
//      used, otherwise T is initialized as the identity matrix;
// pX is the number of variables in the first data set. In practice it corresponds to the dimension 
//      of the block in the matrix T of the interactions from the first data set to itself;
// p is the dimension of our matrixes, the total number of variables;
// ialpha is the alpha of the current combination of iperparameters. Remember that alpha is 1/(1+c) where c
//      is the weight of the collaborative term of coglasso's expression;
// Lambda_star is the matrix of lambda_w and lambda_b multiplied by alpha, where lambda_w is the penalization
//      for the "within" interactions and lambda_b is the penalization for the "between" interactions;
// df is the number of degrees of freedom (?) of matrix T (the number of non-zero off-diagonal elements). 
//      When warm start initialization is possible, the df number of the last available T matrix will be
//      used, otherwise df is initialized to 0 (as T is the identity matrix);
// scr is a parameter coming from the original glasso implementation, which I am thinking to remove. Not 
//      used in the current coglasso implementation;
// index is the index of the current combination of hyperparameters (that goes from number of combinations to 0).
//
void coglasso_sub(Eigen::MatrixXd& S, Eigen::MatrixXd& W, Eigen::MatrixXd& T, int pX, int p, double ialpha, double ic, Eigen::MatrixXd& Lambda_star, int& df, int& converged, bool scr, unsigned int index)
{
    int i, j, k; //initialize indices
    int rss_idx, w_idx;

    int gap_int;
    double gap_ext, gap_act;
    double thol_act = 1e-4;
    double thol_ext = 1e-4;

    int MAX_ITER_EXT = 100;
    int MAX_ITER_INT = 10000;
    int MAX_ITER_ACT = 10000;
    int iter_ext, iter_int, iter_act;



    Eigen::MatrixXi idx_a(p, p); // active set
    Eigen::MatrixXi idx_i(p, p); // The set possibly can join active set
    int* size_a = (int*)malloc(p * sizeof(int)); //sizes of active sets
    double* w1 = (double*)malloc(p * sizeof(double));
    double* ww = (double*)malloc(p * sizeof(double));


    int size_a_prev; //original size of the active set
    int junk_a; //the number of coefficients returning to the inactive set from the active set

    double r; //partial residual
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;

    // The following for loop is just to initialize the sets of active and inactive
    // coefficients from the matrix T given as input to the core function. T can be
    // or the identity matrix or the warm start initialized matrix (the last T matrix
    // that converged). If coefficients are different from 0 they are classified as
    // active, otherwise as inactive.
    //
    for (i = 0; i < p; i++) {

        W(i, i) = S(i, i) + Lambda_star(i,i)/ialpha; 
        // The diagonal elements are set optimal ???? 
        // TODO: check final likelihood equation form to determine wich is the optimal
        // diagonal element (checked for the case of mixed penalties (not yet for 
        // coglasso per se), should be ok like this)
        size_a[i] = 0;
        tmp1 = T(i, i);
        T(i, i) = 0;

        for (j = 0; j < p; j++) {
            if (scr) // NOT USED REALLY
                if (fabs(S(j, i)) <= Lambda_star(j,i)) {
                    idx_i(j, i) = -1;
                    T(j, i) = 0;
                    continue;
                }

            if (T(j, i) != 0) {
                idx_a(size_a[i], i) = j; //initialize the active set
                size_a[i]++;
                idx_i(j, i) = -1; //initialize the inactive set
                T(j, i) = -T(j, i) / tmp1;
            }
            else idx_i(j, i) = 1;
        }
        idx_i(i, i) = -1;
    }
    
    // EXTERNAL WHILE LOOP - updates the matrix T by running the INTERNAL WHILE LOOP
    // column by column:
    // This loop iterates trhough the columns of the matrix of coefficients T with the    
    // for loop that begins few lines below this point. For each column, the INTERNAL 
    // WHILE LOOP is applied until that column does not converge. After all columns have
    // converged, the current iteration of the EXTERNAL WHILE LOOP is over, and the 
    // convergence criterion is computed. It consists of the normalized absolute
    // difference between the matrix T before the current iteration of the external
    // loop and after it is over (and hence after T's update). If this value gap_ext is
    // < 0.0001, the EXTERNAL WHILE LOOP ends. For additional explanations on the 
    // convergence criterion of this loop, please look after the end of the INTERNAL 
    // WHILE LOOP.
    //
    gap_ext = 1;
    iter_ext = 0;
    while (gap_ext > thol_ext && iter_ext < MAX_ITER_EXT) //outer loop
    {
        tmp1 = 0; // What for??? Not used in the loop!
        tmp6 = 0;
        tmp5 = 0;
        for (i = 0; i < p; i++)
        {


            gap_int = 1;
            iter_int = 0;

            for (j = 0; j < p; j++)
                ww[j] = T(j, i);

            // INTERNAL WHILE LOOP - activates inactive coefficients when appropriate, 
            // updates active coefficients and stops when the set of active coefficients 
            // does not change after updating them:
            // This internal while loop updates the ith column of the matrix T it 
            // received from the EXTERNAL WHILE LOOP until convergence. First, it 
            // computes the size of the set of active coefficients (size_a_prev) in the
            // column. Active coefficents are those that after the last update have not
            // been set to 0. Then it decides which coefficients in the column of matrix
            // T must be activated. It does so with the ACTIVATING FOR LOOP. It computes
            // the difference in size of the set of active coefficients before and after
            // the ACTIVATING FOR LOOP (gap_int = size_a_prev - size_a [computed after the 
            // ACTIVATING FOR LOOP]). It then updates the active coefficients until 
            // convergence with the ACTIVE COEFFICIENTS WHILE LOOP. 
            // Convergence of the INTERNAL WHILE LOOP happens when the difference gap_int
            // computed after the ACTIVATING FOR LOOP was 0 (did not manage to activate 
            // any new coefficient), but not before a final round of the ACTIVE COEFFICIENTS WHILE LOOP.
            //
            while (gap_int != 0 && iter_int < MAX_ITER_INT)
            {
                size_a_prev = size_a[i];
                // TODO: generalize following loop to more than two dataset

                // ACTIVATING FOR LOOP - activates inactive coefficients when appropriate:
                // This loop iterates through inactive (idx_i(j, i) != -1) coefficients
                // of the ith column of matrix T and activates them when r is outside 
                // the interval [-ialpha*lambda, ialpha*lambda].
                // In practice it runs for a single iteration the update rule for each 
                // coefficient and it updates it (and activates it) when it does not remain 0.
                // All active coefficients will enter the following ACTIVE COEFFICIENTS WHILE LOOP.
                //
                for (j = 0; j < p; j++)
                {
                    if (idx_i(j, i) != -1)
                    {
                        r = ialpha*S(j, i); 
                        // in theory to exclude T(w_idx, i) * W(w_idx, w_idx) from 
                        // the following summation of subtractions (k != j in the 
                        // summation of subtractions of the formula) we should add
                        // it here as in the following loop, but we know that T(w_idx, i)
                        // if the coefficient is inactive

                        if (j < pX) {
                            for (k = 0; k < size_a[i]; k++)
                            {
                                rss_idx = idx_a(k, i);
                                if (rss_idx < pX) {
                                    r = r - W(j, rss_idx) * T(rss_idx, i);
                                }
                                else {
                                    // r = r + (1 - ialpha) * W(j, rss_idx) * T(rss_idx, i);
                                    r = r - ialpha * (1 - ic) * W(j, rss_idx) * T(rss_idx, i);
                                }
                                
                            }
                        }
                        else {
                            for (k = 0; k < size_a[i]; k++)
                            {
                                rss_idx = idx_a(k, i);
                                if (rss_idx < pX) {
                                    r = r - ialpha * (1 - ic) * W(j, rss_idx) * T(rss_idx, i);
                                }
                                else {
                                    r = r - W(j, rss_idx) * T(rss_idx, i);
                                }
                            }
                        }
                        
                        if (r > Lambda_star(j, i))
                        {
                            w1[j] = (r - Lambda_star(j, i)) / W(j, j);
                            idx_a(size_a[i], i) = j;
                            size_a[i] = size_a[i] + 1;
                            idx_i(j, i) = -1;
                        }

                        else if (r < -Lambda_star(j, i))
                        {
                            w1[j] = (r + Lambda_star(j, i)) / W(j, j);
                            idx_a(size_a[i], i) = j;
                            size_a[i] = size_a[i] + 1;
                            idx_i(j, i) = -1;
                        }

                        else {
                            w1[j] = 0;
                        }

                        T(j, i) = w1[j];
                    }
                }
                
                // if (iter_ext == 1 && i == 1 && iter_int==0) {
                //     Rcout << iter_int << "\n" << S << "\n\n" << W << "\n\n" << T << "\n\n";
                // }
                
                gap_int = size_a[i] - size_a_prev;

                gap_act = 1;
                iter_act = 0;

                // ACTIVE COEFFICIENTS WHILE LOOP - updates active coefficients until convergence:
                // In all iterations of this active coefficients while loop we iterate 
                // through the active coefficients of the ith column with the for loop 
                // indexing active coefficients with j. For every jth active coefficient 
                //   i) we apply the update rule once
                //  ii) we sum the absolute value of the newly obtained coefficient to tmp4
                // iii) we compute the absolute difference between the jth coefficient at 
                //      the previous iteration and we sum it to tmp3 and
                //  iv) we immediately update the jth entry of the ith column of matrix T.
                // Through the iterations of this for loop, tmp3 accumulates the absolute
                // differences of all the updated active coefficients and tmp4 the sum of 
                // the absolute values of the updated coefficients.
                // When this for loop through active coefficients ends, the convergence 
                // criterion tmp3 / tmp4 is computed. If it is less than 0.0001 we say 
                // that the ith column has converged for the current active coefficients while 
                // loop and we exit. In practice we exit when the ratio between the 
                // absolute distances current - previous and the absolute sum of the current
                // estimated column of T is less than 0.0001. Otherwise we update once again.
                // This is the loop where coefficients go to infinity at a certain point.
                // 
                while (gap_act > thol_act && iter_act < MAX_ITER_ACT)
                {
                    tmp3 = 0;
                    tmp4 = 0;
                 
                    for (j = 0; j < size_a[i]; j++)
                    {
                        w_idx = idx_a(j, i);

                        // here we iterate only through the active coefficients (the old 
                        // and the new ones coming from previous loop)
                        if (w_idx != -1)
                        {
                            //tmp_a = w_idx*p;
                            r = ialpha*S(w_idx, i) + T(w_idx, i) * W(w_idx, w_idx);
                            // we exclude T(w_idx, i) * W(w_idx, w_idx) from the following
                            // summation of subtractions (k != j in the summation of 
                            // subtractions of the formula) by adding it here

                            if (w_idx < pX) {
                                for (k = 0; k < size_a[i]; k++)
                                {
                                    rss_idx = idx_a(k, i);
                                    if (rss_idx < pX) {
                                        r = r - W(w_idx, rss_idx) * T(rss_idx, i);
                                    }
                                    else {
                                        // r = r + (1 - ialpha) * W(w_idx, rss_idx) * T(rss_idx, i);
                                        r = r - ialpha * (1 - ic) * W(w_idx, rss_idx) * T(rss_idx, i);
                                    }
                                }
                            }

                            else {
                                for (k = 0; k < size_a[i]; k++)
                                {
                                    rss_idx = idx_a(k, i);
                                    if (rss_idx < pX) {
                                        // r = r + (1 - ialpha) * W(w_idx, rss_idx) * T(rss_idx, i);
                                        r = r - ialpha * (1 - ic) * W(w_idx, rss_idx) * T(rss_idx, i);
                                    }
                                    else {
                                        r = r - W(w_idx, rss_idx) * T(rss_idx, i);
                                    }
                                }
                            }

                            if (r > Lambda_star(w_idx, i)) {
                                w1[w_idx] = (r - Lambda_star(w_idx, i)) / W(w_idx, w_idx);
                                tmp4 += w1[w_idx];
                            }

                            else if (r < -Lambda_star(w_idx, i)) {
                                w1[w_idx] = (r + Lambda_star(w_idx, i)) / W(w_idx, w_idx);
                                tmp4 -= w1[w_idx];
                            }

                            else w1[w_idx] = 0;  

                            // if (index == 6 && iter_ext == 1 && i == 7 && iter_int == 0 && iter_act < 20) {
                            //     Rcout << w1[w_idx] << "\t" << fabs(w1[w_idx] - T(w_idx, i)) << "\n";
                            // }

                            // if (index == 6 && iter_ext == 1 && i == 7 && iter_int == 0 && iter_act < 20) {
                            //     Rcout << w1[w_idx] << "\t" << fabs(w1[w_idx] - T(w_idx, i)) << "\n";
                            // }

                            tmp3 = tmp3 + fabs(w1[w_idx] - T(w_idx, i));

                            T(w_idx, i) = w1[w_idx];

                        }
                    }
                    gap_act = tmp3 / tmp4;
                    // Easy to see how iter by iter the values explode (both negatively 
                    // and positively) while THEY SHOULD CONVERGE instead. Anything wrong 
                    // with the rule? How to converge instead?
                    //if (iter_ext >= 2 && iter_act < 3) {
                    //    Rcout << iter_ext << "\t" << i << "\t" << iter_int << "\t" << iter_act << "\n" << T << "\n\n";
                    //}

                    // if (iter_ext == 1 && i == 1 && iter_int == 0 && iter_act<20) {
                    //     Rcout << iter_act << "\n" << T << "\n\n";
                    // }

                    iter_act++;
                }
                
                //Here we classify as inactive the coefficients that have been updated to 0
                junk_a = 0;
                for (j = 0; j < size_a[i]; j++) {
                    w_idx = idx_a(j, i);
                    if (w1[w_idx] == 0) {
                        junk_a++;
                        idx_i(w_idx, i) = 1;
                        idx_a(j, i) = -1;
                    }
                    else idx_a(j - junk_a, i) = w_idx;
                }
                size_a[i] = size_a[i] - junk_a;

                // Rcout << iter_ext << "\t" << i << "\t" << iter_int << "\n\n";

                iter_int++;
            }

            // The convergence of the INTERNAL WHILE LOOP could be seen also as the
            // convergence of the ith column of the matrix T during the current iteration
            // of the EXTERNAL WHILE LOOP. After the INTERNAL WHILE LOOP converges, (hence
            // the ith column has been updated), we compute the absolute difference between
            // the updated and the previous elements of the ith column and we sum it to 
            // tmp5. We also take tmp4 (containing the sum of the absolute values of the
            // updated ith column of T) and we sum it to tmp6. The ith row and column of 
            // the matrix W are updated. The updated ith row (and column) are obtained done
            // by multiplying (by matrix multiplication) the (just obtained) updated ith 
            // column of T by the W obtained after the previous iteration of the for loop
            // through columns (when the i-1 th row and column have updated). At the end of
            // the for loop through all ith columns (when the INTERNAL WHILE LOOP will hence
            // be applied to all the columns of T), tmp5 will store the sum of absoulute 
            // distances between the updated matrix T and the previous one. Instead, tmp6
            // will contain the sum of the absolute values of the current T matrix after
            // convergence of the internal loop for all the columns.
            // So when this for loop through columns is exited, the criterion gap_ext will 
            // contain the ratio between the absolute difference between the T computed in
            // the previous iteration of the external loop and the current one over the
            // absolute sum of the values of the current T. If this ratio is below 0.0001 the
            // external loop converges.
            //
            //update W Beta
            Eigen::MatrixXd temp = W.transpose() * T.col(i);
            for (j = 0; j < i; j++) {
                W(j, i) = temp(j, 0);
                W(i, j) = temp(j, 0);
            }
            for (j = i + 1; j < p; j++) {
                W(j, i) = temp(j, 0);
                W(i, j) = temp(j, 0);
            }

            for (j = 0; j < p; j++)
                tmp5 = tmp5 + fabs(ww[j] - T(j, i));
            tmp6 = tmp6 + tmp4;
        }
        gap_ext = tmp5 / tmp6;
        // Rcout << gap_ext << "\n";
        //printf("%g\n",gap_ext);
        iter_ext++;

        // Rcout << iter_ext << "\n" << T << "\n\n";
    }

    for (i = 0; i < p; i++) //Compute the final T
    {
        tmp2 = 0;
        tmp2 = W.col(i).transpose() * T.col(i) - W(i, i) * T(i, i);

        tmp1 = 1 / (W(i, i) - tmp2);
        T.col(i) *= -tmp1;
        T(i, i) = tmp1;
    }
    for (i = 0; i < p; i++)
        df += size_a[i];

    // if (isnan(T(0, 0))) Rcout << Lambda_star(0, 0) / ialpha << "\t" << Lambda_star(0, p - 1) / ialpha << "\t" << Lambda_star(0, 0) << "\t" << Lambda_star(0, p - 1) << "\t" << ialpha << "\t" << "Did not converge\n";
    // else Rcout << Lambda_star(0, 0) / ialpha << "\t" << Lambda_star(0, p - 1) / ialpha << "\t" << Lambda_star(0, 0) << "\t" << Lambda_star(0, p - 1) << "\t" << ialpha  << "\n";
    // if (isnan(T(0, 0))) Rcout << Lambda_star(0, 0) / ialpha << "\t" << Lambda_star(0, p - 1) / ialpha << "\t" << Lambda_star(0, 0) << "\t" << Lambda_star(0, p - 1) << "\t" << ialpha << "\n";
    if (!isnan(T(0, 0))) converged = 1;
    // Rcout << converged << "\n";

    // Rcout << T << "\n\n";

    free(size_a);
    free(w1);
    free(ww);
}