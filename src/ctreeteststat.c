
#include "ctree.h"

/**
    Standardizes a statistic t of length pq with mean mu and covariance Sigma
    for variances > tol \n
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_standardize (const double *t, const double *mu, const double *Sigma, 
                    int pq, double tol, double *ans) {

    int i;
    double sd;
    
    for (i = 0; i < pq; i++) {
        sd = Sigma[i*pq + i]; 
        if (sd > tol)
            ans[i] = (t[i] - mu[i])/sqrt(sd);
        else
            ans[i] = 0.0;
    }
}


/**
    Absolute values of standardized statistics
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/

void C_absstandardize (const double *t, const double *mu, const double *Sigma,
                       int pq, double tol, double *ans) {
                      
    C_standardize(t, mu, Sigma, pq, tol, ans);
    C_abs_double(ans, pq);
}


/**
    Maximum absolute values of standardized statistics
    *\param t the vector of statistics
    *\param mu expectations
    *\param Sigma covariance matrix
    *\param pq dimension of t
    *\param tol tolerance for variances
*/

double C_maxabsTestStatistic (const double *t, const double *mu, const double *Sigma,
                              int pq, double tol) {
                           
    int i;
    double sd, maxabs = 0.0, tmp = 0.0;
    
    for (i = 0; i < pq; i++) {
        sd = Sigma[i*pq + i]; 
        if (sd > tol)
            tmp = fabs((t[i] - mu[i])/sqrt(sd));
        if (tmp > maxabs)
            maxabs = tmp;
    }
    return(maxabs);
}
