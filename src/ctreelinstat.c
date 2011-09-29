
#include "ctree.h"

/**
    Linear statistics for conditional inference based on Strasser & Weber (1999)
    *\file Linstat.c
    *\author $Author: hothorn $
    *\date $Date: 2006-08-25 10:53:10 +0200 (Fri, 25 Aug 2006) $
*/

/**
    Computes the linear statistic, formula (1) in the paper\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; a pointer to a REALSXP-vector of length pq
*/
  
void C_Linstat_numeric (const double *x, const double *y, 
                        const int q, const int *weights, const int n,
                        double *ans) {
              
    int i, k;
    double tmp;

    for (k = 0; k < q; k++)
        ans[k] = 0.0;
            
    for (i = 0; i < n; i++) {
                
        /* optimization: weights are often zero */
        if (weights[i] == 0) continue;
 
        tmp = x[i] * weights[i];
        for (k = 0; k < q; k++)  
            ans[k] += tmp * y[k * n + i];
    }
}


void C_Linstat_factor (const int *x, const int p,
                       const double *y, const int q,
                       const int *weights, const int n,
                       double *ans) {
              
    int i, j, k, kp, kn, pq;
    double tmp;


    pq = p * q;
    for (k = 0; k < pq; k++) ans[k] = 0.0;

    for (k = 0; k < q; k++) {

        kn = k * n;
        kp = k * p;
            
        for (i = 0; i < n; i++) {
                
            /* optimization: weights are often zero, treatment contrasts */
            if (weights[i] == 0) continue;
                
            /* FIXME: treatment contrasts? */
            tmp = y[kn + i] * weights[i];
            ans[kp + (x[i] - 1)] += tmp;
        }
    }
}

/**
    Conditional expectation and covariance of the influence function\n
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; an object of class `ExpectCovarInfluence'
*/

void C_ExpInf (const double* y, const int q, const int* weights, 
               const int sw, const int n, double* ans) {

    int i, j;    
    
    for (j = 0; j < q; j++) ans[j] = 0.0;
    
    if (sw <= 1) 
        error("C_ExpInf: sum of weights is less than one");

    /*
     *    Expectation of the influence function
     */

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
    
        if (weights[i] == 0) continue;
    
        for (j = 0; j < q; j++)
            ans[j] += weights[i] * y[j * n + i];
    }

    for (j = 0; j < q; j++)
        ans[j] = ans[j] / sw;

}

/**
    Conditional expectation and covariance of the influence function\n
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param ans return value; an object of class `ExpectCovarInfluence'
*/

void C_CovInf (const double* y, const int q, const int* weights, const int sw, 
               const int n, const double* exp, double *ans) {

    int i, j, k, jq;
    double tmp;
    
    for (j = 0; j < q*q; j++) ans[j] = 0.0;

    if (sw <= 1) 
        error("C_CovInf: sum of weights is less than one");
    
    /*
     *    Covariance of the influence function
     */

    for (i = 0; i < n; i++) {

        if (weights[i] == 0) continue;
     
        for (j = 0; j < q; j++) {
            tmp = weights[i] * (y[j * n + i] - exp[j]);
            jq = j * q;
            for (k = 0; k < q; k++)
                ans[jq + k] += tmp * (y[k * n + i] - exp[k]);
        }
    }

    for (j = 0; j < q*q; j++)
        ans[j] = ans[j] / sw;
}

void C_swx_numeric (const double *x, const int *weights, const int n, 
                    double *swx, double *swx2) {

    int i;
    double tmp;

    swx[0] = 0.0;
    swx2[0] = 0.0;

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
        if (weights[i] == 0) continue;
    
        tmp = weights[i] * x[i];
        swx[0] += tmp;

        /* covariance part */
        swx2[0] += tmp * x[i];
    }
}

void C_swx_factor (const int *x, const int p, const int *weights, const int n,
                   double *swx, double *swx2) {

    int i, j, k;

    for (j = 0; j < p; j++) {
        swx[j] = 0.0;
        for (k = 0; k < p; k++)
            swx2[j * p + k] = 0.0;
    }

    for (i = 0; i < n; i++) {

        /*  observations with zero case weights do not contribute */
        if (weights[i] == 0) continue;

        /* FIXME: treatment contrasts ? */
        swx[x[i] - 1] += (double) weights[i];
    }
    
    /* covariance part */
    for (j = 0; j < p; j++)
        swx2[j * p + j] = swx[j];
}


/**
    Conditional expectation and covariance of the a linear statistic\n
    *\param x values of the transformation
    *\param p dimension of the transformation
    *\param y values of the influence function
    *\param q dimension of the influence function
    *\param weights case weights
    *\param n number of observations
    *\param expcovinf an object of class `ExpectCovarInfluence'
    *\param ans return value; an object of class `ExpectCovar'
*/

void C_ExpCovLinstat (const double *swx, const double *swx2, const int p, const int q, 
                      const int sw, const double *expinf, const double *covinf, 
                      double *explinstat, double *covlinstat) {

    int i, j, k, pq, ip;
    double f1, f2, tmp, dsw;
    double *CT2, *Covy_x_swx;

    pq = p * q;
    
    /*
    *   explinstat: expectation of the linear statistic T
    */

    for (k = 0; k < p; k++) {
        for (j = 0; j < q; j++)
            explinstat[j * p + k] = swx[k] * expinf[j];
    }

    if (sw <= 1.0) 
        error("C_ExpCovLinstat: sum of weights is less than one");

    /* 
    *   covlinstat:  covariance of the linear statistic T
    */

    dsw = (double) sw;
    f1 = dsw / (dsw - 1);
    f2 = (1 / (dsw - 1));

    if (pq == 1) {
        covlinstat[0] = f1 * covinf[0] * swx2[0];
        covlinstat[0] -= f2 * covinf[0] * swx[0] * swx[0];
    } else {
        /* two more helpers needed */
        CT2 = Calloc(pq * pq, double);            /* pq x pq */
        Covy_x_swx = Calloc(pq * q, double);      /* pq x q  */
        
        C_kronecker(covinf, q, q, swx2, p, p, covlinstat);
        C_kronecker(covinf, q, q, swx, p, 1, Covy_x_swx);
        C_kronecker(Covy_x_swx, pq, q, swx, 1, p, CT2);

        for (k = 0; k < (pq * pq); k++)
            covlinstat[k] = f1 * covlinstat[k] - f2 * CT2[k];

        /* clean up */
        Free(CT2); Free(Covy_x_swx);
    }
}

/**
    R-interface to C_ExpCovLinstat\n
    *\param x values of the transformation
    *\param y values of the influence function
    *\param weights case weights
    *\param expcovinf an object of class `ExpectCovarInfluence'
*/

void C_LinstatExpCov (const SEXP x, const SEXP y, const SEXP weights, 
                      int *thisweights, SEXP ans) {
    
    SEXP expcovinf, explinstat, covlinstat, linstat, dim, cx, test;
    double *expinf, *covinf, *swx, *swx2, *dy, *dx;
    int sw = 0, i, n, p, q, pq;
    int *ix;

    /* determine the dimensions and some checks */

    n  = NROW(x);
    q  = NCOL(y);
    dy = REAL(y);

    expinf = Calloc(q, double);
    covinf = Calloc(q * q, double);

    if (NROW(y) != n)
        error("R_ExpCovLinstat: y does not have %d rows", n);
    if (LENGTH(weights) != n) 
        error("R_ExpCovLinstat: vector of case weights does not have %d elements", n);

    if (isReal(x)) {
        p = 1;
        pq = p * q;

        swx = Calloc(p, double);
        swx2 = Calloc(p, double);
        dx = REAL(x);
        SET_VECTOR_ELT(ans, 1, linstat = allocVector(REALSXP, pq));
        
        NA_weights_double(INTEGER(weights), dx, n, thisweights, &sw);
        C_swx_numeric(dx, thisweights, n, swx, swx2);
        C_Linstat_numeric(dx, dy, q, thisweights, n, REAL(linstat));
    } 
    if (isFactor(x) && !isOrdered(x)) { 
        p = C_nlevels(x);
        pq = p * q;

        swx = Calloc(p, double);
        swx2 = Calloc(p * p, double);
        ix = INTEGER(x);
        SET_VECTOR_ELT(ans, 1, linstat = allocVector(REALSXP, pq));

        NA_weights_factor(INTEGER(weights), ix, n, thisweights, &sw);                        
        C_swx_factor(ix, p, thisweights, n, swx, swx2);
        C_Linstat_factor(ix, p, dy, q, thisweights, n, REAL(linstat));
    }
    /*  FIXME: implement scores for ordered */
    if (isInteger(x) || isOrdered(x)) {
        p = 1;
        pq = p * q;

        swx = Calloc(p, double);
        swx2 = Calloc(p, double);
        PROTECT(cx = coerceVector(x, REALSXP));
        dx = REAL(cx);
        SET_VECTOR_ELT(ans, 1, linstat = allocVector(REALSXP, pq));

        NA_weights_double(INTEGER(weights), dx, n, thisweights, &sw);                        
        C_swx_numeric(dx, thisweights, n, swx, swx2);
        C_Linstat_numeric(dx, dy, q, thisweights, n, REAL(linstat));
        UNPROTECT(1);
    }

    SET_VECTOR_ELT(ans, 0, dim = allocVector(INTSXP, 2));    
    INTEGER(dim)[0] = p;
    INTEGER(dim)[1] = q;
    SET_VECTOR_ELT(ans, 2, explinstat = allocVector(REALSXP, pq));
    SET_VECTOR_ELT(ans, 3, covlinstat = allocVector(REALSXP, pq * pq));

    C_ExpInf(dy, q, thisweights, sw, n, expinf);
    C_CovInf(dy, q, thisweights, sw, n, expinf,  covinf);
    C_ExpCovLinstat(swx, swx2, p, q, sw, expinf, covinf, REAL(explinstat), 
                    REAL(covlinstat));

    Free(swx);
    Free(swx2);
    Free(expinf);
    Free(covinf);
}
