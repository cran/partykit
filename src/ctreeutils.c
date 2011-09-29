
#include "ctree.h"

void C_abs_double (double *x, int n) {

    int i;
    for (i = 0; i < n; i++) x[i] = fabs(x[i]);
}
        
double C_max_double (const double *x, const int n) {

   double max = 0.0;
   int i;
         
   for (i = 0; i < n; i++) {
       if (x[i] > max) max = x[i];
   }
   return(max);
}

int C_nlevels (SEXP x) {
   return(LENGTH(getAttrib(x, R_LevelsSymbol)));
}

int NROW (SEXP x) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(LENGTH(x));
    return(INTEGER(a)[0]);
}
    
int NCOL (SEXP x) {
    SEXP a;
    a = getAttrib(x, R_DimSymbol);
    if (a == R_NilValue) return(1);
    return(INTEGER(a)[1]);
}

/**
    Computes the Kronecker product of two matrices\n
    *\param A matrix
    *\param m nrow(A)
    *\param n ncol(A)
    *\param B matrix
    *\param r nrow(B)
    *\param s ncol(B)
    *\param ans return value; a pointer to a REALSXP-vector of length (mr x ns)
*/

void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans) {

    int i, j, k, l, mr, js, ir;
    double y;

    mr = m * r;
    for (i = 0; i < m; i++) {
        ir = i * r;
        for (j = 0; j < n; j++) {
            js = j * s;
            y = A[j*m + i];
            for (k = 0; k < r; k++) {
                for (l = 0; l < s; l++) {
                    ans[(js + l) * mr + ir + k] = y * B[l * r + k];
                }
            }
        }
    }
}  

void NA_weights_double (const int *weights, const double *x, const int n, 
                        int *thisweights, int *sw) {

    int i;
    
    sw[0] = 0;
    for (i = 0; i < n; i++) {
        if (ISNA(x[i])) {
            thisweights[i] = 0;
        } else {
            thisweights[i] = weights[i];
            sw[0] += weights[i];
        }
    }
}

void NA_weights_factor (const int *weights, const int *x, const int n, 
                        int *thisweights, int *sw) {

    int i;
    
    sw[0] = 0;
    for (i = 0; i < n; i++) {
        if (x[i] == NA_INTEGER) {
            thisweights[i] = 0;
        } else {
            thisweights[i] = weights[i];
            sw[0] += weights[i];
        }
    }
}
