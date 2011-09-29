
void C_Linstat_numeric (const double *x, const double *y, 
                        const int q, const int *weights, const int n,
                        double *ans);
void C_Linstat_factor (const int *x, const int p,
                       const double *y, const int q,
                       const int *weights, const int n,
                       double *ans);
void C_ExpInf (const double* y, const int q, const int* weights, const int sw,
               const int n, double* ans);
void C_CovInf (const double* y, const int q, const int* weights, const int sw,
               const int n, const double* exp, double *ans);
void C_ExpCovLinstat (const double *swx, const double *swx2, const int p, const int q, 
                      const int sw, const double *expinf, const double *covinf, 
                      double *explinstat, double *covlinstat);
void C_LinstatExpCov(const SEXP x, const SEXP y, const SEXP weights, 
                     int *thisweights, SEXP ans);
