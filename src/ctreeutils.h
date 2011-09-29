
void C_abs_double (double *x, int n);
double C_max_double (const double *x, const int n);
int C_nlevels (SEXP x);
int NROW (SEXP x);
int NCOL (SEXP x);
void C_kronecker (const double *A, const int m, const int n,
                  const double *B, const int r, const int s,
                  double *ans);
void NA_weights_double (const int *weights, const double *x, const int n, 
                        int *thisweights, int *sw);
void NA_weights_factor (const int *weights, const int *x, const int n, 
                        int *thisweights, int *sw);
