
#include "ctree.h"

void C_split_numeric (SEXP x, SEXP y, int *iweights, int minbucket, 
                      double *expinf, double *covinf, double *breaks) {

     int *orderx, i, j, k, n, q, sw = 0, count;
     double dsw = 0.0, *dx, lastx, thisx;
     double cstat = 0.0, laststat = 0.0;
     double *linstat, *explinstat, *covlinstat, *dy;

     n = NROW(x);
     if (NCOL(x) != 1 || !isReal(x))
         error("");

     /* orderx = order(x) */
     orderx = Calloc(n, int);
     for (i = 0; i < n; i++) orderx[i] = i;
     dx = Calloc(n, double);
     for (i = 0; i < n; i++) dx[i] = REAL(x)[i];
     rsort_with_index(dx, orderx, n);
     Free(dx);
     dx = REAL(x);

     q = NCOL(y);
     dy = REAL(y); 
     
     linstat = Calloc(q, double);
     explinstat = Calloc(q, double);
     covlinstat = Calloc(q * q, double);
     
     for (j = 0; j < q; j++) {
         linstat[j] = 0.0;
         explinstat[j] = 0.0;
         for (k = 0; k < q; k++) 
             covlinstat[j * q + k] = 0.0;
     }

     dsw = 0;
     lastx = NA_REAL;
     breaks[0] = NA_REAL;
     
     for (i = 0; i < n; i++) {
         if (iweights[i] == 0 || ISNA(dx[i])) continue;
         sw += iweights[i];
     }

     C_ExpInf(dy, q, iweights, sw, n, expinf);
     C_CovInf(dy, q, iweights, sw, n, expinf, covinf);
     
      /* check sample size contraint */
     if (sw >= minbucket) {
     for (i = 0; i < (n - 1); i++) {

         /* value of the orderx[i]th largest observation */
         thisx = dx[orderx[i]];
         
         /* nothing do to for zero weights or NAs */
         if (iweights[orderx[i]] == 0 || ISNA(thisx)) continue;

         /* sum of the weights of all x <= thisx */
         dsw += (double) iweights[orderx[i]];

         /* linear statistic for y ~ (x <= thisx) */
         for (k = 0; k < q; k++)
             linstat[k] += iweights[orderx[i]] * dy[k * n + orderx[i]];

         /* check sample size constraints */
         if (dsw > (sw - minbucket)) break;
         if (dsw < minbucket) continue;
         
         /* handle ties: look ahead if next valid x value is tied */
         count = 0;
         for (j = i + 1; j < (n - 1); j++) {
             count = j;
             if (iweights[orderx[j]] > 0) break;
         }
         if (thisx == dx[orderx[count]]) continue;
         lastx = thisx;

         /* compute expectation and covariance for linstat */
         C_ExpCovLinstat(&dsw, &dsw, 1, q, sw, expinf, covinf, explinstat,
                         covlinstat);

         /* max(abs(linstat - explinstat) / sqrt(diag(covlinstat))))) */
         cstat = C_maxabsTestStatistic(linstat, explinstat, covlinstat, q, 0.0);

         /* check for new maximum test statistic and save break point */
         if (cstat > laststat) {
             breaks[0] = lastx;
             laststat = cstat;
         }
     }
     }

     Free(orderx);
     Free(linstat);
     Free(explinstat);
     Free(covlinstat);
}


void C_split_factor (SEXP x, SEXP y, int *iweights, int minbucket, 
                     double *expinf, double *covinf, int *ans) {

    int i, jj, mi, l, j, p, q, *indl, n, k, sw = 0;
    int *ix;
    double *linstat, *explinstat, *covlinstat, swx, *dy;
    double cstat = 0.0, laststat = 0.0;
    
    n = NROW(x);
    ix = INTEGER(x);
    p = C_nlevels(x);
    q = NCOL(y);
    dy = REAL(y);

    for (i = 0; i < n; i++) {
        if (iweights[i] == 0 || ix[i] == NA_INTEGER) continue;
        sw += iweights[i];
    }

    for (l = 1; l < p; l++) ans[l] = NA_INTEGER;

    if (p == 2) {

        swx = 0.0;
        for (i = 0; i < n; i++) {
            if (iweights[i] == 0 || ix[i] == NA_INTEGER) continue;
            if (ix[i] == 1)
                swx += iweights[i];
        }

        /* check sample size contraint */
        if (!(swx > (sw - minbucket)) & !(swx < minbucket)) {
            ans[0] = 1;
            ans[1] = 2;
        }
    } else {

        /* check sample size contraint */
        if (sw >= minbucket) {

        C_ExpInf(dy, q, iweights, sw, n, expinf);
        C_CovInf(dy, q, iweights, sw, n, expinf, covinf);

        laststat = 0.0;
        mi = 1;  
        indl = Calloc(p, int);
        linstat = Calloc(q, double);
        explinstat = Calloc(q, double);
        covlinstat = Calloc(q * q, double);

        /* number of possible binary splits */
        for (l = 1; l < p; l++) mi *= 2;

        for (j = 1; j < mi; j++) { /* go though all splits */
            
            /* indl determines if level k is left or right */
            jj = j;
            for (l = 1; l < p; l++) {
                indl[l] = (jj%2);
                jj /= 2;
            }
            swx = 0.0;
            for (k = 0; k < q; k++)
                linstat[k] = 0.0;
        
            /* compute linear statistic and sum of weights */
            for (i = 0; i < n; i++) {
                if (iweights[i] > 0 && indl[ix[i] - 1]) {
                    swx += (double) iweights[i];
                    for (k = 0; k < q; k++)
                        linstat[k] += iweights[i] * dy[k * n + i];
                }
            }

            /* check sample size constraints */
            if (swx > (sw - minbucket)) continue;
            if (swx < minbucket) continue;

            /* compute expectation and covariance for linstat */
            C_ExpCovLinstat(&swx, &swx, 1, q, sw, expinf, covinf, explinstat, covlinstat);

            /* max(abs(linstat - explinstat) / sqrt(diag(covlinstat))))) */
            cstat = C_maxabsTestStatistic(linstat, explinstat, covlinstat, q, 0.0);

            /* check for new maximum test statistic and save break point */            
            if (cstat > laststat) {
                laststat = cstat;
                for (l = 0; l < p; l++)
                    ans[l] = indl[l] + 1;
            }
        }
        }

        Free(indl);
        Free(linstat);
        Free(explinstat);
        Free(covlinstat);
    }
}                                          
