
/**
    utility functions.
    *\file utils.c
    *\author $Author$
    *\date $Date$
*/

#include "party.h"

/**
    determine in which of the n intervals given by  
    (-Inf, breaks[0]], (breaks[0], breaks[1]], ..., (breaks[n - 1], +Inf)
    the observation x falls (intervals are numbered 0, ..., n).
    *\param x observation
    *\param breaks ordered cutpoints
    *\param n length(breaks)
    *\param right logical: TRUE means intervals are closed on the right, 
                  FALSE closed on the left
*/

int cut(double x, double *breaks, int n, int right) {

    int ret, i;

    ret = NA_INTEGER;
    if (ISNA(x)) return(ret);

    /* x falls in last intervall */
    if (x > breaks[n - 1]) {
        ret = n;
    } else {
        for (i = 0; i < n; i++) {
            if (x <= breaks[i]) {
                ret = i;
                break;
            }
        }
        /* intervals are closed on the left */
        if (!right)
            if (x == breaks[ret]) ret++;
    }
    /* ret is in 0, ..., n */
    return(ret);
}

void C_SampleNoReplace(int *x, int m, int k, int *ans) {
     
     
    int i, j, n = m;
         
    for (i = 0; i < m; i++)
        x[i] = i;
    for (i = 0; i < k; i++) {
        j = n * unif_rand(); 
        ans[i] = x[j];
        x[j] = x[--n];
    }
}
