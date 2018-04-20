
#include "rfweights.h"

SEXP R_rfweights (SEXP fdata, SEXP fnewdata, SEXP weights, SEXP scale) {

    SEXP ans;
    double *dans;
    int *id, *ind, *iweights, *tnsize;
    int Ntree = LENGTH(fdata), Ndata, Nnewdata;
    int OOB = LENGTH(fnewdata) == 0;

    if (TYPEOF(fdata) != VECSXP)
        error("fdata is not a list");    
    if (LENGTH(fdata) == 0)
        return(R_NilValue);
    if (TYPEOF(weights) != VECSXP)
        error("weights is not a list");    
    if (TYPEOF(scale) != LGLSXP || LENGTH(scale) != 1)
        error("scale is not a scalar logical");    
    if (LENGTH(weights) == 0)
        return(R_NilValue);
    Ndata = LENGTH(VECTOR_ELT(fdata, 0));
    if (OOB) {
        Nnewdata = Ndata;
        fnewdata = fdata;
    } else {
        if (LENGTH(fnewdata) == 0)
            return(R_NilValue);
        if (TYPEOF(fnewdata) != VECSXP)
            error("fnewdata is not a list");    
        Nnewdata = LENGTH(VECTOR_ELT(fnewdata, 0));
    }

    PROTECT(ans = allocMatrix(REALSXP, Ndata, Nnewdata));
    dans = REAL(ans);
    
    for (int i = 0; i < Ndata * Nnewdata; i++)
        dans[i] = 0.0;
        
    /* sum of weights for each terminal node id
       because trees can be very large (terminal node size = 1)
       we only once allocate Ndata integers */
    tnsize = Calloc(Ndata, int);
    for (int i = 0; i < Ndata; i++)
        tnsize[i] = 1;
        
    for (int b = 0; b < Ntree; b++) {
        if (TYPEOF(VECTOR_ELT(weights, b)) != INTSXP)
            error("some elements of weights are not integer");
        if (LENGTH(VECTOR_ELT(weights, b)) != Ndata)
            error("some elements of weights have incorrect length");
        if (TYPEOF(VECTOR_ELT(fnewdata, b)) != INTSXP)
            error("some elements of fnewdata are not integer");
        if (LENGTH(VECTOR_ELT(fnewdata, b)) != Nnewdata)
            error("some elements of fnewdata have incorrect length");
        if (TYPEOF(VECTOR_ELT(fdata, b)) != INTSXP)
            error("some elements of fdata are not integer");
        if (LENGTH(VECTOR_ELT(fdata, b)) != Ndata)
            error("some elements of fdata have incorrect length");
        iweights = INTEGER(VECTOR_ELT(weights, b));
        ind = INTEGER(VECTOR_ELT(fnewdata, b));
        id = INTEGER(VECTOR_ELT(fdata, b));
        
        if (LOGICAL(scale)[0]) {
            /* reset to zero */
            for (int i = 0; i < Ndata; i++)
                tnsize[i] = 0;
            /* sum of weights in each terminal node id */
            for (int i = 0; i < Ndata; i++)
                tnsize[id[i] - 1] += iweights[i];
        } /* else: tnsize == 1 */
        
        for (int j = 0; j < Nnewdata; j++) {
            if (OOB & (iweights[j] > 0)) continue;
            for (int i = 0; i < Ndata; i++) {
                if (id[i] == ind[j]) 
                    dans[j * Ndata + i] += (double) iweights[i] / tnsize[id[i] - 1];
            }
        }
    }
    
    Free(tnsize);
    
    UNPROTECT(1);
    return(ans);
}
