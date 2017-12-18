
#include "rfweights.h"

SEXP R_rfweights (SEXP fdata, SEXP fnewdata, SEXP weights) {

    SEXP ans;
    int *id, *ind, *iweights, *ians;
    int Ntree = LENGTH(fdata), Ndata, Nnewdata;
    int OOB = LENGTH(fnewdata) == 0;

    if (TYPEOF(fdata) != VECSXP)
        error("fdata is not a list");    
    if (LENGTH(fdata) == 0)
        return(R_NilValue);
    if (TYPEOF(weights) != VECSXP)
        error("weights is not a list");    
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

    /* make this integer because ans can be really big */    
    PROTECT(ans = allocMatrix(INTSXP, Ndata, Nnewdata));
    ians = INTEGER(ans);
    
    for (int i = 0; i < Ndata * Nnewdata; i++)
        ians[i] = 0;
        
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
        for (int j = 0; j < Nnewdata; j++) {
            if (OOB & (iweights[j] > 0)) continue;
            for (int i = 0; i < Ndata; i++) {
                if (id[i] == ind[j]) 
                    ians[j * Ndata + i] += iweights[i];
            }
        }
    }
    
    UNPROTECT(1);
    return(ans);
}
