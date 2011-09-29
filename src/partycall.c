
/**
    .Call interfaces
    *\file utils.c
    *\author $Author$
    *\date $Date$
*/

#include "party.h"
#include "partysplit.h"
#include "partynode.h"

/**
    determine the child node observations obs has to go into based
    on one split (either primary or surrogate)
    *\param split a split object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs integer vector of observation numbers
*/
                
SEXP R_kidids_split(SEXP split, SEXP data, SEXP vmatch, SEXP obs) {

    SEXP ans;
    int i, tmp;   
        
    PROTECT(ans = allocVector(INTSXP, LENGTH(obs)));
    for (i = 0; i < LENGTH(ans); i++) {
        tmp = kidid_split(split, data, vmatch, INTEGER(obs)[i] - 1);
        if (tmp != NA_INTEGER) {
            INTEGER(ans)[i] = tmp + 1;
        } else {
            INTEGER(ans)[i] = NA_INTEGER;
        }
        
    }
    UNPROTECT(1);
    return(ans);
}

/**
    determine the terminal node id for observations obs.
    *\param node a node object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs integer vector of observation numbers
*/

SEXP R_fitted_node(SEXP node, SEXP data, SEXP vmatch, SEXP obs, SEXP perm) {

    SEXP ans, perms, p;
    int i, *tmp;

    /* we might want to do random splitting */
    GetRNGstate();

    if (perm != R_NilValue) {
        PROTECT(perms = allocVector(VECSXP, LENGTH(perm)));
        tmp = Calloc(LENGTH(obs), int);
        for (i = 0; i < LENGTH(perm); i++) {
            SET_VECTOR_ELT(perms, i, p = allocVector(INTSXP, LENGTH(obs)));
            C_SampleNoReplace(tmp, LENGTH(obs), LENGTH(obs), INTEGER(p));
        }
    }
         
    PROTECT(ans = allocVector(INTSXP, LENGTH(obs)));
    for (i = 0; i < LENGTH(ans); i++)
        INTEGER(ans)[i] = fitted_node(node, data, vmatch, perm, perms, 
                                      INTEGER(obs)[i] - 1);

    PutRNGstate();

    UNPROTECT(1);
    if (perm != R_NilValue) 
        UNPROTECT(1);
    return(ans);
}
