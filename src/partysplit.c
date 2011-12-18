
/**
    infrastructure for splits
    *\file utils.c
    *\author $Author$
    *\date $Date$
*/

#include "party.h"
#include "partysplit.h"


void init_partysplit(SEXP varid, SEXP breaks, SEXP index, SEXP right,
             SEXP prob, SEXP info, SEXP split) {
             
    if (LENGTH(split) != LENGTH_SPLIT)
        error("split is not a list with %d elements", LENGTH_SPLIT);
        
    if (isInteger(varid))
        SET_VECTOR_ELT(split, VARID_SPLIT, varid);

    /* FIXME: sort first? */
    SET_VECTOR_ELT(split, BREAKS_SPLIT, coerceVector(breaks, REALSXP));

    SET_VECTOR_ELT(split, INDEX_SPLIT, index);
    SET_VECTOR_ELT(split, RIGHT_SPLIT, right);
    SET_VECTOR_ELT(split, PROB_SPLIT, prob);
    SET_VECTOR_ELT(split, INFO_SPLIT, info);
}

int varid_split(SEXP split) {

    return(INTEGER(VECTOR_ELT(split, VARID_SPLIT))[0]);
}

SEXP breaks_split(SEXP split) {

    return(VECTOR_ELT(split, BREAKS_SPLIT));
}

SEXP index_split(SEXP split) {

    return(VECTOR_ELT(split, INDEX_SPLIT));
}

int right_split(SEXP split) {

    return(LOGICAL(VECTOR_ELT(split, RIGHT_SPLIT))[0]);
}

/**
    extract or setup a probability distribution on the kids 
    for random splitting (in case of missing values)
    *\param split a split object
*/

SEXP prob_split(SEXP split) {

    SEXP prob, index, breaks;
    double sum = 0.0;
    int i;

    /* someone already specified the distribution */
    prob = VECTOR_ELT(split, PROB_SPLIT);
    if (prob != R_NilValue)
        return(prob);
        
    /* set it up from index */
    index = index_split(split);
    
    /* index not explicitly given */
    if (index == R_NilValue) {
        /* index = 1:(length(breaks) + 1) */
        breaks = breaks_split(split);
        if (breaks == R_NilValue)
            error("prob, index and breaks are missing");
        SET_VECTOR_ELT(split, INDEX_SPLIT, 
                       index = allocVector(INTSXP, LENGTH(breaks) + 1));
        for (i = 0; i <= LENGTH(breaks); i++)
            INTEGER(index)[i] = i + 1;
    }
    /* probability distribution with support index[!is.na(index)] */
    SET_VECTOR_ELT(split, PROB_SPLIT, 
                   prob = allocVector(REALSXP, LENGTH(index)));
    for (i = 0; i < LENGTH(index); i++) {
        REAL(prob)[i] = (double) (INTEGER(index)[i] != NA_INTEGER);
        sum += REAL(prob)[i];
    }
    for (i = 0; i < LENGTH(index); i++)
        REAL(prob)[i] = REAL(prob)[i] / sum;
    return(prob);
}

SEXP info_split(SEXP split) {

    return(VECTOR_ELT(split, INFO_SPLIT));
}

SEXP split_data(SEXP split, SEXP data, SEXP vmatch) {

    if (vmatch == R_NilValue)
        return(VECTOR_ELT(data, varid_split(split) - 1));
    return(VECTOR_ELT(data, INTEGER(vmatch)[varid_split(split) - 1] - 1));
}

double x2d(SEXP x, int obs) {

    double ret = NA_REAL;

    /* FIXME: is this enough checking? */
    if (isReal(x)) { ret = REAL(x)[obs];
    } else {
        if (INTEGER(x)[obs] != NA_INTEGER)
            ret = (double) INTEGER(x)[obs];
    }
    /* FIXME: take care of other storage modes (char for example) */
    /* if (ISNA(ret)) error("can't coerce x to REAL or INTEGER"); */
    return(ret);
}

/**
    determine the child node observation obs has to go into based
    on one split (either primary or surrogate)
    *\param split a split object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs observation number (starting with 0)
*/

int kidid_split(SEXP split, SEXP data, SEXP vmatch, int obs) {

    SEXP x, breaks;
    int ret = NA_INTEGER;

    /* get the variable (all observations) */
    x = split_data(split, data, vmatch);

    breaks = breaks_split(split);

    /* x is a factor and needs no breaks */    
    if (breaks == R_NilValue) {
        /* FIXME: is this enough checking? */
        if (isReal(x)) 
            error("x is not an integer or factor (variable %d) \n", varid_split(split));
        ret = INTEGER(x)[obs];
        if (ret != NA_INTEGER) ret = ret - 1;
    } else {
        /* determine number of interval */
        ret = cut(x2d(x, obs), REAL(breaks), LENGTH(breaks), 
                  right_split(split));
    }
    /* use index (if there) to assign child nodes */
    if (ret != NA_INTEGER) {
        if (index_split(split) != R_NilValue) {
           ret = INTEGER(index_split(split))[ret];
           if (ret != NA_INTEGER) ret = ret - 1;
        }
    }
    /* ret is in 0, ..., LENGTH(index) - 1 */
    return(ret);
}
