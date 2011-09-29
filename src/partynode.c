
/**
    infrastructure for nodes
    *\file utils.c
    *\author $Author$
    *\date $Date$
*/


#include "party.h"
#include "partynode.h"
#include "partysplit.h"

int id_node(SEXP node) {
    return(INTEGER(VECTOR_ELT(node, ID_NODE))[0]);
}

SEXP split_node(SEXP node) {
    return(VECTOR_ELT(node, SPLIT_NODE));
}

SEXP kids_node(SEXP node) {
    return(VECTOR_ELT(node, KIDS_NODE));
}

SEXP surrogates_node(SEXP node) {
    return(VECTOR_ELT(node, SURROGATES_NODE));
}

SEXP info_node(SEXP node) {
    return(VECTOR_ELT(node, INFO_NODE));
}

int is_terminal_node(SEXP node) {
    if (kids_node(node) == R_NilValue) {
        if (split_node(node) != R_NilValue)
            error("kids != split");
        return(1);
    }
    return(0);
}

/**
    determine the child node observation obs has to go into based
    on one all available splits (primary, surrogate, random).
    *\param node a node object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs observation number (starting with 0)
*/
    
int kidid_node(SEXP node, SEXP data, SEXP vmatch, int obs) {

    SEXP primary, surrogates, prob;
    int ret = NA_INTEGER;
    int i;
    double *dprob;

    primary = split_node(node);
    surrogates = surrogates_node(node);

    /* perform primary split */
    ret = kidid_split(primary, data, vmatch, obs);

    /* surrogate / random splits if needed */
    if (ret == NA_INTEGER) {
        /* surrogate splits first */
        if (LENGTH(surrogates) >= 1) {
            for (i = 0; i < LENGTH(surrogates); i++) {
                if (ret != NA_INTEGER) break;
                ret = kidid_split(VECTOR_ELT(surrogates, i), data, vmatch, obs);
            }
        }

        /* random splits when necessary */
        if (ret == NA_INTEGER) {
            prob = prob_split(primary);
            /* the following equals sample(index, 1, prob) */
            dprob = Calloc(LENGTH(prob) - 1, double);
            dprob[0] = REAL(prob)[0];
            for (i = 1; i < LENGTH(prob) - 1; i++)
                dprob[i] = REAL(prob)[i] + dprob[i - 1];
            ret = cut(unif_rand(), dprob, LENGTH(prob) - 1, 1);
            Free(dprob);
        }
    }

    /* just in case ... */
    if (ret == NA_INTEGER)
        error("failed to predict kidid from node %d for observation %d\n", 
              id_node(node), obs);

    /* ret is in 0, ..., LENGTH(kidid_split(node)) - 1 */
    return(ret);
}

/**
    determine the terminal node id for observation obs.
    *\param node a node object
    *\param data a list
    *\param vmatch an integer for permuting variables
    *\param obs observation number (starting with 0)
*/

int fitted_node(SEXP node, SEXP data, SEXP vmatch, SEXP perm, SEXP perms, int obs) {

    int kidid, ret, i;
    
    /* return the id of the terminal node */
    if (is_terminal_node(node))
        return(id_node(node));
    
    /* permute obs */
    if (perm != R_NilValue) {
        for (i = 0; i < LENGTH(perm); i++) {
            if (varid_split(split_node(node)) == INTEGER(perm)[i])
                obs = INTEGER(VECTOR_ELT(perms, i))[obs];
        }
    }
    
    /* determine next child node */
    kidid = kidid_node(node, data, vmatch, obs);
    
    /* send observation down to node kidid (recursively!) */
    ret = fitted_node(VECTOR_ELT(kids_node(node), kidid), 
                      data, vmatch, perm, perms, obs);

    /* ret is in 1, ..., #nodes */
    return(ret);
}
