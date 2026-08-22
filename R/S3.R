
### check if methods were called directly, prepared in 1.3-0
chkS3 <- function(caller, callee) {

    ### hold fire
    return(NULL)

    if (is.null(callee)) return(NULL)

    callee <- rev(as.character(as.list(callee)[[1]]))[1L]
    genfun <- strsplit(callee, "\\.")[[1L]][1L]

    ### called directly from .Globalenv
    if (is.null(caller)) {
        .Deprecated(new = genfun, old = callee, package = "partykit",
                    msg = "calling partykit methods directly is deprecated, please call the generic")
        return(NULL)
    }

    callerchr <- paste(deparse(caller), sep = "--", collapse = "--")
    if (length(grep("UseMethod", callerchr)))
        return(NULL)

    caller <- rev(as.character(as.list(caller)[[1]]))[1L]
    if (caller != genfun && ### it wasn't the generic
        inherits(try(getFromNamespace(caller, ns = "partykit"), silent = TRUE), 
                 "try-error")) ### nor a partykit fct
        .Deprecated(new = genfun, old = callee, package = "partykit",
                    msg = "calling partykit methods directly is deprecated, please call the generic")
    return(NULL)
}
