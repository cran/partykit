
logLik.constparty <- function(object, newdata, weights, perm = NULL, ...) {

    y <- object$fitted[["(response)"]]
    if (missing(newdata)) {
        fitted <- if (is.null(perm)) {
            object$fitted[["(fitted)"]]
        } else {
            ### no need to watch vmatch because newdata is always mf
            if (!is.null(perm)) {
                vnames <- names(object$data)
                if (is.character(perm)) {
                    stopifnot(all(perm %in% vnames))
                    perm <- match(perm, vnames)
                } else {
                    ### perm is a named list of factors coding strata
                    ### (for varimp(..., conditional = TRUE)
                    stopifnot(all(names(perm) %in% vnames))
                    stopifnot(all(sapply(perm, is.factor)))
                    tmp <- vector(mode = "list", length = length(vnames))
                    tmp[match(names(perm), vnames)] <- perm
                    perm <- tmp
                }
            }
            fitted_node(node_party(object), data = object$data, perm = perm)
        }
        pr <- predict_party(object, id = fitted, newdata = object$data,
                            type = ifelse(inherits(y, "factor"), "prob", "response"), ...)
    } else {
        pr <- predict(object, newdata = newdata, 
                      type = ifelse(inherits(y, "factor"), "prob", "response"), ...)
    }
    ll <- switch(class(y)[1], 
           "integer" = {
               -(y - pr)^2
           },
           "numeric" = {
               -(y - pr)^2
           },
          "factor" = {
              log(pmax(pr[cbind(1:length(y), unclass(y))], 
                        sqrt(.Machine$double.eps)))
          },
          "ordered" = {
              log(pmax(pr[cbind(1:length(y), unclass(y))], 
                       sqrt(.Machine$double.eps)))
          },
          "Surv" = stop("not yet implemented"),
          stop("not yet implemented")   
    )

    if (missing(weights)) weights <- data_party(object)[["(weights)"]]
    if (is.null(weights)) return(sum(ll) / length(y))
    return(sum(weights * ll) / sum(weights))
}

miscls <- function(object, newdata, weights, perm = NULL, ...) {

    y <- object$fitted[["(response)"]]
    stopifnot(is.factor(y))
    if (missing(newdata)) {
        fitted <- if (is.null(perm)) {
            object$fitted[["(fitted)"]]
        } else {
            ### no need to watch vmatch because newdata is always mf
            if (!is.null(perm)) {
                vnames <- names(object$data)
                if (is.character(perm)) {
                    stopifnot(all(perm %in% vnames))
                    perm <- match(perm, vnames)
                } else {
                    ### perm is a named list of factors coding strata
                    ### (for varimp(..., conditional = TRUE)
                    stopifnot(all(names(perm) %in% vnames))
                    stopifnot(all(sapply(perm, is.factor)))
                    tmp <- vector(mode = "list", length = length(vnames))
                    tmp[match(names(perm), vnames)] <- perm
                    perm <- tmp
                }
            }
            fitted_node(node_party(object), data = object$data, perm = perm)
        }
        pr <- predict_party(object, id = fitted, newdata = object$data,
                            type = "response", ...)
    } else {
        pr <- predict(object, newdata = newdata, type = "response", ...)
    }
    ll <- unclass(y) != unclass(pr)

    if (missing(weights)) weights <- data_party(object)[["(weights)"]]
    if (is.null(weights)) return(sum(ll) / length(y))
    return(sum(weights * ll) / sum(weights))
}

varimp <- function(object, nperm = 1L, ...)
    UseMethod("varimp")

varimp.constparty <- function(object, nperm = 1L, risk = c("loglik", "misclassification"), 
                              conditions = NULL, mincriterion = 0, ...) {

    if (!is.function(risk)) {
        risk <- match.arg(risk)
        ### risk is _NEGATIVE_ log-likelihood
        risk <- switch(risk, "loglik" = function(...) -logLik(...),
                             "misclassification" = miscls)
    }

    if (mincriterion > 0) 
        stop("mincriterion not yet implemented") ### use nodeprune

    psplitids <- unique(do.call("c", 
        nodeapply(node_party(object), 
                  ids = nodeids(node_party(object)),
                  FUN = function(x) split_node(x)$varid)))
    vnames <- names(object$data)
    psplitvars <- vnames[psplitids]
    ret <- numeric(length(psplitvars))
    names(ret) <- psplitvars

    for (vn in psplitvars) {
        cvn <- conditions[[vn]]
        if (is.null(cvn)) {
            perm <- vn
        } else {
            blocks <- .get_psplits(object, cvn) 
            if (length(blocks) == 0) blocks <- factor(rep(1, nrow(object$data)))
            perm <- vector(mode = "list", length = 1)
            names(perm) <- vn
            perm[[vn]] <- blocks
           }
        for (p in 1:nperm)
            ret[vn] <- ret[vn] + risk(object, perm = perm, ...)
    }
    ret <- ret / nperm - risk(object, ...)

    ret
}

gettree <- function(object, tree = 1L, ...)
    UseMethod("gettree")

gettree.cforest <- function(object, tree = 1L, ...) {
    ft <- object$fitted
    ft[["(weights)"]] <- object$weights[[tree]]
    ret <- party(object$nodes[[tree]], data = object$data, fitted = ft)
    ret$terms <- object$terms
    class(ret) <- c("constparty", class(ret))
    ret
}

.create_cond_list <- function(object, threshold) {

    d <- object$data
    response <- names(d)[attr(object$terms, "response")]
    xnames <- all.vars(object$terms)
    xnames <- xnames[xnames != response]

    ret <- lapply(xnames, function(x) {
        tmp <- ctree(as.formula(paste(x, "~", paste(xnames[xnames != x], collapse = "+"))),
                     data = d, control = ctree_control(teststat = "quad", testtype = "Univariate",
                                                       stump = TRUE))
        pval <- info_node(node_party(tmp))$criterion["p.value",]
        pval[is.na(pval)] <- 1
        ### make the meaning of threshold equal to partykit
        ret <- names(pval)[(1 - pval) > threshold] 
        if (length(ret) == 0) return(NULL) 
        return(ret)
    })
    names(ret) <- xnames
    return(ret)
}

.get_psplits <- function(object, xnames) {

    d <- object$data
    ret <- lapply(xnames, function(xn) {
        id <- which(colnames(d) == xn)
        psplits <- nodeapply(node_party(object), 
            ids = nodeids(node_party(object)),
            FUN = function(x) {
                if (is.null(x)) return(NULL)
                if (is.terminal(x)) return(NULL)
                if (split_node(x)$varid == id)
                    return(split_node(x))
                return(NULL)
            })
        psplits <- psplits[!sapply(psplits, is.null)]
        if (length(psplits) > 0)
            return(do.call("interaction", lapply(psplits, kidids_split, data = d))[, drop = TRUE])
        return(NULL)
    })
    ret <- ret[!sapply(ret, is.null)]
    if (length(ret) > 0) {
        if (length(ret) == 1) return(factor(ret[[1]], exclude = NULL))
        ### get rid of empty levels quickly; do.call("interaction", ret)
        ### explodes
        for (i in 2:length(ret))
            ret[[1]] <- factor(interaction(ret[[1]], ret[[i]])[, drop = TRUE], exclude = NULL)
        return(ret[[1]])
    }
    return(NULL)
}

varimp.cforest <- function(object, nperm = 1L, OOB = TRUE, risk = c("loglik", "misclassification"), 
                           conditional = FALSE, threshold = .2,    
                           applyfun = NULL, cores = NULL, ...) {

    ret <- matrix(NA, nrow = length(object$nodes), ncol = ncol(object$data))
    colnames(ret) <- names(object$data)

    if (conditional) {
        conditions <- .create_cond_list(object, threshold)
    } else {
        conditions <- NULL
    }

    ## apply infrastructure 
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply  
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.set.seed = TRUE, mc.cores = cores)
        }
    }

    vi <- applyfun(1:length(object$nodes), function(b) {
        tree <- gettree(object, b)
        if (OOB) {
            oobw <- as.integer(object$weights[[b]] == 0)
            vi <- varimp(tree, nperm = nperm, risk = risk, conditions = conditions, 
                         weights = oobw, ...)
        } else {
            vi <- varimp(tree, nperm = nperm, risk = risk, conditions = conditions, 
                         ...)
        }
        return(vi)
    })

    for (b in 1:length(object$nodes))
        ret[b, match(names(vi[[b]]), colnames(ret))] <- vi[[b]]

    ret <- colMeans(ret, na.rm = TRUE)
    ret[!sapply(ret, is.na)]
}
