
.ctree_select <- function(...)
    function(model, trafo, data, subset, weights, whichvar, ctrl) {
        args <- list(...)
        ctrl[names(args)] <- args
        .select(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .ctree_test)
    }

.ctree_split <- function(...)
    function(model, trafo, data, subset, weights, whichvar, ctrl) {
        args <- list(...)
        ctrl[names(args)] <- args
        .split(model, trafo, data, subset, weights, whichvar, ctrl, FUN = .ctree_test)
    }

.ctree_test <- function(model, trafo, data, subset, weights, j, SPLITONLY = FALSE, ctrl) {

    ix <- data$zindex[[j]] ### data[[j, type = "index"]]
    iy <- data$yxindex ### data[["yx", type = "index"]]
    Y <- model$estfun

    if (!is.null(iy)) {
        stopifnot(NROW(levels(iy)) == (NROW(Y) - 1))
        return(.ctree_test_2d(data = data, j = j, Y = Y, iy = iy,
                              subset = subset, weights = weights,
                              SPLITONLY = SPLITONLY, ctrl = ctrl))
    }

    stopifnot(NROW(Y) == length(ix))

    NAyx <- data$yxmissings ### data[["yx", type = "missings"]]
    NAz <- data$missings[[j]] ### data[[j, type = "missings"]]
    if (ctrl$MIA && (ctrl$splittest || SPLITONLY)) {
        subsetNArm <- subset[!(subset %in% NAyx)]
    } else {
        subsetNArm <- subset[!(subset %in% c(NAyx, NAz))]
    }
    ### report by Kevin Ummel: _all_ obs being missing lead to
    ### subset being ignored completely
    if (length(subsetNArm) == 0) 
        return(list(statistic = NA, p.value = NA))

    return(.ctree_test_1d(data = data, j = j, Y = Y, subset = subsetNArm,
                          weights = weights, SPLITONLY = SPLITONLY, ctrl = ctrl))
}

.partysplit <- function(varid, breaks = NULL, index = NULL, right = TRUE,
                        prob = NULL, info = NULL) {
    ret <- list(varid = varid, breaks = breaks, index = index, right = right,
                prob = prob, info = info)
    class(ret) <- "partysplit"
    ret
}

.ctree_test_1d <- function(data, j, Y, subset, weights, SPLITONLY = FALSE, ctrl) {

    x <- data[[j]]
    MIA <- FALSE
    if (ctrl$MIA) {
        NAs <- data$missings[[j]] ### data[[j, type = "missings"]]
        MIA <- (length(NAs) > 0)
    }

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- data$zindex[[j]] ### data[[j, type = "index"]]

    scores <- data[[j, type = "scores"]]
    ORDERED <- is.ordered(x) || is.numeric(x)

    ux <- Xleft <- Xright <- NULL

    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        if (is.numeric(x)) {
            X <- data$zindex[[j]] ###data[[j, type = "index"]]
            ux <- levels(X)
        }
        if (MIA) {
            Xlev <- attr(X, "levels")
            Xleft <- X + 1L
            Xleft[NAs] <- 1L
            Xright <- X
            Xright[NAs] <- as.integer(length(Xlev) + 1L)
            attr(Xleft, "levels") <- c(NA, Xlev)
            attr(Xright, "levels") <- c(Xlev, NA)
        }
    } else {
        MAXSELECT <- FALSE
        if (is.numeric(x)) {
            if (storage.mode(x) == "double") {
                X <- x
            } else {
                X <- as.double(x) ### copy when necessary
            }
        }
        MIA <- FALSE
    }
    cluster <- data[["(cluster)"]]

    .ctree_test_internal(x = x, X = X, ix = NULL, Xleft = Xleft, Xright = Xright,
                         ixleft = NULL, ixright = NULL, ux = ux, scores = scores,
                         j = j, Y = Y, iy = NULL, subset = subset, weights = weights,
                         cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY,
                         MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}


.ctree_test_2d <- function(data, Y, iy, j, subset, weights, SPLITONLY = FALSE, ctrl) {

    x <- data[[j]]
    ix <- data$zindex[[j]] ### data[[j, type = "index"]]
    ux <- attr(ix, "levels")

    MIA <- FALSE
    if (ctrl$MIA) MIA <- any(ix[subset] == 0)

    ### X for (ordered) factors is always dummy matrix
    if (is.factor(x) || is.ordered(x))
        X <- integer(0)

    scores <- data[[j, type = "scores"]]
    ORDERED <- is.ordered(x) || is.numeric(x)

    if (ctrl$splittest || SPLITONLY) {
        MAXSELECT <- TRUE
        X <- integer(0)

        if (MIA) {
            Xlev <- attr(ix, "levels")
            ixleft <- ix + 1L
            ixright <- ix
            ixright[ixright == 0L] <- as.integer(length(Xlev) + 1L)
            attr(ixleft, "levels") <- c(NA, Xlev)
            attr(ixright, "levels") <- c(Xlev, NA)
            Xleft <- Xright <- X
        }
    } else {
        MAXSELECT <- FALSE
        MIA <- FALSE
        if (is.numeric(x))
            X <- matrix(c(0, as.double(attr(ix, "levels"))), ncol = 1)
    }
    cluster <- data[["(cluster)"]]

    .ctree_test_internal(x = x, X = X, ix = ix, Xleft = Xleft, Xright = Xright,
                         ixleft = ixleft, ixright = ixright, ux = ux, scores = scores,
                         j = j, Y = Y, iy = iy, subset = subset, weights = weights,
                         cluster = cluster, MIA = MIA, SPLITONLY = SPLITONLY,
                         MAXSELECT = MAXSELECT, ORDERED = ORDERED, ctrl = ctrl)
}



.ctree_test_internal <- function(x, X, ix, Xleft, Xright, ixleft, ixright, ux, scores, j, Y , iy,
                                 subset, weights, cluster, MIA, SPLITONLY, MAXSELECT, ORDERED, ctrl) {

    if (SPLITONLY) {
        nresample <- 0L
        varonly <- TRUE
        pvalue <- FALSE
        teststat <- ctrl$splitstat
    } else {
        nresample <- ifelse("MonteCarlo" %in% ctrl$testtype,
                        ctrl$nresample, 0L)
        pvalue <- !("Teststatistic" %in% ctrl$testtype)
        if (ctrl$splittest) {
            if (ctrl$teststat != ctrl$splitstat)
                warning("Using different test statistics for testing and splitting")
            teststat <- ctrl$splitstat
            if (nresample == 0 && pvalue)
               stop("MonteCarlo approximation mandatory for splittest = TRUE")
        } else {
           teststat <- ctrl$teststat
        }
        varonly <- "MonteCarlo" %in% ctrl$testtype &&
                   teststat == "maxtype"
    }

    ### see libcoin
    if (MAXSELECT) {
        if (!is.null(cluster)) varonly <- FALSE
    } else {
        if (is.ordered(x) && !ctrl$splittest)
            varonly <- FALSE
    }

    ### if (MIA) use tst as fallback
    ### compute linear statistic + expecation and covariance
    lev <- LinStatExpCov(X = X, Y = Y, ix = ix, iy = iy, subset = subset,
                         weights = weights, block = cluster,
                         nresample = nresample, varonly = varonly, checkNAs = FALSE)

    ### in some cases, estfun() might return NAs which we don't check
    if (any(is.na(lev$LinearStatistic))) {
        if (!is.null(iy)) {
            Ytmp <- Y[iy[subset] + 1L,]
        } else {
            Ytmp <- Y[subset,]
        }
        cc <- complete.cases(Ytmp)
        if (!all(cc)) { ### only NAs left
            if (SPLITONLY) return(NULL)
            return(list(statistic = NA, p.value = NA))
        }
        lev <- LinStatExpCov(X = X, Y = Y, ix = ix, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             nresample = nresample, varonly = varonly, checkNAs = TRUE)
    }

    if (!MAXSELECT) {
        if (is.ordered(x) && !ctrl$splittest)
            lev <- libcoin::lmult(matrix(scores, nrow = 1), lev)
    }

    ### check if either X or Y were unique
    if (all(lev$Variance < ctrl$tol)) {
        if (SPLITONLY) return(NULL)
        return(list(statistic = NA, p.value = NA))
    }

    ### compute test statistic and log(1 - p-value)
    tst <- doTest(lev, teststat = teststat, pvalue = pvalue,
                  lower = TRUE, log = TRUE, ordered = ORDERED,
                  maxselect = MAXSELECT,
                  minbucket = ctrl$minbucket, pargs = ctrl$pargs)

    if (MIA) {
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xleft, Y = Y, ix = ixleft, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             nresample = nresample, varonly = varonly, checkNAs = FALSE)
        ### compute test statistic and log(1 - p-value)
        tstleft <- doTest(lev, teststat = teststat, pvalue = pvalue,
                          lower = TRUE, log = TRUE, ordered = ORDERED,
                          minbucket = ctrl$minbucket, pargs = ctrl$pargs)
        ### compute linear statistic + expecation and covariance
        lev <- LinStatExpCov(X = Xright, Y = Y, ix = ixright, iy = iy, subset = subset,
                             weights = weights, block = cluster,
                             nresample = nresample, varonly = varonly, checkNAs = FALSE)
        ### compute test statistic and log(1 - p-value)
        tstright <- doTest(lev, teststat = teststat, pvalue = pvalue,
                           lower = TRUE, log = TRUE, ordered = ORDERED,
                           minbucket = ctrl$minbucket, pargs = ctrl$pargs)
    }

    if (!SPLITONLY) {
        if (MIA) {
            tst <- tstleft
            if (tst$TestStatistic < tstright$TestStatistic)
                tst <- tstright
        }
        return(list(statistic = log(pmax(tst$TestStatistic, .Machine$double.eps)),
                    p.value = tst$p.value))
    }

    ret <- NULL
    if (MIA && !any(is.na(tst$index))) {
        if (ORDERED) {
            if (tstleft$TestStatistic >= tstright$TestStatistic) {
                if (all(tst$index == 1)) { ### case C
                    ret <- .partysplit(as.integer(j), breaks = Inf,
                                      index = 1L:2L, prob = as.double(0:1))
                } else {
                    sp <- tstleft$index - 1L ### case A
                    if (!is.ordered(x)) {
                        ### interpolate split-points, see https://arxiv.org/abs/1611.04561
                        if (ctrl$intersplit & sp < length(ux)) {
                            sp <- (ux[sp] + ux[sp + 1]) / 2 ### <FIXME> use weighted mean here? </FIXME>
                        } else {
                            sp <- ux[sp]  ### X <= sp vs. X > sp
                        }
                    }
                    ret <- .partysplit(as.integer(j), breaks = sp,
                                      index = 1L:2L, prob = as.double(rev(0:1)))
                }
            } else {
                ### case C was handled above (tstleft = tstright in this case)
                sp <- tstright$index ### case B
                if (!is.ordered(x)) {
                    ### interpolate split-points, see https://arxiv.org/abs/1611.04561
                    if (ctrl$intersplit & sp < length(ux)) {
                        sp <- (ux[sp] + ux[sp + 1]) / 2
                    } else {
                        sp <- ux[sp]  ### X <= sp vs. X > sp
                    }
                }
                ret <- .partysplit(as.integer(j), breaks = sp,
                                  index = 1L:2L, prob = as.double(0:1))
            }
        } else {
            sp <- tstleft$index[-1L] ### tstleft = tstright for unordered factors
            if (length(unique(sp)) == 1L) { ### case C
                ret <- .partysplit(as.integer(j), index = as.integer(tst$index) + 1L)
            } else { ### always case A
                ret <- .partysplit(as.integer(j),
                                  index = as.integer(sp) + 1L,
                                  prob = as.double(rev(0:1)))
            }
        }
    } else {
        sp <- tst$index
        if (all(is.na(sp))) return(NULL)
        if (ORDERED) {
            if (!is.ordered(x))
                ### interpolate split-points, see https://arxiv.org/abs/1611.04561
                if (ctrl$intersplit & sp < length(ux)) {
                    sp <- (ux[sp] + ux[sp + 1]) / 2
                } else {
                    sp <- ux[sp]  ### X <= sp vs. X > sp
                }
                ret <- .partysplit(as.integer(j), breaks = sp,
                                  index = 1L:2L)
        } else {
            ret <- .partysplit(as.integer(j),
                              index = as.integer(sp) + 1L)
        }
    }
    return(ret)

}

ctree_control <- function
(
    teststat = c("quadratic", "maximum"),
    splitstat = c("quadratic", "maximum"), ### much better for q > 1, max was default
    splittest = FALSE,
    testtype = c("Bonferroni", "MonteCarlo",
                 "Univariate", "Teststatistic"),
    pargs = GenzBretz(),
    nmax = c("yx" = Inf, "z" = Inf),
    alpha = 0.05,
    mincriterion = 1 - alpha,
    logmincriterion = log(mincriterion),
    minsplit = 20L,
    minbucket = 7L,
    minprob = 0.01,
    stump = FALSE,
    lookahead = FALSE,	### try trafo() for daugther nodes before implementing the split
    MIA = FALSE,	### DOI: 10.1016/j.patrec.2008.01.010
    nresample = 9999L,
    tol = sqrt(.Machine$double.eps),
    maxsurrogate = 0L,
    numsurrogate = FALSE,
    mtry = Inf,
    maxdepth = Inf,
    multiway = FALSE,
    splittry = 2L,
    intersplit = FALSE,
    majority = FALSE,
    caseweights = TRUE,
    applyfun = NULL,
    cores = NULL,
    saveinfo = TRUE,
    update = NULL,
    splitflavour = c("ctree", "exhaustive")
) {

    testtype <- match.arg(testtype, several.ok = TRUE)
    if (length(testtype) == 4) testtype <- testtype[1]
    ttesttype <- testtype
    if (length(testtype) > 1) {
        stopifnot(all(testtype %in% c("Bonferroni", "MonteCarlo")))
        ttesttype <- "MonteCarlo"
    }

    if (MIA && maxsurrogate > 0)
        warning("Mixing MIA splits with surrogate splits does not make sense")

    if (MIA && majority)
        warning("Mixing MIA splits with majority does not make sense")

    splitstat <- match.arg(splitstat)
    teststat <- match.arg(teststat)

    if (!caseweights)
        stop("only caseweights currently implemented in ctree")

    splitflavour <- match.arg(splitflavour)

    c(extree_control(criterion = ifelse("Teststatistic" %in% testtype,
                                      "statistic", "p.value"),
                     logmincriterion = logmincriterion, minsplit = minsplit,
                     minbucket = minbucket, minprob = minprob,
                     nmax = nmax, stump = stump, lookahead = lookahead,
                     mtry = mtry, maxdepth = maxdepth, multiway = multiway,
                     splittry = splittry, maxsurrogate = maxsurrogate,
                     numsurrogate = numsurrogate,
                     majority = majority, caseweights = caseweights,
                     applyfun = applyfun, saveinfo = saveinfo,  ### always
                     selectfun = .ctree_select(),
                     splitfun = if (splitflavour == "ctree") .ctree_split() else .objfun_test(),
                     svselectfun = .ctree_select(),
                     svsplitfun =.ctree_split(minbucket = 0),
                     bonferroni = "Bonferroni" %in% testtype,
                     update = update),
      list(teststat = teststat, splitstat = splitstat, splittest = splittest, pargs = pargs,
           testtype = ttesttype, nresample = nresample, tol = tol,
           intersplit = intersplit, MIA = MIA))
}

ctree <- function(formula, data, subset, weights, na.action = na.pass, offset, cluster,
                  control = ctree_control(...), ytrafo = NULL, converged = NULL, scores = NULL,
                  doFit = TRUE, ...) {

    ## set up model.frame() call
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights",
                 "offset", "cluster", "scores"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$yx <- "none"
    if (is.function(ytrafo)) {
        if (all(c("y", "x") %in% names(formals(ytrafo))))
            mf$yx <- "matrix"
    }
    mf$nmax <- control$nmax
    ## evaluate model.frame
    mf[[1L]] <- quote(partykit::extree_data)

    d <- eval(mf, parent.frame())
    subset <- .start_subset(d)

    weights <- model.weights(model.frame(d))

    if (is.function(ytrafo)) {
        if (is.null(control$update))
            control$update <- TRUE
        nf <- names(formals(ytrafo))
        if (all(c("data", "weights", "control") %in% nf))
            ytrafo <- ytrafo(data = d, weights = weights, control = control)
        nf <- names(formals(ytrafo))
        stopifnot(all(c("subset", "weights", "info", "estfun", "object") %in% nf) ||
                  all(c("y", "x", "weights", "offset", "start") %in% nf))
    } else {
        if (is.null(control$update))
            control$update <- FALSE
        stopifnot(length(d$variables$x) == 0)
        mfyx <- model.frame(d, yxonly = TRUE)
        mfyx[["(weights)"]] <- mfyx[["(offset)"]] <- NULL
        yvars <- names(mfyx)
        for (yvar in yvars) {
            sc <- d[[yvar, "scores"]]
            if (!is.null(sc))
                attr(mfyx[[yvar]], "scores") <- sc
        }
        Y <- .y2infl(mfyx, response = d$variables$y, ytrafo = ytrafo)
        if (!is.null(iy <- d[["yx", type = "index"]])) {
            Y <- rbind(0, Y)
        }
        ytrafo <- function(subset, weights, info, estfun, object, ...)
            list(estfun = Y, unweighted = TRUE)
            ### unweighted = TRUE prevents estfun / w in extree_fit
    }
    if (is.function(converged)) {
        stopifnot(all(c("data", "weights", "control") %in% names(formals(converged))))
        converged <- converged(d, weights, control = control)
    } else {
        converged <- TRUE
    }

    update <- function(subset, weights, control, doFit = TRUE)
        extree_fit(data = d, trafo = ytrafo, converged = converged, partyvars = d$variables$z,
                   subset = subset, weights = weights, ctrl = control, doFit = doFit)
    if (!doFit) return(list(d = d, update = update))
    tree <- update(subset = subset, weights = weights, control = control)
    trafo <- tree$trafo
    tree <- tree$nodes

    mf <- model.frame(d)
    if (is.null(weights)) weights <- rep(1, nrow(mf))

    fitted <- data.frame("(fitted)" = fitted_node(tree, mf),
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = mf, fitted = fitted,
                 info = list(call = match.call(), control = control))
    ret$update <- update
    ret$trafo <- trafo
    class(ret) <- c("constparty", class(ret))

    ### doesn't work for Surv objects
    # ret$terms <- terms(formula, data = mf)
    ret$terms <- d$terms$all
    ### need to adjust print and plot methods
    ### for multivariate responses
    ### if (length(response) > 1) class(ret) <- "party"
    return(ret)
}


.logrank_trafo <- function(...)
    return(coin::logrank_trafo(...))

### convert response y to influence function h(y)
.y2infl <- function(data, response, ytrafo = NULL) {

    if (length(response) == 1) {
        if (!is.null(ytrafo[[response]])) {
            yfun <- ytrafo[[response]]
            rtype <- "user-defined"
        } else {
            rtype <- class(data[[response]])[1]
            if (rtype == "integer") rtype <- "numeric"
            if (rtype == "AsIs") rtype <- "numeric"
        }
        response <- data[[response]]

        infl <- switch(rtype,
            "user-defined" = yfun(response),
            "factor" = {
                X <- model.matrix(~ response - 1)
                if (nlevels(response) > 2) return(X)
                return(X[,-1, drop = FALSE])
            },
            "ordered" = {
                sc <- attr(response, "scores")
                if (is.null(sc)) sc <- 1L:nlevels(response)
                sc <- as.numeric(sc)
                return(matrix(sc[as.integer(response)], ncol = 1))
            },
            "numeric" = response,
            "Surv" = .logrank_trafo(response)
        )
    } else {
        ### multivariate response
        infl <- lapply(response, .y2infl, data = data)
        tmp <- do.call("cbind", infl)
        attr(tmp, "assign") <- rep(1L:length(infl), sapply(infl, NCOL))
        infl <- tmp
    }
    if (!is.matrix(infl)) infl <- matrix(infl, ncol = 1)
    storage.mode(infl) <- "double"
    return(infl)
}

sctest.constparty <- function(object, node = NULL, ...)
{

    if(is.null(node)) {
        ids <- nodeids(object, terminal = FALSE) ### all nodes
    } else {
        ids <- node
    }

    rval <- nodeapply(object, ids, function(n) {
        crit <- info_node(n)$criterion
        if (is.null(crit)) return(NULL)
        ret <- crit[c("statistic", "p.value"),,drop = FALSE]
        ret
    })
    names(rval) <- ids
    if(length(ids) == 1L)
        return(rval[[1L]])
    return(rval)
}
