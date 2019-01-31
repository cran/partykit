
### constructor for forest objects
constparties <- function(nodes, data, weights, fitted = NULL, terms = NULL, info = NULL) {

    stopifnot(all(sapply(nodes, function(x) inherits(x, "partynode"))))
    stopifnot(inherits(data, "data.frame"))
    stopifnot(inherits(weights, "list"))

    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot(nrow(data) == 0L | nrow(data) == nrow(fitted))
        if (nrow(data) == 0L)
            stopifnot("(response)" %in% names(fitted))
    } else {
        stopifnot(nrow(data) > 0L)
        stopifnot(!is.null(terms))
        fitted <- data.frame("(response)" = model.response(model.frame(terms, data = data, 
                                                                       na.action = na.pass)),
                             check.names = FALSE)
    }

    ret <- list(nodes = nodes, data = data, weights = weights, fitted = fitted)
    class(ret) <- c("constparties", "parties")

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
        ret$terms <- terms
    }

    if (!is.null(info))
        ret$info <- info

    ret
}

.perturb <- function(replace = FALSE, fraction = .632) {
    ret <- function(prob) {
        if (replace) {
            rw <- rmultinom(1, size = length(prob), prob = prob)
        } else {
            rw <- integer(length(prob))
            i <- sample(1:length(prob), ceiling(fraction * length(prob)), prob = prob)
            rw[i] <- 1L
        }
        as.integer(rw)
    }
    ret
}

cforest <- function
(
    formula,
    data,   
    weights,
    subset, 
    offset, 
    cluster,
    strata,
    na.action = na.pass,
    control = ctree_control(
        teststat = "quad", testtype = "Univ", mincriterion = 0,
        saveinfo = FALSE, ...),
    ytrafo = NULL, 
    scores = NULL, 
    ntree = 500L, 
    perturb = list(replace = FALSE, fraction = 0.632),
    mtry = ceiling(sqrt(nvar)), 
    applyfun = NULL,
    cores = NULL, 
    trace = FALSE,
    ...
) {
   
    ### get the call and the calling environment for .urp_tree
    call <- match.call(expand.dots = FALSE)
    oweights <- NULL
    if (!missing(weights))
        oweights <- weights
    m <- match(c("formula", "data", "subset", "na.action", "offset", "cluster", 
                 "scores", "ytrafo", "control"), names(call), 0L)
    ctreecall <- call[c(1L, m)]
    ctreecall$doFit <- FALSE
    if (!is.null(oweights))
        ctreecall$weights <- 1:NROW(oweights)
    ctreecall$control <- control ### put ... into ctree_control()
    ctreecall[[1L]] <- quote(partykit::ctree)
    tree <- eval(ctreecall, parent.frame())

    if (is.null(control$update))
        control$update <- is.function(ytrafo)

    d <- tree$d
    updatefun <- tree$update

    nvar <- sum(d$variables$z > 0)
    control$mtry <- mtry
    control$applyfun <- lapply
 
    strata <- d[["(strata)"]]
    if (!is.null(strata)) {
        if (!is.factor(strata)) stop("strata is not a single factor")
    }
    
    probw <- NULL
    iweights <- model.weights(model.frame(d))
    if (!is.null(oweights)) {
        if (is.matrix(oweights)) {
            weights <- oweights[iweights,,drop = FALSE]
        } else {
            weights <- oweights[iweights]
        }
    } else {
        weights <- NULL
    }
    rm(oweights)
    rm(iweights)
    N <- nrow(model.frame(d))
    rw <- NULL
    if (!is.null(weights)) {
        if (is.matrix(weights)) {
            if (ncol(weights) == ntree && nrow(weights) == N) {
                rw <- unclass(as.data.frame(weights))
                rw <- lapply(rw, function(w) 
                    rep(1:length(w), w))
                weights <- integer(0)
            } else {
                stop(sQuote("weights"), "argument incorrect")
            }
        } else {
            probw <- weights / sum(weights)
        }
    } else {
        weights <- integer(0)
    }

    idx <- .start_subset(d)
    if (is.null(rw)) {
        if (is.null(strata)) {
            size <- N
            if (!perturb$replace) size <- floor(size * perturb$fraction)
            rw <- replicate(ntree, 
                            sample(idx, size = size, 
                                   replace = perturb$replace, prob = probw),
                            simplify = FALSE)
        } else {
            frac <- if (!perturb$replace) perturb$fraction else 1
            rw <- replicate(ntree, function() 
                  do.call("c", tapply(idx, strata, 
                          function(i) 
                              sample(i, size = length(i) * frac, 
                                     replace = perturb$replace, prob = probw[i]))))
        }
    }

    ## apply infrastructure for determining split points
    ## use RNGkind("L'Ecuyer-CMRG") to make this reproducible
    if (is.null(applyfun)) {
        applyfun <- if(is.null(cores)) {
            lapply  
        } else {
            function(X, FUN, ...)
                parallel::mclapply(X, FUN, ..., mc.set.seed = TRUE, mc.cores = cores)
        }
    }

    trafo <- updatefun(sort(rw[[1]]), integer(0), control, doFit = FALSE)
    if (trace) pb <- txtProgressBar(style = 3) 
    forest <- applyfun(1:ntree, function(b) {
        if (trace) setTxtProgressBar(pb, b/ntree)
        ret <- updatefun(sort(rw[[b]]), integer(0), control)
        # trafo <<- ret$trafo
        ret$nodes
    })
    if (trace) close(pb)

    fitted <- data.frame(idx = 1:N)  
    mf <- model.frame(d)
    fitted[[2]] <- mf[, d$variables$y, drop = TRUE]
    names(fitted)[2] <- "(response)"
    if (length(weights) > 0)
        fitted[["(weights)"]] <- weights

    ### turn subsets in weights (maybe we can avoid this?)
    rw <- lapply(rw, function(x) as.integer(tabulate(x, nbins = length(idx))))

    control$applyfun <- applyfun

    ret <- constparties(nodes = forest, data = mf, weights = rw,
                        fitted = fitted, terms = d$terms$all,
                        info = list(call = match.call(), control = control))
    ret$trafo <- trafo
    ret$predictf <- d$terms$z
    class(ret) <- c("cforest", class(ret))

    return(ret)
}

predict.cforest <- function(object, newdata = NULL, type = c("response", "prob", "weights", "node"), 
                            OOB = FALSE, FUN = NULL, simplify = TRUE, scale = TRUE, ...) {

    responses <- object$fitted[["(response)"]]
    forest <- object$nodes
    nd <- object$data
    vmatch <- 1:ncol(nd)
    NOnewdata <- TRUE
    if (!is.null(newdata)) {
        factors <- which(sapply(nd, is.factor))
        xlev <- lapply(factors, function(x) levels(nd[[x]]))
        names(xlev) <- names(nd)[factors]
        nd <- model.frame(object$predictf, ### all variables W/O response
                          data = newdata, na.action = na.pass, xlev = xlev)
        OOB <- FALSE
        vmatch <- match(names(object$data), names(nd))
        NOnewdata <- FALSE
    }
    nam <- rownames(nd)

    type <- match.arg(type)

    ### return terminal node ids for data or newdata
    if (type == "node")
        return(lapply(forest, fitted_node, data = nd, vmatch = vmatch, ...))

    ### extract weights
    rw <- object$weights

    w <- 0L

    applyfun <- lapply
    if (!is.null(object$info))
        applyfun <- object$info$control$applyfun

    fdata <- lapply(forest, fitted_node, data = object$data, ...)
    if (NOnewdata && OOB) {
        fnewdata <- list()
    } else {
        fnewdata <- lapply(forest, fitted_node, data = nd, vmatch = vmatch, ...)
    }

    w <- .rfweights(fdata, fnewdata, rw, scale)

#    for (b in 1:length(forest)) {
#        ids <- nodeids(forest[[b]], terminal = TRUE)
#        fnewdata <- fitted_node(forest[[b]], nd, vmatch = vmatch, ...)
#        fdata <- fitted_node(forest[[b]], object$data, ...)
#        tw <- rw[[b]]
#        pw <- sapply(ids, function(i) tw * (fdata == i))
#        ret <- pw[, match(fnewdata, ids), drop = FALSE]
#        ### obs which are in-bag for this tree don't contribute
#        if (OOB) ret[,tw > 0] <- 0
#        w <- w + ret
#    }
#
#    #w <- Reduce("+", bw)
#    if (!is.matrix(w)) w <- matrix(w, ncol = 1)

    if (type == "weights") {
        ret <- w
        colnames(ret) <- nam
        rownames(ret) <- rownames(responses)
        return(ret)
    }
    
    pfun <- function(response) {

        if (is.null(FUN)) {

            rtype <- class(response)[1]
            if (rtype == "ordered") rtype <- "factor"
            if (rtype == "integer") rtype <- "numeric"

            FUN <- switch(rtype,
                "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
                "factor" = if (type == "response") .pred_factor_response else .pred_factor,
                "numeric" = if (type == "response") .pred_numeric_response else .pred_ecdf)
        }

        ret <- vector(mode = "list", length = ncol(w))
        for (j in 1:ncol(w))
            ret[[j]] <- FUN(response, w[,j])
        ret <- as.array(ret)
        dim(ret) <- NULL
        names(ret) <- nam
         
        if (simplify)
            ret <- .simplify_pred(ret, names(ret), names(ret))
        ret
    }
    if (!is.data.frame(responses)) {
        ret <- pfun(responses)
    } else {
        ret <- lapply(responses, pfun)
        if (all(sapply(ret, is.atomic)))
            ret <- as.data.frame(ret)
        names(ret) <- colnames(responses)
    }
    ret
}

model.frame.cforest <- function(formula, ...) {
    class(formula) <- "party"
    model.frame(formula, ...)
}
