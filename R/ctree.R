
### calculate quad p-values
.MPinv <- function (X, tol = sqrt(.Machine$double.eps)) 
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    Xplus <- if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2:1])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
    list(Xplus = Xplus, rank = sum(Positive))
}

.pX2 <- function(lin, exp, cov, pval = TRUE) {
    if (length(lin) == 1) {
        if (cov < .Machine$double.eps) return(c(-Inf, -Inf))
        X2 <- ((lin - exp)^2) / cov
        df <- 1
    } else {
        tmp <- matrix(lin - exp, ncol = 1)
        Xplus <- .MPinv(matrix(cov, ncol = length(lin)))
        X2 <- crossprod(tmp, Xplus$Xplus) %*% tmp
        df <- Xplus$rank
    }
    
    X2 <- pmax(0, X2) ## Z: X2 may be slightly negative...
    
    if (pval)
        return(c(log(X2), pchisq(X2, df = df, lower.tail = TRUE, 
                                 log.p = TRUE)))
    return(c(log(X2), NA))
}

### calculate max-T p-value
.pmaxT <- function(lin, exp, cov, pval = TRUE) {

    if (length(lin) == 1) {
        if (cov < .Machine$double.eps) return(c(-Inf, -Inf))
        maxT <- abs(lin - exp) / sqrt(cov)
        v <- 1
        V <- matrix(1)
    } else {
        v <- diag(V <- matrix(cov, ncol = length(lin)))
        lin <- as.vector(lin)[v > 0]
        exp <- as.vector(exp)[v > 0]
        V <- V[v > 0, v > 0, drop = FALSE]
        v <- v[v > 0]
        if (length(v) == 0) return(c(-Inf, -Inf))
        maxT <- as.vector(max(abs(lin - exp) / sqrt(v)))
        if (is.na(maxT)) return(c(-Inf, -Inf))
    }
    if (pval) 
        return(c(log(maxT), log(pmvnorm(lower = rep(-maxT, length(v)),
                                upper = rep(maxT, length(v)),
                                sigma = cov2cor(V)))))
    return(c(log(maxT), NA))
}

### surrogate splits
.csurr <- function(split, data, inp, weights, ctrl) {

    ### <FIXME> surrogate splits for multiway splits </FIXME>?
    stopifnot(length(unique(split)) == 2)
    response <- as.factor(split)
    response <- model.matrix(~ response - 1)
    if (ncol(response) == 2) response <- response[, -1, drop = FALSE]
    storage.mode(response) <- "double"

    lin <- .Call("R_LinstatExpCov", data, inp, response, weights)
    p <- sapply(lin[inp], function(x) do.call(".pX2", x[-1]))
    colnames(p) <- colnames(data)[inp]
    rownames(p) <- c("teststat", "pval")
    crit <- p["pval",]

    ret <- vector(mode = "list", length = min(sum(inp), ctrl$maxsurrogate))

    for (i in 1:length(ret)) {
        isel <- which.max(crit)
        isel <- which(inp)[isel]
        x <- data[[isel]]
        sp <- .Call("R_split", x, response, weights, as.integer(0))
        if (any(is.na(sp))) next
        if (length(sp) == 1) {
            ret[[i]] <- partysplit(as.integer(isel), breaks = sp, index = 1:2)
        } else {
            ret[[i]] <- partysplit(as.integer(isel), index = sp)
        }
        tmp <- kidids_split(ret[[i]], data, obs = weights > 0)
        tmps <- split[weights > 0]
        tab <- table(tmp, tmps)
        if (tab[1, 1] < tab[1, 2]) {
            indx <- ret[[i]]$index
            ret[[i]]$index[indx == 1] <- 2
            ret[[i]]$index[indx == 2] <- 1
        }
        crit[which.max(crit)] <- -Inf
    }
    ret <- ret[!sapply(ret, is.null)]
    if (length(ret) == 0) ret <- NULL
    return(ret)
}

### set up new node for conditional inference tree
.cnode <- function(id = 1, data, response, inputs, weights, ctrl, cenv = NULL) {

    if (is.null(cenv)) {
        cenv <- new.env()
        depth <- 0
    } else {
        depth <- get("depth", envir = cenv)
        if (depth >= ctrl$maxdepth)
            return(partynode(as.integer(id)))
    }
    weights <- as.integer(weights)
    if (sum(weights) < ctrl$minsplit) return(partynode(as.integer(id)))
    if (id > 1 && ctrl$stump) return(partynode(as.integer(id)))

    inp <- inputs
    if (ctrl$mtry < Inf) {
        mtry <- min(sum(inp), ctrl$mtry)
        s <- sample(which(inp), mtry)
        inp <- logical(length(inp))
        inp[s] <- TRUE
    } 
    lin <- .Call("R_LinstatExpCov", data, inp, response, weights)
    p <- sapply(lin[inp], function(x) do.call(ctrl$cfun, x[-1]))
    crit <- p[1,,drop = FALSE]
    p <- p[-1,,drop = FALSE]
    colnames(p) <- colnames(data)[inp]

    mb <- ctrl$minbucket
    mp <- ctrl$minprob
    storage.mode(mb) <- "integer"

    count <- 1
    thissplit <- NULL
    while(count <= ctrl$splittry) {
        if (any(crit > ctrl$mincriterion)) {
            isel <- which.max(crit)
            isel <- which(inp)[isel]
        } else {
            return(partynode(as.integer(id), info = exp(p)))
        }
        x <- data[[isel]]
        swp <- ceiling(sum(weights) * mp)
        if (mb < swp) mb <- as.integer(swp)

        if ((ctrl$multiway && ctrl$maxsurrogate == 0) && is.factor(x)) {
            if (all(table(x[rep(1:length(x), weights)]) > mb)) {
                thissplit <- partysplit(as.integer(isel), index = 1:levels(x))
                break()
            }
        } else {
            sp <- .Call("R_split", x, response, weights, mb)
            if (!any(is.na(sp))) {
                if (length(sp) == 1) {
                    thissplit <- partysplit(as.integer(isel), breaks = sp)
                } else {
                    ### deal with empty levels -> NA in sp
                    if (is.factor(x)) 
                        sp[table(rep(x, weights)) == 0] <- NA
                    thissplit <- partysplit(as.integer(isel), index = sp)
                }
                break()
            }
        }
        crit[which.max(crit)] <- -Inf
        count <- count + 1
    }
    if (is.null(thissplit))
        return(partynode(as.integer(id), info = exp(p)))

    ret <- partynode(as.integer(id))
    ret$split <- thissplit
    ret$info <- exp(p)
    thissurr <- NULL
    kidids <- kidids_node(ret, data)

    if (ctrl$maxsurrogate > 0) {
        inp <- inputs
        inp[isel] <- FALSE
        w <- weights
        xna <- is.na(x)
        w[xna] <- 0L
        ret$surrogates <- .csurr(kidids, data, inp, w, ctrl)
        kidids[xna] <- kidids_node(ret, data, obs = xna)
    }

    kids <- vector(mode = "list", length = max(kidids)) ## Z: was 1:max(kidids)
    nextid <- id + 1
    for (k in 1:max(kidids)) {
        w <- weights
        w[kidids != k] <- 0
        assign("depth", depth + 1, envir = cenv)
        kids[[k]] <- .cnode(nextid, data, response, inputs, w, ctrl, cenv)
        nextid <- max(nodeids(kids[[k]])) + 1
    }
    ret$kids <- kids

    return(ret)
}

ctree_control <- function(teststat = c("quad", "max"),
    testtype = c("Bonferroni", "Univariate", "Teststatistic"),
    mincriterion = 0.95, minsplit = 20L, minbucket = 7L, minprob = 0.01,
    stump = FALSE, maxsurrogate = 0L, mtry = Inf, maxdepth = Inf, 
    multiway = FALSE, splittry = 2L) {

    teststat <- match.arg(teststat)
    if (teststat == "max") stopifnot(require("mvtnorm"))
    testtype <- match.arg(testtype)
    list(teststat = teststat,
         testtype = testtype, mincriterion = log(mincriterion),
         minsplit = minsplit, minbucket = minbucket, 
         minprob = minprob, stump = stump, mtry = mtry, 
         maxdepth = maxdepth, multiway = multiway, splittry = splittry,
         maxsurrogate = maxsurrogate)
}

ctree <- function(formula, data, weights, subset, na.action = na.pass, 
                  control = ctree_control(...), ...) {

    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mf), 0)
    mf <- mf[c(1, m)]
    
    ### only necessary for extended model formulae 
    ### e.g. multivariate responses
    if (require("Formula")) {
        formula <- Formula(formula)
    } else {
        if (length(formula[[2]]) > 1)
            stop("Package ", sQuote("Formula"),
                 " not available for handling extended model formula ",
                 sQuote("formula"))
    }
    mf$formula <- formula
    mf$drop.unused.levels <- FALSE
    mf$na.action <- na.action
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    response <- names(mf)[1]
    if (inherits(formula, "Formula"))
        response <- names(model.part(formula, mf, lhs = 1))
    weights <- model.weights(mf)
    dat <- mf[, colnames(mf) != "(weights)"]
    ret <- .ctree_fit(dat, response, weights = weights, ctrl = control)
    ret$terms <- terms(formula, data = mf)
    ### need to adjust print and plot methods
    ### for multivariate responses
    if (length(response) > 1) class(ret) <- "party"
    return(ret)
}

.ctree_fit <- function(data, response, weights = NULL, 
                      ctrl = ctree_control()) {

    inputs <- !(colnames(data) %in% response)

    infl <- .y2infl(data, response)

    if (is.null(weights))
        weights <- rep(1, nrow(data))
    storage.mode(weights) <- "integer"

    ctrl$cfun <- function(...) {
        if (ctrl$teststat == "quad")
            p <- .pX2(..., pval = (ctrl$testtype != "Teststatistic"))
        if (ctrl$teststat == "max")
            p <- .pmaxT(..., pval = (ctrl$testtype != "Teststatistic"))
        names(p) <- c("teststat", "pval")

        if (ctrl$testtype == "Bonferroni")
            p["pval"] <- p["pval"] * min(sum(inputs), ctrl$mtry)
        crit <-  p["teststat"]
        if (ctrl$testtype != "Teststatistic")
        crit <- p["pval"]
        c(crit, p)
    }
    tree <- .cnode(1L, data, infl, inputs, weights, ctrl)
    fitted <- data.frame("(fitted)" = fitted_node(tree, data), 
                         "(weights)" = weights,
                         check.names = FALSE)
    fitted[[3]] <- data[, response, drop = length(response) == 1]
    names(fitted)[3] <- "(response)"
    ret <- party(tree, data = data, fitted = fitted)
    class(ret) <- c("constparty", class(ret))
    return(ret)
}

.logrank_trafo <- function(x, ties.method = c("logrank", "HL")) {
    ties.method <- match.arg(ties.method)
    time <- x[,1]
    event <- x[,2]
    n <- length(time)
    ot <- order(time, event)
    rt <- rank(time, ties.method = "max")
    mt <- rank(time, ties.method = "min") - 1
    fact <- switch(ties.method, "logrank" = event / (n - mt),
                                "HL" = event/(n - rt + 1)
                  )   
    event - cumsum(fact[ot])[rt]
}

### convert response y to influence function h(y)
.y2infl <- function(data, response) {

    if (length(response) == 1) {
        response <- data[[response]]
        rtype <- class(response)[1]
        if (rtype == "integer") rtype <- "numeric"

        infl <- switch(rtype,
            "factor" = { 
                X <- model.matrix(~ response - 1)
                if (nlevels(response) > 2) return(X)
                return(X[,-1, drop = FALSE])
            },
            "ordered" = (1:nlevels(response))[as.integer(response)],
            "numeric" = response,
            "Surv" = .logrank_trafo(response)
        )
    } else {
        ### multivariate response
        infl <- lapply(response, .y2infl, data = data)
        tmp <- do.call("cbind", infl)
        attr(tmp, "assign") <- rep(1:length(infl), sapply(infl, NCOL))
        infl <- tmp
    }
    storage.mode(infl) <- "double"
    return(infl)
}
