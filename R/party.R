## FIXME: data in party
##   - currently assumed to be a data.frame
##   - potentially empty
##   - the following are all assumed to work:
##     dim(data), names(data)
##     sapply(data, class), lapply(data, levels)
##   - potentially these need to be modified if data/terms
##     should be able to deal with data bases

party <- function(node, data, fitted = NULL, terms = NULL, names = NULL, info = NULL) {

    stopifnot(inherits(node, "partynode"))
    stopifnot(inherits(data, "data.frame"))
    ### make sure all split variables are there 
    ids <- nodeids(node)[!nodeids(node) %in% nodeids(node, terminal = TRUE)]
    varids <- unique(unlist(nodeapply(node, ids = ids, FUN = function(x) 
        varid_split(split_node(x)))))
    stopifnot(varids %in% 1:ncol(data))
    
    if(!is.null(fitted)) {
        stopifnot(inherits(fitted, "data.frame"))
        stopifnot("(fitted)" == names(fitted)[1])
        stopifnot(nrow(data) == 0 | nrow(data) == nrow(fitted))

        nt <- nodeids(node, terminal = TRUE)
        stopifnot(all(fitted[["(fitted)"]] %in% nt))

        node <- as.partynode(node, from = 1L)
        nt2 <- nodeids(node, terminal = TRUE)
        fitted[["(fitted)"]] <- nt2[match(fitted[["(fitted)"]], nt)]
    } else {
        node <- as.partynode(node, from = 1L)
    }
    
    party <- list(node = node, data = data, fitted = fitted, 
                  terms = NULL, names = NULL, info = info)
    class(party) <- "party"

    if(!is.null(terms)) {
        stopifnot(inherits(terms, "terms"))
	party$terms <- terms
    }

    if (!is.null(names)) {
        n <- length(nodeids(party, terminal = FALSE))
        if (length(names) != n)
            stop("invalid", " ", sQuote("names"), " ", "argument")
        party$names <- names
    }

    party
}

length.party <- function(x)
    length(nodeids(x))

names.party <- function(x)
    .names_party(x)

"names<-.party" <- function(x, value) {
     n <- length(nodeids(x, terminal = FALSE))
     if (!is.null(value) && length(value) != n)
         stop("invalid", " ", sQuote("names"), " ", "argument")
     x$names <- value
     x
}

.names_party <- function(party) {
    names <- party$names
    if (is.null(names))
        names <- as.character(nodeids(party, terminal = FALSE))
    names
}

node_party <- function(party) {
    stopifnot(inherits(party, "party"))
    party$node
}

is.constparty <- function(party) {
    stopifnot(inherits(party, "party"))
    if (!is.null(party$fitted)) 
        return(all(c("(fitted)", "(response)") %in% names(party$fitted)))
    return(FALSE)
}

as.constparty <- function(obj, ...) {
    if (is.constparty(obj)) {
        ret <- obj
        class(ret) <- c("constparty", class(obj))
        return(ret)
    }
    stop("cannot coerce object of class", " ", sQuote(class(obj)), 
          " ", "to", " ", sQuote("constparty"))
}

"[.party" <- "[[.party" <- function(x, i, ...) {
    if (is.character(i) && !is.null(names(x)))
        i <- which(names(x) %in% i)
    stopifnot(length(i) == 1 & is.numeric(i))
    stopifnot(i <= length(x) & i >= 1)
    i <- as.integer(i)
    dat <- data_party(x, i)
    if (!is.null(x$fitted)) {
        findx <- which("(fitted)" == names(dat))[1]
        fit <- dat[,findx:ncol(dat), drop = FALSE]
        dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
        if (ncol(dat) == 0)
            dat <- x$data
    } else {
        fit <- NULL
        dat <- x$data
    }
    nam <- names(x)[nodeids(x, from = i, terminal = FALSE)]

    recFun <- function(node) {
        if (id_node(node) == i) return(node)
        kid <- sapply(kids_node(node), id_node)
        return(recFun(node[[max(which(kid <= i))]]))
    }
    node <- recFun(node_party(x))

    ret <- party(node = node, data = dat, fitted = fit, 
                 terms = x$terms, names = nam, info = x$info)
    class(ret) <- class(x)
    ret
}

nodeids <- function(obj, ...)
    UseMethod("nodeids")

nodeids.partynode <- function(obj, from = NULL, terminal = FALSE, ...) {

    if(is.null(from)) from <- id_node(obj)

    id <- function(node, record = TRUE, terminal = FALSE) {
      if(!record) return(NULL)
      if(!terminal)
          return(id_node(node))
      else
          if(is.terminal(node)) return(id_node(node)) else return(NULL)
    }

    rid <- function(node, record = TRUE, terminal = FALSE) {  
        myid <- id(node, record = record, terminal = terminal)
        if(is.terminal(node)) return(myid)
        kids <- kids_node(node)    
        kids_record <- if(record)  
            rep(TRUE, length(kids))
        else
            sapply(kids, id_node) == from
        return(c(myid,
            unlist(lapply(1:length(kids), function(i)
                rid(kids[[i]], record = kids_record[i], terminal = terminal)))
        ))
    }

    return(rid(obj, from == id_node(obj), terminal))
}

nodeids.party <- function(obj, from = NULL, terminal = FALSE, ...)
    nodeids(node_party(obj), from = from, terminal = terminal, ...)

nodeapply <- function(obj, ids = 1, FUN = NULL, ...)
    UseMethod("nodeapply")

nodeapply.party <- function(obj, ids = 1, FUN = NULL, by_node = TRUE, ...) {

    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)

    if(is.null(FUN)) FUN <- function(x, ...) x

    if (length(ids) == 0)
        return(NULL)

    if (by_node) {
        rval <- nodeapply(node_party(obj), ids = ids, FUN = FUN, ...)
    } else {
        rval <- lapply(ids, function(i) FUN(obj[[i]], ...))
    }

    names(rval) <- names(obj)[ids]
    return(rval)
}

nodeapply.partynode <- function(obj, ids = 1, FUN = NULL, ...) {

    stopifnot(isTRUE(all.equal(ids, round(ids))))
    ids <- as.integer(ids)

    if(is.null(FUN)) FUN <- function(x, ...) x

    if (length(ids) == 0)
        return(NULL)

    rval <- vector(mode = "list", length = length(ids))
    rval_id <- rep(0, length(ids))
    i <- 1
	
    recFUN <- function(node, ...) {
        if(id_node(node) %in% ids) {
            rval_id[i] <<- id_node(node)
            rval[[i]] <<- FUN(node, ...)
            i <<- i + 1
        }
        kids <- kids_node(node)
        if(length(kids) > 0) {
            for(j in 1:length(kids)) recFUN(kids[[j]])
        }
        invisible(TRUE)
    }
    foo <- recFUN(obj)
    rval <- rval[match(rval_id, ids)]
    return(rval)
}

predict.party <- function(object, newdata = NULL, ...)
{
    ### compute fitted node ids first
    fitted <- if(is.null(newdata)) object$fitted[["(fitted)"]] else {

        terminal <- nodeids(object, terminal = TRUE)
        inner <- 1:max(terminal)
        inner <- inner[-terminal]

        primary_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
            varid_split(split_node(node))
        })
        surrogate_vars <- nodeapply(object, ids = inner, by_node = TRUE, FUN = function(node) {
            surr <- surrogates_node(node)
            if(is.null(surr)) return(NULL) else return(sapply(surr, varid_split))
        })
        vnames <- names(object$data)

        ## ## FIXME: the is.na() call takes loooong on large data sets
        ## unames <- if(any(sapply(newdata, is.na))) 
        ##     vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
        ## else 
        ##     vnames[unique(unlist(primary_vars))]
	unames <- vnames[unique(unlist(c(primary_vars, surrogate_vars)))]
	
        vclass <- structure(lapply(object$data, class), .Names = vnames)
        ndnames <- names(newdata)
        ndclass <- structure(lapply(newdata, class), .Names = ndnames)
        if(all(unames %in% ndnames) &&
           all(unlist(lapply(unames, function(x) vclass[[x]] == ndclass[[x]])))) {
            vmatch <- match(vnames, ndnames)
            fitted_node(node_party(object), newdata, vmatch)
        } else {
            if (!is.null(object$terms)) {
                mf <- model.frame(delete.response(object$terms), newdata)
                fitted_node(node_party(object), mf, match(vnames, names(mf)))
            } else
                stop("") ## FIXME: write error message
        }
    }
    ### compute predictions
    predict_party(object, fitted, newdata, ...)
}

predict_party <- function(party, id, newdata = NULL, ...)
    UseMethod("predict_party")

### do nothing expect returning the fitted ids
predict_party.default <- function(party, id, newdata = NULL, ...) {

    if (length(list(...)) > 1) 
        warning("argument(s)", " ", sQuote(names(list(...))), " ", "have been ignored")

    ## get observation names: either node names or
    ## observation names from newdata
    nam <- if(is.null(newdata)) names(party)[id] else rownames(newdata)
    if(length(nam) != length(id)) nam <- NULL

    ## special case: fitted ids
    return(structure(id, .Names = nam))
}

predict_party.constparty <- function(party, id, newdata = NULL,
    type = c("response", "prob", "node"), FUN = NULL, simplify = TRUE, ...)
{
    ## extract fitted information
    response <- party$fitted[["(response)"]]
    weights <- party$fitted[["(weights)"]]
    fitted <- party$fitted[["(fitted)"]]
    if (is.null(weights)) weights <- rep(1, NROW(response))

    ## get observation names: either node names or
    ## observation names from newdata
    nam <- if(is.null(newdata)) names(party)[id] else rownames(newdata)
    if(length(nam) != length(id)) nam <- NULL

    ## match type
    type <- match.arg(type)

    ## special case: fitted ids
    if(type == "node")
      return(structure(id, .Names = nam))

    ### multivariate response
    if (is.data.frame(response)) {
        ret <- lapply(response, function(r) {
            ret <- .predict_party_constparty(node_party(party), fitted = fitted, 
                response = r, weights, id = id, type = type, FUN = FUN, ...)
            if (simplify) .simplify_pred(ret, id, nam) else ret
        })
        if (all(sapply(ret, is.atomic)))
            ret <- as.data.frame(ret)
        names(ret) <- colnames(response)
        return(ret)
    }

    ### univariate response
    ret <- .predict_party_constparty(node_party(party), fitted = fitted, response = response, 
        weights = weights, id = id, type = type, FUN = FUN, ...)
    if (simplify) .simplify_pred(ret, id, nam) else ret[as.character(id)]
}

### functions for node prediction based on fitted / response
.pred_Surv <- function(y, w)
    survival:::survfit(y ~ 1, weights = w, subset = w > 0)

.pred_Surv_response <- function(y, w)
    .median_survival_time(.pred_Surv(y, w))
                    
.pred_factor <- function(y, w) {
    lev <- levels(y)
    sumw <- tapply(w, y, sum)
    sumw[is.na(sumw)] <- 0
    prob <- sumw / sum(w)
    names(prob) <- lev
    return(prob)
}

.pred_factor_response <- function(y, w) {
    prob <- .pred_factor(y, w)
    return(factor(which.max(prob), levels = 1:nlevels(y),
                  labels = levels(y), 
                  ordered = is.ordered(y)))
    return(prob) 
}
                    
.pred_numeric <- function(y, w) weighted.mean(y, w, na.rm = TRUE)

### workhorse: compute predictions based on fitted / response data
.predict_party_constparty <- function(node, fitted, response, weights,
    id = id, type = c("response", "prob"), FUN = NULL, ...) {

    if (is.null(FUN)) {

        rtype <- class(response)[1]
        if (rtype == "ordered") rtype <- "factor"    
        if (rtype == "integer") rtype <- "numeric"

        FUN <- switch(rtype,
            "Surv" = if (type == "response") .pred_Surv_response else .pred_Surv,
            "factor" = if (type == "response") .pred_factor_response else .pred_factor,
            "numeric" = {
                if (type == "prob")
                    stop(sQuote("type = \"prob\""), " ", "is not available")
                .pred_numeric
           })
    }
      
    ## empirical distribution in each leaf
    if (all(id %in% fitted)) {
        tab <- tapply(1:NROW(response), fitted, 
                      function(i) FUN(response[i], weights[i]), simplify = FALSE)
    } else {
        ### id may also refer to inner nodes
        tab <- as.array(lapply(sort(unique(id)), function(i) {
            index <- fitted %in% nodeids(node, i, terminal = TRUE)
            FUN(response[index], weights[index])
        }))
        names(tab) <- as.character(sort(unique(id)))
    }
    tn <- names(tab)
    dim(tab) <- NULL
    names(tab) <- tn

    tab
}


### simplify structure of predictions
.simplify_pred <- function(tab, id, nam) {

    if (all(sapply(tab, length) == 1) & all(sapply(tab, is.atomic))) {
        ret <- do.call("c", tab)
        names(ret) <- names(tab)
        ret <- if (is.factor(tab[[1]]))
            factor(ret[as.character(id)], levels = 1:length(levels(tab[[1]])),
		   labels = levels(tab[[1]]), ordered = is.ordered(tab[[1]]))
        else 
            ret[as.character(id)]
        names(ret) <- nam
    } else if (length(unique(sapply(tab, length))) == 1 & 
               all(sapply(tab, is.numeric))) {
        ret <- matrix(unlist(tab), nrow = length(tab), byrow = TRUE)
        colnames(ret) <- names(tab[[1]])
        rownames(ret) <- names(tab)
        ret <- ret[as.character(id),, drop = FALSE]
        rownames(ret) <- nam
    } else {
        ret <- tab[as.character(id)]
        names(ret) <- nam
    }
    ret
}

data_party <- function(party, id = 1L)
    UseMethod("data_party")

data_party.default <- function(party, id = 1L) {
    
    extract <- function(id) {
        if(is.null(party$fitted))
            if(nrow(party$data) == 0) return(NULL)
        else
            stop("cannot subset data without fitted ids")

        ### which terminal nodes follow node number id?
        nt <- nodeids(party, id, terminal = TRUE)
        wi <- party$fitted[["(fitted)"]] %in% nt

        ret <- if (nrow(party$data) == 0)
            subset(party$fitted, wi)
        else
            subset(cbind(party$data, party$fitted), wi)
        ret
    }
    if (length(id) > 1)
        return(lapply(id, extract))
    else 
        return(extract(id))
}

width.party <- function(x, ...) {
  width(node_party(x), ...)
}

depth.party <- function(x, root = FALSE, ...) {
  depth(node_party(x), root = root, ...)
}

.list.rules.party <- function(x, i = NULL, ...) {
    if (is.null(i)) i <- nodeids(x, terminal = TRUE)
    if (length(i) > 1) {
        ret <- sapply(i, .list.rules.party, x = x)
        names(ret) <- if (is.character(i)) i else names(x)[i]
        return(ret)
    }
    if (is.character(i) && !is.null(names(x)))
        i <- which(names(x) %in% i)
    stopifnot(length(i) == 1 & is.numeric(i))
    stopifnot(i <= length(x) & i >= 1)
    i <- as.integer(i)
    dat <- data_party(x, i)  
    if (!is.null(x$fitted)) {
        findx <- which("(fitted)" == names(dat))[1]  
        fit <- dat[,findx:ncol(dat), drop = FALSE]   
        dat <- dat[,-(findx:ncol(dat)), drop = FALSE]
        if (ncol(dat) == 0)
            dat <- x$data
    } else {
        fit <- NULL  
        dat <- x$data
    }

    rule <- c()

    recFun <- function(node) {
        if (id_node(node) == i) return(NULL)   
        kid <- sapply(kids_node(node), id_node)
        whichkid <- max(which(kid <= i))
        split <- split_node(node)
        ivar <- varid_split(split)
        svar <- names(dat)[ivar]
        index <- index_split(split)
        if (is.factor(dat[, svar])) {
            slevels <- levels(dat[, svar])[index == whichkid]
            srule <- paste(svar, " %in% c(\"", 
                paste(slevels, collapse = "\", \"", sep = ""), "\")",
                sep = "")
        } else {
            if (is.null(index)) index <- 1:length(kid)
            breaks <- cbind(c(-Inf, breaks_split(split)), 
                            c(breaks_split(split), Inf))
            sbreak <- breaks[index == whichkid,]
            right <- right_split(split)
            srule <- c()
            if (is.finite(sbreak[1]))
                srule <- c(srule, 
                    paste(svar, ifelse(right, ">", ">="), sbreak[1]))
            if (is.finite(sbreak[2]))
                srule <- c(srule, 
                    paste(svar, ifelse(right, "<=", ">"), sbreak[2]))
            srule <- paste(srule, collapse = " & ")
        }
        rule <<- c(rule, srule)
        return(recFun(node[[whichkid]]))
    }
    node <- recFun(node_party(x))
    paste(rule, collapse = " & ")
}
