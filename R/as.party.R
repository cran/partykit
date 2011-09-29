as.party <- function(obj, ...)
    UseMethod("as.party")

as.party.rpart <- function(obj, ...) {

    ff <- obj$frame
    n  <- nrow(ff)
    if (n==1) return(partynode(as.integer(1)))  # special case of no splits

    is.leaf <- (ff$var == "<leaf>")
    vnames <- ff$var[!is.leaf]  #the variable names for the primary splits

    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    splitindex <- list()
    splitindex$primary <- numeric(n)
    splitindex$primary[!is.leaf] <- index[c(!is.leaf, FALSE)]
    splitindex$surrogate <- lapply(1:n, function(i) {
        prim <- splitindex$primary[i]
        if (prim < 1 || ff[i, "nsurrogate"] == 0) return(NULL)
        else return(prim + ff[i, "ncompete"] + 1:ff[i, "nsurrogate"])
    })
    
    mf <- model.frame(obj)

    rpart_fitted <- function() {
        y <- model.response(mf)
        weights <- model.weights(mf)
        ret <- data.frame("(fitted)" = obj$where, "(response)" = y, check.names = FALSE)
        if (!is.null(weights)) ret[["(weights)"]] <- weights
        ret
    }
    fitted <- rpart_fitted()

    rpart_kids <- function(i) {
        if (is.leaf[i]) return(NULL)
        else return(c(i + 1, 
            which((cumsum(!is.leaf[-(1:i)]) + 1) == cumsum(is.leaf[-(1:i)]))[1] + 1 + i))
    }

    rpart_onesplit <- function(j) {
        if (j < 1) return(NULL)
        ### numeric
        if (abs(obj$split[j, "ncat"]) == 1) {
            ret <- partysplit(varid = which(rownames(obj$split)[j] == names(mf)),
                      breaks = as.double(obj$split[j, "index"]),
                      right = FALSE,
                      index = if(obj$split[j, "ncat"] > 0) 2:1)
        } else {
            index <- obj$csplit[obj$split[j, "index"],]
            ### csplit has columns 1:max(nlevels) for all factors
            index <- index[1:nlevels(mf[, rownames(obj$split)[j]])]
            index[index == 2] <- NA ### level not present in split
            index[index == 3] <- 2  ### 1..left, 3..right
            ret <- partysplit(varid = which(rownames(obj$split)[j] == names(mf)),
                      index = as.integer(index))
        }
        ret
    }
                      
    rpart_split <- function(i)
        rpart_onesplit(splitindex$primary[i])
    
    rpart_surrogates <- function(i)
        lapply(splitindex$surrogate[[i]], rpart_onesplit)

    rpart_node <- function(i) {
        if (is.null(rpart_kids(i))) return(partynode(as.integer(i)))
        nd <- partynode(as.integer(i), split = rpart_split(i),
	           kids = lapply(rpart_kids(i), rpart_node),
	           surrogates = rpart_surrogates(i))

        ### determine majority for (non-random) splitting
        left <- nodeids(kids_node(nd)[[1]], terminal = TRUE)
        right <- nodeids(kids_node(nd)[[2]], terminal = TRUE)
        nd$split$prob <- c(0, 0)
        nl <- sum(fitted[["(fitted)"]] %in% left)
        nr <- sum(fitted[["(fitted)"]] %in% right)
        nd$split$prob <- if (nl > nr) c(1, 0) else c(0, 1)
        nd$split$prob <- as.double(nd$split$prob)
        return(nd)
    }

    node <- rpart_node(1)

    rval <- party(node = node, data = mf[0,], fitted = fitted,
      terms = obj$terms, info = list(method = "rpart"))
    class(rval) <- c("constparty", class(rval))
    return(rval)
}

model.frame.rpart <- function(formula, ...) {
  mf <- formula$call
  mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L))]
  if (is.null(mf$na.action)) mf$na.action <- na.rpart
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  env <- if (is.null(environment(formula$terms))) environment(formula$terms) 
             else parent.frame()
  mf <- eval(mf, env)
  return(mf)
}

as.party.J48 <- function(obj, ...) {

  ## construct metadata
  mf <- model.frame(obj)
  mf_class <- sapply(mf, class)
  mf_levels <- lapply(mf, levels)

  x <- rJava::.jcall(obj$classifier, "S", "graph")
  x <- RWeka:::parse_Weka_digraph(x, plainleaf = TRUE)
  nodes <- x$nodes
  edges <- x$edges
  is.leaf <- x$nodes[, "splitvar"] == ""

  j48_kids <- function(i) {
    if (is.leaf[i]) return(NULL)
      else return(which(nodes[,"name"] %in% edges[nodes[i,"name"] == edges[,"from"], "to"]))
  }

  j48_split <- function(i) {
    if(is.leaf[i]) return(NULL)
    
    var_id <- which(nodes[i, "splitvar"] == names(mf))
    ##
    edges <- edges[nodes[i,"name"] == edges[,"from"], "label"]
    split <- Map(c, sub("^([[:punct:]]+).*$", "\\1", edges), sub("^([[:punct:]]+) *", "", edges))
    ## ## for J48 the following suffices
    ## split <- strsplit(edges[nodes[i,"name"] == edges[,"from"], "label"], " ")

    if(mf_class[var_id] %in% c("ordered", "factor")) {
      stopifnot(all(sapply(split, head, 1) == "="))
      stopifnot(all(sapply(split, tail, 1) %in% mf_levels[[var_id]]))
      
      split <- partysplit(varid = as.integer(var_id),
        index = match(mf_levels[[var_id]], sapply(split, tail, 1)))
    } else {
      breaks <- unique(as.numeric(sapply(split, tail, 1)))
      breaks <- if(mf_class[var_id] == "integer") as.integer(breaks) else as.double(breaks) ## FIXME: check
      
      stopifnot(length(breaks) == 1 && !is.na(breaks))
      stopifnot(all(sapply(split, head, 1) %in% c("<=", ">")))
      
      split <- partysplit(varid = as.integer(var_id),
        breaks = breaks, right = TRUE,
	index = if(split[[1]][1] == ">") 2:1)
    }
    return(split)
  }

  j48_node <- function(i) {
    if(is.null(j48_kids(i))) return(partynode(as.integer(i)))
    partynode(as.integer(i), split = j48_split(i), kids = lapply(j48_kids(i), j48_node))
  }

  node <- j48_node(1)

  j48 <- party(node = node,
               data = mf[0,],
               fitted = data.frame("(fitted)" = fitted_node(node, mf),
	                           "(response)" = model.response(mf),
				   check.names = FALSE),
               terms = obj$terms,
	       info = list(method = "J4.8"))

  class(j48) <- c("constparty", class(j48))
  return(j48)
}
