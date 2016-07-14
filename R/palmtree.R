utils::globalVariables(c(".tree", ".lm", ".weights"))

palmtree <- function(formula, data, weights = NULL, family = NULL,
  lmstart = NULL, abstol = 0.001, maxit = 100, 
  dfsplit = TRUE, verbose = FALSE, plot = FALSE, ...)
{
  ## remember call
  cl <- match.call()

  ## process family
  if(!is.null(family)) {
    if(!inherits(family, "family")) {
      if(is.character(family)) family <- get(family)
      if(is.function(family)) family <- family()
    }
  }
  
  ## formula processing (full, tree, regression)
  ff <- Formula::as.Formula(formula)
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 3L))
  rf <- formula(ff, lhs = 1L, rhs = 1L)
  ## intercept in lmtree?
  intercept <- attr(terms(rf), "intercept") > 0L
  if(!intercept) rf <- update(rf, . ~ . + 1)
  ## without tree
  rf0 <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
    lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
  ## with tree
  rf <- if(intercept) update(rf, . ~ .tree / .) else update(rf, . ~ .tree:.)
  rf <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
    lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)

  ## weights
  data$.weights <- if(is.null(weights)) rep(1, nrow(data)) else weights

  ## initialization
  iteration <- 0L
  if (is.null(lmstart)) {
    lmstart <- if(is.null(family)) {
      lm(rf0, data = data, weights = .weights)
    } else {
      glm(rf0, data = data, weights = .weights, family = family)
    }
    lmstart <- predict(lmstart)
  }
  data$.lm <- lmstart  

  continue <- TRUE
  oldloglik <- -Inf

  ## iterate between palm and lmtree estimation
  while (continue) {
    iteration <- iteration + 1L

    ## lmtree
    tree <- if(is.null(family)) {
      lmtree(tf, data = data, offset = .lm, weights = .weights, dfsplit = FALSE, ...)
    } else {
      glmtree(tf, data = data, offset = .lm, weights = .weights, dfsplit = FALSE, family = family, ...)    
    }
    if(plot) plot(tree)
    data$.tree <- factor(predict(tree, newdata = data, type = "node"))

    ## estimate full lm model but force all coefficients from the
    ## .tree (and the overall intercept) to zero for the prediction
    palm <- if(length(levels(data$.tree)) > 1L) {
      if(is.null(family)) {
        lm(rf, data = data, weights = .weights)
      } else {
        glm(rf, data = data, weights = .weights, family = family)      
      }
    } else {
      if(is.null(family)) {
        lm(rf0, data = data, weights = .weights)
      } else {
        glm(rf0, data = data, weights = .weights, family = family)      
      }
    }
    b <- palm$coefficients
    if(length(levels(data$.tree)) > 1L) {
      palm$coefficients[substr(names(b), 1L, 5L) %in% c(if(intercept) "(Inte" else NULL, ".tree")] <- 0
    } else {
      palm$coefficients[names(coef(tree))] <- 0
    }
    data$.lm <- predict(palm, newdata = data)
    palm$coefficients <- b

    ## iteration information
    newloglik <- logLik(palm)
    continue <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if(verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    palm = palm,
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    lmstart = lmstart, 
    abstol = abstol,
    intercept = intercept,
    family = family,
    mob.control = list(...)
  )
  class(result) <- "palmtree"
  return(result)
}

coef.palmtree <- function(object, model = c("all", "tree", "palm"), ...)
{
  model <- match.arg(model, c("all", "tree", "palm"))
  if(model == "all")  return(coef(object$palm, ...))
  if(model == "tree") return(coef(object$tree, ...))
  
  cf_all  <- coef(object$palm)
  cf_tree <- coef(object$tree, drop = FALSE)
  excl <- which(substr(names(cf_all), 1L, 5L) %in% c(if(object$intercept) "(Inte" else NULL, ".tree"))
  excl <- c(excl, which(names(cf_all) %in% colnames(cf_tree)))
  return(cf_all[-excl])
}

plot.palmtree <- function(x, ...) {
  plot(x$tree, ...)
}

logLik.palmtree <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$dfsplit
  dfsplit <- as.integer(dfsplit) * (length(object$tree) - length(nodeids(object$tree, terminal = TRUE)))
  structure(object$loglik, df = object$df + dfsplit, nobs = object$nobs, class = "logLik")
}

print.palmtree <- function(x, title = NULL, ...)
{
  if(is.null(title)) {
    title <- if(is.null(x$family)) {
      "Partially additive linear model tree"
    } else {
      sprintf("Partially additive generalized linear model tree (family: %s)", x$family$family)
    }
  }
  print(x$tree, title = title, ...)
  if(length(coef(x$palm)[-grep(".tree", names(coef(x$palm)))]) > 1L) {
    cat("\nLinear fixed effects (from palm model):\n")
    print(coef(x, model = "palm"))
  }
  invisible(x)
}

predict.palmtree <- function(object, newdata = NULL, type = "response", ...) { 
  if(is.null(newdata)) {
    newdata <- object$data
  }
  if(type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
    newdata$.tree <- factor(newdata$.tree,
    			    labels = levels(object$data$.tree))
    predict(object$palm, newdata = newdata, type = type, ...)
  }
}


dgp <- function(nobs = 1000, nreg = 5, creg = 0.4, ptreat = 0.5, sd = 1,
  coef = c(1, 0.25, 0.25, 0, 0, -0.25), eff = 1)
{
  d <- mvtnorm::rmvnorm(nobs, mean = rep(0, nreg), sigma = diag(1 - creg, nreg) + creg)
  colnames(d) <- paste0("x", 1:nreg)
  d <- as.data.frame(d)
  d$a <- rbinom(nobs, size = 1, prob = ptreat)
  d$err <- rnorm(nobs, mean = 0, sd = sd)

  gopt <- function(d) {
    as.numeric(d$x1 > -0.545) * as.numeric(d$x2 < 0.545)
  }
  d$y <- coef[1] + as.matrix(d[, paste0("x", 1:5)]) %*% coef[-1] - eff * (d$a - gopt(d))^2 + d$err
  d$a <- factor(d$a)
  return(d)
}
