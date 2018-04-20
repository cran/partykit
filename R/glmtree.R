## simple wrapper function to specify fitter and return class
glmtree <- function(formula, data, subset, na.action, weights, offset, cluster,
  family = gaussian, epsilon = 1e-8, maxit = 25, ...)
{
  ## use dots for setting up mob_control
  control <- mob_control(...)

  ## keep call
  cl <- match.call(expand.dots = TRUE)

  ## extend formula if necessary
  f <- Formula::Formula(formula)
  if(length(f)[2L] == 1L) {
    attr(f, "rhs") <- c(list(1), attr(f, "rhs"))
    formula[[3L]] <- formula(f)[[3L]]
  } else {
    f <- NULL
  }

  ## process family
  if(inherits(family, "family")) {
    fam <- TRUE
  } else {
    fam <- FALSE
    if(is.character(family)) family <- get(family)
    if(is.function(family)) family <- family()
  }

  ## distinguish whether glm should be fixed for case weights or not
  glmfit0 <- function(y, x, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
    estfun = FALSE, object = FALSE, caseweights = TRUE)
  {
    glmfit(y = y, x = x, start = start, weights = weights, offset = offset, cluster = cluster, ...,
      estfun = estfun, object = object, caseweights = control$caseweights)
  }

  ## call mob
  m <- match.call(expand.dots = FALSE)
  if(!is.null(f)) m$formula <- formula
  m$fit <- glmfit0
  m$control <- control
  m$epsilon <- epsilon
  m$maxit <- maxit
  if("..." %in% names(m)) m[["..."]] <- NULL
  if(!fam) m$family <- family
  m[[1L]] <- as.call(quote(partykit::mob))
  rval <- eval(m, parent.frame())
  
  ## extend class and keep original call
  rval$info$call <- cl
  rval$info$family <- family$family
  class(rval) <- c("glmtree", class(rval))
  return(rval)
}

## actual fitting function for mob()
glmfit <- function(y, x, start = NULL, weights = NULL, offset = NULL, cluster = NULL, ...,
  estfun = FALSE, object = FALSE, caseweights = TRUE)
{
  ## catch control arguments
  args <- list(...)
  ctrl <- list()
  for(n in c("epsilon", "maxit")) {
    if(n %in% names(args)) {
      ctrl[[n]] <- args[[n]]
      args[[n]] <- NULL
    }
  }
  args$control <- do.call("glm.control", ctrl)
  
  ## add intercept-only regressor matrix (if missing)
  ## NOTE: does not have terms/formula
  if(is.null(x)) x <- matrix(1, nrow = NROW(y), ncol = 1L,
    dimnames = list(NULL, "(Intercept)"))
  
  ## call glm fitting function
  args <- c(list(x = x, y = y, start = start, weights = weights, offset = offset), args)
  z <- do.call("glm.fit", args)

  ## degrees of freedom
  df <- z$rank
  if(z$family$family %in% c("gaussian", "Gamma", "inverse.gaussian")) df <- df + 1
  if(substr(z$family$family, 1L, 5L) != "quasi") objfun <- z$aic/2 - df else objfun <- z$deviance

  ## list structure
  rval <- list(
    coefficients = z$coefficients,
    objfun = objfun,
    estfun = NULL,
    object = NULL
  )

  ## add estimating functions (if desired)
  if(estfun) {
    wres <- as.vector(z$residuals) * z$weights
    dispersion <- if(substr(z$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) {
      1
    } else {
      ## for case weights: fix dispersion estimate
      if(!is.null(weights) && caseweights) {
        sum(wres^2/weights, na.rm = TRUE)/sum(z$weights, na.rm = TRUE)
      } else {
        sum(wres^2, na.rm = TRUE)/sum(z$weights, na.rm = TRUE)
      }
    }
    rval$estfun <- wres * x[, !is.na(z$coefficients), drop = FALSE]/dispersion
  }

  ## add model (if desired)
  if(object) {
    class(z) <- c("glm", "lm")
    z$offset <- if(is.null(offset)) 0 else offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- attr(x, "xlevels")    

    cl <- as.call(expression(glm))
    cl$formula <- attr(x, "formula")	
    if(!is.null(offset)) cl$offset <- attr(x, "offset")
    z$call <- cl
    z$terms <- attr(x, "terms")

    ## for case weights: change degrees of freedom
    if(!is.null(weights) && caseweights) {
      z$df.null     <- z$df.null     - sum(weights > 0) + sum(weights)
      z$df.residual <- z$df.residual - sum(weights > 0) + sum(weights)
    }

    rval$object <- z
  }

  return(rval)
}

## methods
print.glmtree <- function(x,
  title = NULL, objfun = NULL, ...)
{
  if(is.null(title)) title <- sprintf("Generalized linear model tree (family: %s)", x$info$family)
  if(is.null(objfun)) objfun <- if(substr(x$info$family, 1L, 5L) != "quasi") "negative log-likelihood" else "deviance"
  print.modelparty(x, title = title, objfun = objfun, ...)
}

predict.glmtree <- function(object, newdata = NULL, type = "response", ...)
{
  ## FIXME: possible to get default?
  if(is.null(newdata) & !identical(type, "node")) stop("newdata has to be provided")
  predict.modelparty(object, newdata = newdata, type = type, ...)
}

plot.glmtree <- function(x, terminal_panel = node_bivplot,
  tp_args = list(), tnex = NULL, drop_terminal = NULL, ...)
{
  nreg <- if(is.null(tp_args$which)) x$info$nreg else length(tp_args$which)
  if(nreg < 1L & missing(terminal_panel)) {
    plot.constparty(as.constparty(x),
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  } else {
    if(is.null(tnex)) tnex <- if(is.null(terminal_panel)) 1L else 2L * nreg
    if(is.null(drop_terminal)) drop_terminal <- !is.null(terminal_panel)
    plot.modelparty(x, terminal_panel = terminal_panel,
      tp_args = tp_args, tnex = tnex, drop_terminal = drop_terminal, ...)
  }
}
