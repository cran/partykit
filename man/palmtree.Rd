\name{palmtree}
\alias{palmtree}
\alias{print.palmtree}
\alias{plot.palmtree}
\alias{coef.palmtree}
\alias{logLik.palmtree}
\alias{predict.palmtree}

\title{Partially Additive (Generalized) Linear Model Trees}

\description{
  Model-based recursive partitioning based on (generalized)
  linear models with some local (i.e., leaf-specific) and
  some global (i.e., constant throughout the tree) regression
  coefficients.
}

\usage{
palmtree(formula, data, weights = NULL, family = NULL,
  lmstart = NULL, abstol = 0.001, maxit = 100, 
  dfsplit = TRUE, verbose = FALSE, plot = FALSE, ...)
}

\arguments{
  \item{formula}{formula specifying the response variable and a three-part
    right-hand-side describing the local (i.e., leaf-specific) regressors,
    the global regressors (i.e., with constant coefficients throughout the tree),
    and partitioning variables, respectively. For details see below.}
  \item{data}{data.frame to be used for estimating the model tree.}
  \item{weights}{numeric. An optional numeric vector of weights. (Note that
    this is passed with standard evaluation, i.e., it is not enough to pass
    the name of a column in \code{data}.)}
  \item{family}{either \code{NULL} so that \code{lm}/\code{lmtree} are used
    or family specification for \code{glm}/\code{glmtree}.
    See \code{\link[stats]{glm}} documentation for families.}
  \item{lmstart}{numeric. A vector of length \code{nrow(data)}, to be used as
    an offset in estimation of the first tree. \code{NULL} by default,
    which results in an initialization with the global model.}
  \item{abstol}{numeric. The convergence criterion used for estimation of the model.
    When the difference in log-likelihoods of the model from two consecutive
    iterations is smaller than \code{abstol}, estimation of the model tree has
    converged.} 
  \item{maxit}{numeric. The maximum number of iterations to be performed in
    estimation of the model tree.}
  \item{dfsplit}{logical or numeric. \code{as.integer(dfsplit)} is the
    degrees of freedom per selected split employed when extracting
    the log-likelihood.}
  \item{verbose}{Should the log-likelihood value of the estimated
    model be printed for every iteration of the estimation?}
  \item{plot}{Should the tree be plotted at every iteration of the estimation?
    Note that selecting this option slows down execution of the function.}
  \item{\dots}{Additional arguments to be passed to \code{lmtree()} or \code{glmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  
\details{
  Partially additive (generalized) linear model (PALM) trees learn a tree
  where each terminal node is associated with different regression coefficients
  while adjusting for additional global regression effects. This allows for
  detection of subgroup-specific coefficients with respect to selected covariates,
  while keeping the remaining regression coefficients constant throughout the tree.
  The estimation algorithm iterates between (1) estimation of the tree given
  an offset of the global effects, and (2) estimation of the global regression effects
  given the tree structure.
  
  To specify all variables in the model a \code{formula} such as
  \code{y ~ x1 + x2 | x3 | z1 + z2 + z3} is used, where \code{y} is the
  response, \code{x1} and \code{x2} are the regressors in every node of the
  tree, \code{x3} has a global regression coefficients, and \code{z1} to \code{z3}
  are the partitioning variables considered for growing the tree.

  The code is still under development and might change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{call}{the matched call.}
  \item{tree}{The final \code{lmtree}/\code{glmtree}.}
  \item{palm}{The final \code{lm}/\code{glm} model.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.ranef} and \code{.treeresponse}
    from the last iteration.}
  \item{nobs}{Number of observations.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{df}{Degrees of freedom.}
  \item{dfsplit}{degrees of freedom per selected split as specified with the \code{dfsplit}
    argument.}
  \item{iterations}{The number of iterations used to estimate the \code{palmtree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{lmstart}{Offset in estimation of the first tree as specified in the \code{lmstart} 
    argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{intercept}{Logical specifying if an intercept was computed.}
  \item{family}{The \code{family} object used.}
  \item{mob.control}{A list containing control parameters passed to
    \code{lmtree()}, as specified with \dots.}
}

\seealso{\code{\link[stats]{lm}}, \code{\link[stats]{glm}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\references{ 
Sies A, Van Mechelen I (2015).
  Comparing Four Methods for Estimating Tree-Based Treatment Regimes.
  \emph{Unpublished Manuscript}.
}

\examples{
## one DGP from Sies and Van Mechelen (2015)
dgp <- function(nobs = 1000, nreg = 5, creg = 0.4, ptreat = 0.5, sd = 1,
  coef = c(1, 0.25, 0.25, 0, 0, -0.25), eff = 1)
{
  d <- mvtnorm::rmvnorm(nobs,
    mean = rep(0, nreg),
    sigma = diag(1 - creg, nreg) + creg)
  colnames(d) <- paste0("x", 1:nreg)
  d <- as.data.frame(d)
  d$a <- rbinom(nobs, size = 1, prob = ptreat)
  d$err <- rnorm(nobs, mean = 0, sd = sd)

  gopt <- function(d) {
    as.numeric(d$x1 > -0.545) * as.numeric(d$x2 < 0.545)
  }
  d$y <- coef[1] + drop(as.matrix(d[, paste0("x", 1:5)]) \%*\% coef[-1]) -
    eff * (d$a - gopt(d))^2 + d$err
  d$a <- factor(d$a)
  return(d)
}
set.seed(1)
d <- dgp()

## estimate PALM tree with correctly specified global (partially
## additive) regressors and all variables considered for partitioning
palm <- palmtree(y ~ a | x1 + x2 + x5 | x1 + x2 + x3 + x4 + x5, data = d)
print(palm)
plot(palm)

## query coefficients
coef(palm, model = "tree")
coef(palm, model = "palm")
coef(palm, model = "all")
}

\keyword{tree}
