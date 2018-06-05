\name{prune.modelparty}
\alias{prune.modelparty}
\alias{prune.lmtree}

\title{Post-Prune \code{modelparty} Objects}

\usage{
prune.modelparty(tree, type = "AIC", ...)
}

\description{
Post-pruning of \code{modelparty} objects based on information
criteria like AIC, BIC, or related user-defined criteria.
}

\arguments{
  \item{tree}{object of class \code{modelparty}.}
  \item{type}{pruning type. Can be \code{"AIC"}, \code{"BIC"} or a
    user-defined function (details below).}
  \item{\dots}{additional arguments.}
}

\details{
  In \code{\link{mob}}-based model trees, pre-pruning based on p-values
  is used by default and often no post-pruning is necessary in such trees.
  However, if pre-pruning is switched off (by using a large \code{alpha})
  or does is not sufficient (e.g., possibly in large samples) the \code{prune}
  method can be used for subsequent post-pruning based on information criteria.

  The function \code{prune.modelparty} can be called directly but it is also
  registered as a method for the generic \code{\link[rpart]{prune}} function
  from the \pkg{rpart} package. Thus, if \pkg{rpart} is attached,
  \code{prune(tree, type = "AIC", ...)} also works (see examples below).

  To customize the post-pruning strategy,
  \code{type} can be set to a \code{function(objfun, df, nobs)}
  which either returns \code{TRUE} to signal that a current node can be pruned
  or \code{FALSE}. All supplied arguments are of length two: \code{objfun} is the sum of objective
  function values in the current node and its child nodes, respectively.
  \code{df} is the degrees of freedom in the current node and its child nodes,
  respectively. \code{nobs} is vector with the number of observations in the
  current node and the total number of observations in the dataset, respectively. 
  
  For \code{"AIC"} and \code{"BIC"} \code{type} is transformed so that AIC
  or BIC are computed. However, this assumes that the \code{objfun} used in \code{tree}
  is actually the negative log-likelihood. The degrees of freedom assumed for a split 
  can be set via the \code{dfsplit} argument in \code{\link{mob_control}} when computing
  the \code{tree} or manipulated later by changing the value of \code{tree$info$control$dfsplit}.
}

\value{
An object of class \code{modelparty} where the associated tree is either the
same as the original or smaller.
}


\seealso{
\code{\link[rpart]{prune}}, \code{\link{lmtree}}, \code{\link{glmtree}}, \code{\link{mob}}
}
\examples{
set.seed(29)
n <- 1000
d <- data.frame(
  x = runif(n),
  z = runif(n),
  z_noise = factor(sample(1:3, size = n, replace = TRUE))
)
d$y <- rnorm(n, mean = d$x * c(-1, 1)[(d$z > 0.7) + 1], sd = 3)

## glm versus lm / logLik versus sum of squared residuals
fmla <- y ~ x | z + z_noise
lm_big <- lmtree(formula = fmla, data = d, maxdepth = 3, alpha = 1)
glm_big <- glmtree(formula = fmla, data = d, maxdepth = 3, alpha = 1)

AIC(lm_big)
AIC(glm_big)

## load rpart for prune() generic
## (otherwise: use prune.modelparty directly)
if (require("rpart")) {

## pruning
lm_aic <- prune(lm_big, type = "AIC")
lm_bic <- prune(lm_big, type = "BIC")

width(lm_big)
width(lm_aic)
width(lm_bic)

glm_aic <- prune(glm_big, type = "AIC")
glm_bic <- prune(glm_big, type = "BIC")

width(glm_big)
width(glm_aic)
width(glm_bic)

}
}
