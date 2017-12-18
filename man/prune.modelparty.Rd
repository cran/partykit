\name{prune.modelparty}
\alias{prune.modelparty}

\title{Post prune \code{modelparty} objects.}

\usage{
\method{prune}{modelparty}(tree, type = "AIC", ...)
}

\description{
This function allows for post pruning of modelparty objects. 
Not thate this is usually not neccessary if p-values are used for stopping.
}

\arguments{
  \item{tree}{object of class \code{modelparty}.}
  \item{type}{pruning type. Can be \code{"AIC"}, \code{"BIC"} or a user 
    defined function (details below).}
  \item{\dots}{additional arguments.}
}

\details{
  \code{type} can be set to a \code{function(objfun, df, nobs)}
  which either returns \code{TRUE} to signal that a current node can be pruned
  or \code{FALSE}. All supplied arguments are of length two: \code{objfun} is the sum of objective
  function values in the current node and its child nodes, respectively.
  \code{df} is the degrees of freedom in the current node and its child nodes,
  respectively. \code{nobs} is vector with the number of observations in the
  current node and the total number of observations in the dataset, respectively.  
}

\value{
An object of class \code{modelparty} that is smaller or equal the original.
}

\author{
Heidi Seibold
}


\seealso{
\code{\link{lmtree}}, \code{\link{glmtree}}, \code{\link{mob}}
}
\examples{
set.seed(29)
n <- 1000
x <- runif(n)
z <- runif(n)
y <- rnorm(n, mean = x * c(-1, 1)[(z > 0.7) + 1], sd = 3)
z_noise <- factor(sample(1:3, size = n, replace = TRUE))
d <- data.frame(y = y, x = x, z = z, z_noise = z_noise)

fmla <- as.formula("y ~ x | z + z_noise")


## glm versus lm / logLik versus sum of squared residuals
lm_big <- lmtree(formula = fmla, data = d, maxdepth = 2, alpha = 1)
glm_big <- glmtree(formula = fmla, data = d, maxdepth = 2, alpha = 1)

AIC(lm_big)
AIC(glm_big)


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
