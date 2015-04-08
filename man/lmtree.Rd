\name{lmtree}
\alias{lmtree}

\alias{plot.lmtree}
\alias{predict.lmtree}
\alias{print.lmtree}

\title{Linear Model Trees}

\description{
  Model-based recursive partitioning based on least squares
  regression.
}

\usage{
lmtree(formula, data, subset, na.action, weights, offset, cluster, \dots)
}

\arguments{
  \item{formula}{symbolic description of the model (of type
    \code{y ~ z1 + \dots + zl} or \code{y ~ x1 + \dots + xk | z1 + \dots + zl};
    for details see below).}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights. By default these are
    treated as case weights but the default can be changed in
    \code{\link{mob_control}}.}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the model \code{y ~ x1 + \dots + xk} (i.e., only when
    \code{x} variables are specified).}
  \item{cluster}{optional vector (typically numeric or factor) with a
    cluster ID to be employed for clustered covariances in the parameter
    stability tests.}
  \item{\dots}{optional control parameters passed to
    \code{\link{mob_control}}.}
}

\details{
Convenience interface for fitting MOBs (model-based recursive partitions) via
the \code{\link{mob}} function. \code{lmtree} internally sets up a model
\code{fit} function for \code{mob}, using either \code{\link[stats]{lm.fit}}
or \code{\link[stats]{lm.wfit}} (depending on whether weights are used or not).
Then \code{mob} is called using the residual sum of squares as the objective
function.

Compared to calling \code{mob} by hand, the implementation tries to avoid
unnecessary computations while growing the tree. Also, it provides a more
elaborate plotting function.
}

\value{
  An object of class \code{lmtree} inheriting from \code{\link{modelparty}}.
  The \code{info} element of the overall \code{party} and the individual
  \code{node}s contain various informations about the models.
}

\references{ 
Zeileis A, Hothorn T, Hornik K (2008).
  Model-Based Recursive Partitioning.
  \emph{Journal of Computational and Graphical Statistics}, \bold{17}(2), 492--514.
}

\seealso{\code{\link{mob}}, \code{\link{mob_control}}, \code{\link{glmtree}}}

\examples{
if(require("mlbench")) {

## Boston housing data
data("BostonHousing", package = "mlbench")
BostonHousing <- transform(BostonHousing,
  chas = factor(chas, levels = 0:1, labels = c("no", "yes")),
  rad = factor(rad, ordered = TRUE))

## linear model tree
bh_tree <- lmtree(medv ~ log(lstat) + I(rm^2) | zn +
  indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing, minsize = 40)

## printing whole tree or individual nodes
print(bh_tree)
print(bh_tree, node = 7)

## plotting
plot(bh_tree)
plot(bh_tree, tp_args = list(which = "log(lstat)"))
plot(bh_tree, terminal_panel = NULL)

## estimated parameters
coef(bh_tree)
coef(bh_tree, node = 9)
summary(bh_tree, node = 9)

## various ways for computing the mean squared error (on the training data)
mean((BostonHousing$medv - fitted(bh_tree))^2)
mean(residuals(bh_tree)^2)
deviance(bh_tree)/sum(weights(bh_tree))
deviance(bh_tree)/nobs(bh_tree)

## log-likelihood and information criteria
logLik(bh_tree)
AIC(bh_tree)
BIC(bh_tree)
## (Note that this penalizes estimation of error variances, which
## were treated as nuisance parameters in the fitting process.)

## different types of predictions
bh <- BostonHousing[c(1, 10, 50), ]
predict(bh_tree, newdata = bh, type = "node")
predict(bh_tree, newdata = bh, type = "response")
predict(bh_tree, newdata = bh, type = function(object) summary(object)$r.squared)

}


if(require("AER")) {

## Demand for economics journals data
data("Journals", package = "AER")
Journals <- transform(Journals,
  age = 2000 - foundingyear,
  chars = charpp * pages)

## linear regression tree (OLS)
j_tree <- lmtree(log(subs) ~ log(price/citations) | price + citations +
  age + chars + society, data = Journals, minsize = 10, verbose = TRUE)

## printing and plotting
j_tree
plot(j_tree)

## coefficients and summary
coef(j_tree, node = 1:3)
summary(j_tree, node = 1:3)

}


if(require("AER")) {

## Beauty and teaching ratings data
data("TeachingRatings", package = "AER")

## linear regression (WLS)
## null model
tr_null <- lm(eval ~ 1, data = TeachingRatings, weights = students,
  subset = credits == "more")
## main effects
tr_lm <- lm(eval ~ beauty + gender + minority + native + tenure + division,
  data = TeachingRatings, weights = students, subset = credits == "more")
## tree
tr_tree <- lmtree(eval ~ beauty | minority + age + gender + division + native + tenure,
   data = TeachingRatings, weights = students, subset = credits == "more",
   caseweights = FALSE)

## visualization
plot(tr_tree)

## beauty slope coefficient
coef(tr_lm)[2]
coef(tr_tree)[, 2]

## R-squared
1 - deviance(tr_lm)/deviance(tr_null)
1 - deviance(tr_tree)/deviance(tr_null)
}

}
\keyword{tree}
