\name{mob}
\alias{mob}
\alias{modelparty}

\alias{coef.modelparty}
\alias{deviance.modelparty}
\alias{fitted.modelparty}
\alias{formula.modelparty}
\alias{getCall.modelparty}
\alias{logLik.modelparty}
\alias{model.frame.modelparty}
\alias{nobs.modelparty}
\alias{plot.modelparty}
\alias{predict.modelparty}
\alias{print.modelparty}
\alias{residuals.modelparty}
\alias{summary.modelparty}
\alias{weights.modelparty}

\alias{refit.modelparty}
\alias{sctest.modelparty}

\title{Model-based Recursive Partitioning}

\description{
  MOB is an algorithm for model-based recursive partitioning yielding
  a tree with fitted models associated with each terminal node.
}

\usage{
mob(formula, data, subset, na.action, weights, offset, cluster,
  fit, control = mob_control(), \dots)
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
    cluster ID to be passed on to the \code{fit} function and employed
    for clustered covariances in the parameter stability tests.}
  \item{fit}{function. A function for fitting the model within each node.
    For details see below.}
  \item{control}{A list with control parameters as returned by
    \code{\link{mob_control}}.}
  \item{\dots}{Additional arguments passed to the \code{fit} function.}
}

\details{
Model-based partitioning fits a model tree using two groups of variables:
(1) The model variables which can be just a (set of) response(s) \code{y} or
additionally include regressors \code{x1}, \dots, \code{xk}. These are used
for estimating the model parameters.
(2) Partitioning variables \code{z1}, \dots, \code{zl}, which are used for
recursively partitioning the data. The two groups of variables are either specified
as \code{y ~ z1 + \dots + zl} (when there are no regressors) or
\code{y ~ x1 + \dots + xk | z1 + \dots + zl} (when the model part contains regressors).
Both sets of variables may in principle be overlapping.

To fit a tree model the following algorithm is used.
  \enumerate{
    \item \code{fit} a model to the \code{y} or \code{y} and \code{x} variables
          using the observations in the current node
    \item Assess the stability of the model parameters with respect to each
          of the partitioning variables \code{z1}, \dots, \code{zl}. If
	  there is some overall instability, choose the variable \code{z}
	  associated with the smallest \eqn{p} value for partitioning, otherwise
	  stop.
    \item Search for the locally optimal split in \code{z} by minimizing the
          objective function of the model. Typically, this will be
	  something like \code{\link{deviance}} or the negative \code{\link{logLik}}.
    \item Refit the \code{model} in both kid subsamples and repeat from step 2.
  }
  
More details on the conceptual design of the algorithm can be found in 
Zeileis, Hothorn, Hornik (2008) and some illustrations are provided in
\code{vignette("MOB")}.
For specifying the \code{fit} function two approaches are possible:

(1) It can be a function \code{fit(y, x = NULL, start = NULL, weights = NULL,
offset = NULL, \dots)}. The arguments \code{y}, \code{x}, \code{weights}, \code{offset}
will be set to the corresponding elements in the current node of the tree.
Additionally, starting values will sometimes be supplied via \code{start}.
Of course, the \code{fit} function can choose to ignore any arguments that are
not applicable, e.g., if the are no regressors \code{x} in the model or if
starting values or not supported. The returned object needs to have a class
that has associated \code{\link[stats]{coef}}, \code{\link[stats]{logLik}}, and
\code{\link[sandwich]{estfun}} methods for extracting the estimated parameters,
the maximized log-likelihood, and the empirical estimating function (i.e.,
score or gradient contributions), respectively.

(2) It can be a function \code{fit(y, x = NULL, start = NULL, weights = NULL,
offset = NULL, \dots, estfun = FALSE, object = FALSE)}. The arguments have the
same meaning as above but the returned object needs to have a different structure.
It needs to be a list with elements \code{coefficients} (containing the estimated
parameters), \code{objfun} (containing the minimized objective function),
\code{estfun} (the empirical estimating functions), and \code{object} (the
fitted model object). The elements \code{estfun}, or \code{object} should be
\code{NULL} if the corresponding argument is set to \code{FALSE}.

Internally, a function of type (2) is set up by \code{mob()} in case a function
of type (1) is supplied. However, to save computation time, a function of type
(2) may also be specified directly.

For the fitted MOB tree, several standard methods are provided such as
\code{print}, \code{predict}, \code{residuals}, \code{logLik}, \code{deviance},
\code{weights}, \code{coef} and \code{summary}. Some of these rely on reusing the
corresponding methods for the individual model objects in the terminal nodes.
Functions such as \code{coef}, \code{print}, \code{summary} also take a
\code{node} argument that can specify the node IDs to be queried.
Some examples are given below.

More details can be found in \code{vignette("mob", package = "partykit")}.
An overview of the connections to other functions in the package is provided
by Hothorn and Zeileis (2015).
}

\value{
  An object of class \code{modelparty} inheriting from \code{\link{party}}.
  The \code{info} element of the overall \code{party} and the individual
  \code{node}s contain various informations about the models.
}

\references{ 
Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.

Zeileis A, Hothorn T, Hornik K (2008).
  Model-Based Recursive Partitioning.
  \emph{Journal of Computational and Graphical Statistics}, \bold{17}(2), 492--514.
}

\seealso{\code{\link{mob_control}}, \code{\link{lmtree}}, \code{\link{glmtree}}}

\examples{
if(require("mlbench")) {

## Pima Indians diabetes data
data("PimaIndiansDiabetes", package = "mlbench")

## a simple basic fitting function (of type 1) for a logistic regression
logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  glm(y ~ 0 + x, family = binomial, start = start, ...)
}

## set up a logistic regression tree
pid_tree <- mob(diabetes ~ glucose | pregnant + pressure + triceps + insulin +
  mass + pedigree + age, data = PimaIndiansDiabetes, fit = logit)
## see lmtree() and glmtree() for interfaces with more efficient fitting functions

## print tree
print(pid_tree)

## print information about (some) nodes
print(pid_tree, node = 3:4)

## visualization
plot(pid_tree)

## coefficients and summary
coef(pid_tree)
coef(pid_tree, node = 1)
summary(pid_tree, node = 1)

## average deviance computed in different ways
mean(residuals(pid_tree)^2)
deviance(pid_tree)/sum(weights(pid_tree))
deviance(pid_tree)/nobs(pid_tree)

## log-likelihood and information criteria
logLik(pid_tree)
AIC(pid_tree)
BIC(pid_tree)

## predicted nodes
predict(pid_tree, newdata = head(PimaIndiansDiabetes, 6), type = "node")
## other types of predictions are possible using lmtree()/glmtree()
}
}
\keyword{tree}
