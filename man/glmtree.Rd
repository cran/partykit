\name{glmtree}
\alias{glmtree}

\alias{plot.glmtree}
\alias{predict.glmtree}
\alias{print.glmtree}

\title{Generalized Linear Model Trees}

\description{
  Model-based recursive partitioning based on generalized linear models.
}

\usage{
glmtree(formula, data, subset, na.action, weights, offset, cluster,
  family = gaussian, epsilon = 1e-8, maxit = 25, \dots)
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
  \item{family}{specification of a family for \code{\link[stats]{glm}}.}
  \item{epsilon, maxit}{control parameters passed to
    \code{\link[stats]{glm.control}}.}
  \item{\dots}{optional control parameters passed to
    \code{\link{mob_control}}.}
}

\details{
Convenience interface for fitting MOBs (model-based recursive partitions) via
the \code{\link{mob}} function. \code{glmtree} internally sets up a model
\code{fit} function for \code{mob}, using \code{\link[stats]{glm.fit}}.
Then \code{mob} is called using the negative log-likelihood as the objective
function.

Compared to calling \code{mob} by hand, the implementation tries to avoid
unnecessary computations while growing the tree. Also, it provides a more
elaborate plotting function.
}

\value{
  An object of class \code{glmtree} inheriting from \code{\link{modelparty}}.
  The \code{info} element of the overall \code{party} and the individual
  \code{node}s contain various informations about the models.
}

\references{ 
Zeileis A, Hothorn T, Hornik K (2008).
  Model-Based Recursive Partitioning.
  \emph{Journal of Computational and Graphical Statistics}, \bold{17}(2), 492--514.
}

\seealso{\code{\link{mob}}, \code{\link{mob_control}}, \code{\link{lmtree}}}

\examples{
if(require("mlbench")) {

## Pima Indians diabetes data
data("PimaIndiansDiabetes", package = "mlbench")

## recursive partitioning of a logistic regression model
pid_tree2 <- glmtree(diabetes ~ glucose | pregnant +
  pressure + triceps + insulin + mass + pedigree + age,
  data = PimaIndiansDiabetes, family = binomial)

## printing whole tree or individual nodes
print(pid_tree2)
print(pid_tree2, node = 1)

## visualization
plot(pid_tree2)
plot(pid_tree2, tp_args = list(cdplot = TRUE))
plot(pid_tree2, terminal_panel = NULL)

## estimated parameters
coef(pid_tree2)
coef(pid_tree2, node = 5)
summary(pid_tree2, node = 5)

## deviance, log-likelihood and information criteria
deviance(pid_tree2)
logLik(pid_tree2)
AIC(pid_tree2)
BIC(pid_tree2)

## different types of predictions
pid <- head(PimaIndiansDiabetes)
predict(pid_tree2, newdata = pid, type = "node")
predict(pid_tree2, newdata = pid, type = "response")
predict(pid_tree2, newdata = pid, type = "link")

}
}
\keyword{tree}
