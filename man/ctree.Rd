\name{ctree}
\alias{ctree}
\alias{sctest.constparty}

\title{Conditional Inference Trees}

\description{
  Recursive partitioning for continuous, censored, ordered, nominal and
  multivariate response variables in a conditional inference framework. 
}

\usage{
ctree(formula, data, subset, weights, na.action = na.pass, offset, cluster, 
    control = ctree_control(\dots), ytrafo = NULL, 
    converged = NULL, scores = NULL, doFit = TRUE, \dots)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Only non-negative integer valued weights are
                  allowed.}
  \item{offset}{ an optional vector of offset values.}
  \item{cluster}{ an optional factor indicating independent clusters. 
                  Highly experimental, use at your own risk.}
  \item{na.action}{a function which indicates what should happen when the data
                   contain missing value.}
  \item{control}{a list with control parameters, see
                 \code{\link{ctree_control}}.}
  \item{ytrafo}{an optional named list of functions to be applied to the response 
                variable(s) before testing their association with the explanatory 
                variables. Note that this transformation is only
                performed once for the root node and does not take weights into account.
                Alternatively, \code{ytrafo} can be a function of \code{data} and
                \code{weights}. In this case, the transformation is computed for
                every node with corresponding weights. This feature is experimental
                and the user interface likely to change.}
  \item{converged}{an optional function for checking user-defined criteria
                   before splits are implemented. This is not to be used and
                   very likely to change.}
  \item{scores}{an optional named list of scores to be attached to ordered
                factors.}
  \item{doFit}{a logical, if \code{FALSE}, the tree is not fitted.}
  \item{\dots}{arguments passed to \code{\link{ctree_control}}.}
}

\details{
  Function \code{partykit::ctree} is a reimplementation of (most of)
  \code{party::ctree} employing the new \code{\link{party}} infrastructure
  of the \pkg{partykit} infrastructure. The vignette \code{vignette("ctree", package = "partykit")}
  explains internals of the different implementations.
  
  Conditional inference trees estimate a regression relationship by binary recursive
  partitioning in a conditional inference framework. Roughly, the algorithm
  works as follows: 1) Test the global null hypothesis of independence between
  any of the input variables and the response (which may be multivariate as well). 
  Stop if this hypothesis cannot be rejected. Otherwise select the input
  variable with strongest association to the response. This
  association is measured by a p-value corresponding to a test for the
  partial null hypothesis of a single input variable and the response.
  2) Implement a binary split in the selected input variable. 
  3) Recursively repeate steps 1) and 2). 

  The implementation utilizes a unified framework for conditional inference,
  or permutation tests, developed by Strasser and Weber (1999). The stop
  criterion in step 1) is either based on multiplicity adjusted p-values 
  (\code{testtype = "Bonferroni"} in \code{\link{ctree_control}})
  or on the univariate p-values (\code{testtype = "Univariate"}). In both cases, the
  criterion is maximized, i.e., 1 - p-value is used. A split is implemented 
  when the criterion exceeds the value given by \code{mincriterion} as
  specified in \code{\link{ctree_control}}. For example, when 
  \code{mincriterion = 0.95}, the p-value must be smaller than
  $0.05$ in order to split this node. This statistical approach ensures that
  the right-sized tree is grown without additional (post-)pruning or cross-validation.
  The level of \code{mincriterion} can either be specified to be appropriate
  for the size of the data set (and \code{0.95} is typically appropriate for
  small to moderately-sized data sets) or could potentially be treated like a
  hyperparameter (see Section~3.4 in Hothorn, Hornik and Zeileis, 2006).
  The selection of the input variable to split in
  is based on the univariate p-values avoiding a variable selection bias
  towards input variables with many possible cutpoints. The test statistics
  in each of the nodes can be extracted with the \code{sctest} method.
  (Note that the generic is in the \pkg{strucchange} package so this either
  needs to be loaded or \code{sctest.constparty} has to be called directly.)
  In cases where splitting stops due to the sample size (e.g., \code{minsplit}
  or \code{minbucket} etc.), the test results may be empty.

  Predictions can be computed using \code{\link{predict}}, which returns predicted means,
  predicted classes or median predicted survival times and 
  more information about the conditional
  distribution of the response, i.e., class probabilities
  or predicted Kaplan-Meier curves. For observations
  with zero weights, predictions are computed from the fitted tree 
  when \code{newdata = NULL}.

  By default, the scores for each ordinal factor \code{x} are
  \code{1:length(x)}, this may be changed for variables in the formula 
  using \code{scores = list(x = c(1, 5, 6))}, for example.

  For a general description of the methodology see Hothorn, Hornik and
  Zeileis (2006) and Hothorn, Hornik, van de Wiel and Zeileis (2006).
}

\value{
  An object of class \code{\link{party}}.
}

\references{ 
Hothorn T, Hornik K, Van de Wiel MA, Zeileis A (2006).
  A Lego System for Conditional Inference.
  \emph{The American Statistician}, \bold{60}(3), 257--263.

Hothorn T, Hornik K, Zeileis A (2006).
  Unbiased Recursive Partitioning: A Conditional Inference Framework.
  \emph{Journal of Computational and Graphical Statistics}, \bold{15}(3), 651--674. 

Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.

Strasser H, Weber C (1999).
  On the Asymptotic Theory of Permutation Statistics.
  \emph{Mathematical Methods of Statistics}, \bold{8}, 220--250.
}

\examples{
### regression
airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq)
airct
plot(airct)
mean((airq$Ozone - predict(airct))^2)

### classification
irisct <- ctree(Species ~ .,data = iris)
irisct
plot(irisct)
table(predict(irisct), iris$Species)

### estimated class probabilities, a list
tr <- predict(irisct, newdata = iris[1:10,], type = "prob")

### survival analysis
if (require("TH.data") && require("survival") && 
    require("coin") && require("Formula")) {

  data("GBSG2", package = "TH.data")
  (GBSG2ct <- ctree(Surv(time, cens) ~ ., data = GBSG2))
  predict(GBSG2ct, newdata = GBSG2[1:2,], type = "response")	  
  plot(GBSG2ct)

  ### with weight-dependent log-rank scores
  ### log-rank trafo for observations in this node only (= weights > 0)
  h <- function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
      if (is.null(weights)) weights <- rep(1, NROW(y))
      s <- logrank_trafo(y[weights > 0,,drop = FALSE])
      r <- rep(0, length(weights))
      r[weights > 0] <- s
      list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }

  ### very much the same tree
  (ctree(Surv(time, cens) ~ ., data = GBSG2, ytrafo = h))
}

### multivariate responses
airct2 <- ctree(Ozone + Temp ~ ., data = airq)
airct2
plot(airct2)
}
\keyword{tree}
