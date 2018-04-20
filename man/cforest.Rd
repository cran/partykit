\name{cforest}
\alias{cforest}
\alias{gettree}
\alias{gettree.cforest}
\alias{predict.cforest}
\encoding{latin1}

\title{Conditional Random Forests}

\description{
  An implementation of the random forest and bagging ensemble algorithms
  utilizing conditional inference trees as base learners.
}

\usage{
cforest(formula, data, weights, subset, offset, cluster, strata,
        na.action = na.pass,
	control = ctree_control(teststat = "quad", testtype = "Univ",
            mincriterion = 0, saveinfo = FALSE, \dots),
        ytrafo = NULL, scores = NULL, ntree = 500L,
        perturb = list(replace = FALSE, fraction = 0.632),
        mtry = ceiling(sqrt(nvar)), applyfun = NULL, cores = NULL,
        trace = FALSE, \dots)
\method{predict}{cforest}(object, newdata = NULL,
        type = c("response", "prob", "weights", "node"),
        OOB = FALSE, FUN = NULL, simplify = TRUE, scale = TRUE, \dots)
\method{gettree}{cforest}(object, tree = 1L, \dots)
}

\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{
    an optional vector of weights to be used in the fitting
    process. Non-negative integer valued weights are
    allowed as well as non-negative real weights.
    Observations are sampled (with or without replacement)
    according to probabilities \code{weights / sum(weights)}.
    The fraction of observations to be sampled (without replacement)
    is computed based on the sum of the weights if all weights
    are integer-valued and based on the number of weights greater zero
    else. Alternatively, \code{weights} can be a double matrix defining
    case weights for all \code{ncol(weights)} trees in the forest directly.
    This requires more storage but gives the user more control.}
  \item{offset}{ an optional vector of offset values.}
  \item{cluster}{ an optional factor indicating independent clusters.
                  Highly experimental, use at your own risk.}
  \item{strata}{ an optional factor for stratified sampling.}
  \item{na.action}{a function which indicates what should happen when the data
                   contain missing value.}
  \item{control}{a list with control parameters, see
                 \code{\link{ctree_control}}. The default values correspond to those
                 of the default values used by \code{\link[party]{cforest}} from the
                 \code{party} package. \code{saveinfo = FALSE} leads to less
                 memory hungry representations of trees. Note that arguments
                 \code{mtry}, \code{cores} and \code{applyfun} in
                 \code{\link{ctree_control}} are ignored for \code{\link{cforest}},
                 because they are already set.}
  \item{ytrafo}{an optional named list of functions to be applied to the response
                variable(s) before testing their association with the explanatory
                variables. Note that this transformation is only
                performed once for the root node and does not take weights into account
                (which means, the forest bootstrap or subsetting is ignored,
                which is almost certainly not a good idea).
                Alternatively, \code{ytrafo} can be a function of \code{data} and
                \code{weights}. In this case, the transformation is computed for
                every node and the corresponding weights. This feature is experimental
                and the user interface likely to change.}
  \item{scores}{an optional named list of scores to be attached to ordered
                factors.}
  \item{ntree}{
    Number of trees to grow for the forest.}
  \item{perturb}{
    a list with arguments \code{replace} and \code{fraction} determining which type of
    resampling with \code{replace = TRUE} referring to the n-out-of-n bootstrap and
    \code{replace = FALSE} to sample splitting. \code{fraction} is the number of observations to draw without
    replacement. }
  \item{mtry}{
    number of input variables randomly sampled as candidates
    at each node for random forest like algorithms. Bagging, as special case
    of a random forest without random input variable sampling, can
    be performed by setting \code{mtry} either equal to \code{Inf} or
    manually equal to the number of input variables.}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
                  \code{function(X, FUN, \dots)}. It is used for computing the variable selection criterion.
                  The default is to use the basic \code{lapply}
                  function unless the \code{cores} argument is specified (see below).}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
               \code{\link[parallel]{mclapply}} with the desired number of \code{cores}.}

  \item{trace}{a logical indicating if a progress bar shall be printed while
          the forest grows.}
  \item{object}{
    An object as returned by \code{cforest}}
  \item{newdata}{
    An optional data frame containing test data.}
  \item{type}{
    a character string denoting the type of predicted value
          returned, ignored when argument \code{FUN} is given.  For
          \code{"response"}, the mean of a numeric response, the predicted
          class for a categorical response or the median survival time for a
          censored response is returned.  For \code{"prob"} the matrix of
          conditional class probabilities (\code{simplify = TRUE}) or a list
          with the conditional class probabilities for each observation
          (\code{simplify = FALSE}) is returned for a categorical response.
          For numeric and censored responses, a list with the empirical
          cumulative distribution functions and empirical survivor functions
          (Kaplan-Meier estimate) is returned when \code{type = "prob"}.
          \code{"weights"} returns an integer vector of prediction weights.
          For \code{type = "where"}, a list of terminal node ids for each of the
          trees in the forest ist returned.}
  \item{OOB}{
    a logical defining out-of-bag predictions (only if \code{newdata = NULL}).}
  \item{FUN}{
    a function to compute summary statistics. Predictions for each node have to be
    computed based on arguments \code{(y, w)} where \code{y} is the response and
    \code{w} are case weights.}
  \item{simplify}{
    a logical indicating whether the resulting list
                   of predictions should be converted to a suitable
                   vector or matrix (if possible).}
  \item{scale}{a logical indicating scaling of the nearest neighbor weights
               by the sum of weights in the corresponding terminal node of
               each tree. In the simple regression forest, predicting
               the conditional mean by nearest neighbor weights will be
               equivalent to (but slower!) the aggregation of means.}
  \item{tree}{ an integer, the number of the tree to extract from the forest.}
  \item{\dots}{ additional arguments. }
}

\details{
  This implementation of the random forest (and bagging) algorithm differs
  from the reference implementation in \code{\link[randomForest]{randomForest}}
  with respect to the base learners used and the aggregation scheme applied.

  Conditional inference trees, see \code{\link{ctree}}, are fitted to each
  of the \code{ntree} perturbed samples of the learning sample. Most of the hyper parameters in
  \code{\link{ctree_control}} regulate the construction of the conditional inference trees.

  Hyper parameters you might want to change are:

  1. The number of randomly preselected variables \code{mtry}, which is fixed
     to the square root of the number of input variables.

  2. The number of trees \code{ntree}. Use more trees if you have more variables.

  3. The depth of the trees, regulated by \code{mincriterion}. Usually unstopped and unpruned
     trees are used in random forests. To grow large trees, set \code{mincriterion} to a small value.


  The aggregation scheme works by averaging observation weights extracted
  from each of the \code{ntree} trees and NOT by averaging predictions directly
  as in \code{\link[randomForest]{randomForest}}.
  See Hothorn et al. (2004) and Meinshausen (2006) for a description.

  Predictions can be computed using \code{\link{predict}}. For observations
  with zero weights, predictions are computed from the fitted tree
  when \code{newdata = NULL}.

  Ensembles of conditional inference trees have not yet been extensively
  tested, so this routine is meant for the expert user only and its current
  state is rather experimental. However, there are some things available
  in \code{\link{cforest}} that can't be done with \code{\link[randomForest]{randomForest}},
  for example fitting forests to censored response variables (see Hothorn et al., 2004, 2006a) or to
  multivariate and ordered responses. Using the rich \code{partykit} infrastructure allows
  additional functionality in \code{cforest}, such as parallel tree growing and probabilistic
  forecasting (for example via quantile regression forests). Also plotting of single trees from
  a forest is much easier now.

  Unlike \code{\link[party]{cforest}}, \code{cforest} is entirely written in R which
  makes customisation much easier at the price of longer computing times. However, trees
  can be grown in parallel with this R only implemention which renders speed less of an issue.
  Note that the default values are different from those used in package \code{party}, most
  importantly the default for mtry is now data-dependent. \code{predict(, type = "node")} replaces
  the \code{\link[party]{where}} function and \code{predict(, type = "prob")} the
  \code{\link[party]{treeresponse}} function.

  Moreover, when predictors vary in their scale of measurement of number
  of categories, variable selection and computation of variable importance is biased
  in favor of variables with many potential cutpoints in \code{\link[randomForest]{randomForest}},
  while in \code{\link{cforest}} unbiased trees and an adequate resampling scheme
  are used by default. See Hothorn et al. (2006b) and Strobl et al. (2007)
  as well as Strobl et al. (2009).
}

\value{
  An object of class \code{cforest}.
}

\references{
Breiman L (2001).
  Random Forests.
  \emph{Machine Learning}, \bold{45}(1), 5--32.

Hothorn T, Lausen B, Benner A, Radespiel-Troeger M (2004).
  Bagging Survival Trees.
  \emph{Statistics in Medicine}, \bold{23}(1), 77--91.

Hothorn T, B?hlmann P, Dudoit S, Molinaro A, Van der Laan MJ (2006a).
  Survival Ensembles.
  \emph{Biostatistics}, \bold{7}(3), 355--373.

Hothorn T, Hornik K, Zeileis A (2006b).
  Unbiased Recursive Partitioning: A Conditional Inference Framework.
  \emph{Journal of Computational and Graphical Statistics}, \bold{15}(3), 651--674.

Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.

Meinshausen N (2006).
  Quantile Regression Forests.
  \emph{Journal of Machine Learning Research}, \bold{7}, 983--999.

Strobl C, Boulesteix AL, Zeileis A, Hothorn T (2007).
  Bias in Random Forest Variable Importance Measures: Illustrations, Sources and a Solution.
  \emph{BMC Bioinformatics}, \bold{8}, 25.
  \url{http://www.biomedcentral.com/1471-2105/8/25}

Strobl C, Malley J, Tutz G (2009).
  An Introduction to Recursive Partitioning: Rationale, Application, and Characteristics of
  Classification and Regression Trees, Bagging, and Random Forests.
  \emph{Psychological Methods}, \bold{14}(4), 323--348.
}

\examples{
## basic example: conditional inference forest for cars data
cf <- cforest(dist ~ speed, data = cars)

## prediction of fitted mean and visualization
nd <- data.frame(speed = 4:25)
nd$mean  <- predict(cf, newdata = nd, type = "response")
plot(dist ~ speed, data = cars)
lines(mean ~ speed, data = nd)

## predict quantiles (aka quantile regression forest)
myquantile <- function(y, w) quantile(rep(y, w), probs = c(0.1, 0.5, 0.9))
p <- predict(cf, newdata = nd, type = "response", FUN = myquantile)
colnames(p) <- c("lower", "median", "upper")
nd <- cbind(nd, p)

## visualization with conditional (on speed) prediction intervals
plot(dist ~ speed, data = cars, type = "n")
with(nd, polygon(c(speed, rev(speed)), c(lower, rev(upper)),
  col = "lightgray", border = "transparent"))
points(dist ~ speed, data = cars)
lines(mean ~ speed, data = nd, lwd = 1.5)
lines(median ~ speed, data = nd, lty = 2, lwd = 1.5)
legend("topleft", c("mean", "median", "10\% - 90\% quantile"),
  lwd = c(1.5, 1.5, 10), lty = c(1, 2, 1),
  col = c("black", "black", "lightgray"), bty = "n")

### we may also use predicted conditional (on speed) densities
mydensity <- function (y, w) approxfun(density(y, weights = w/sum(w))[1:2], rule = 2)
pd <- predict(cf, newdata = nd, type = "response", FUN = mydensity)

## visualization in heatmap (instead of scatterplot)
## with fitted curves as above
dist <- -10:150
dens <- t(sapply(seq_along(pd), function(i) pd[[i]](dist)))
image(nd$speed, dist, dens, xlab = "speed", col = rev(gray.colors(9)))
lines(mean ~ speed, data = nd, lwd = 1.5)
lines(median ~ speed, data = nd, lty = 2, lwd = 1.5)
lines(lower ~ speed, data = nd, lty = 2)
lines(upper ~ speed, data = nd, lty = 2)

\dontrun{

### honest (i.e., out-of-bag) cross-classification of
### true vs. predicted classes
data("mammoexp", package = "TH.data")
table(mammoexp$ME, predict(cforest(ME ~ ., data = mammoexp, ntree = 50),
                           OOB = TRUE, type = "response"))

### fit forest to censored response
if (require("TH.data") && require("survival")) {

    data("GBSG2", package = "TH.data")
    bst <- cforest(Surv(time, cens) ~ ., data = GBSG2, ntree = 50)

    ### estimate conditional Kaplan-Meier curves
    print(predict(bst, newdata = GBSG2[1:2,], OOB = TRUE, type = "prob"))

    print(gettree(bst))
}
}
}

\keyword{tree}
