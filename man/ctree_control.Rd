\name{ctree_control}
\alias{ctree_control}
\title{ Control for Conditional Inference Trees }
\description{

  Various parameters that control aspects of the `ctree' fit.

}
\usage{
ctree_control(teststat = c("quadratic", "maximum"),
    splitstat = c("quadratic", "maximum"),
    splittest = FALSE,
    testtype = c("Bonferroni", "MonteCarlo", "Univariate", "Teststatistic"),
    pargs = GenzBretz(),
    nmax = c(yx = Inf, z = Inf), alpha = 0.05, mincriterion = 1 - alpha,
    logmincriterion = log(mincriterion), minsplit = 20L, minbucket = 7L,
    minprob = 0.01, stump = FALSE, lookahead = FALSE, MIA = FALSE, nresample = 9999L,
    tol = sqrt(.Machine$double.eps),maxsurrogate = 0L, numsurrogate = FALSE,
    mtry = Inf, maxdepth = Inf,
    multiway = FALSE, splittry = 2L, intersplit = FALSE, majority = FALSE,
    caseweights = TRUE, applyfun = NULL, cores = NULL, saveinfo = TRUE,
    update = NULL, splitflavour = c("ctree", "exhaustive"))
}
\arguments{
  \item{teststat}{ a character specifying the type of the test statistic
                       to be applied for variable selection. }
  \item{splitstat}{ a character specifying the type of the test statistic
                       to be applied for splitpoint selection. Prior to
                       version 1.2-0, \code{maximum} was implemented only.}
  \item{splittest}{ a logical changing linear (the default \code{FALSE}) to
                    maximally selected statistics for
                    variable selection. Currently needs \code{testtype = "MonteCarlo"}.}
  \item{testtype}{ a character specifying how to compute the distribution of
                   the test statistic. The first three options refer to
                   p-values as criterion, \code{Teststatistic} uses the raw
                   statistic as criterion. \code{Bonferroni} and
                   \code{Univariate} relate to p-values from the asymptotic
                   distribution (adjusted or unadjusted).
                   Bonferroni-adjusted Monte-Carlo p-values are computed
                   when both \code{Bonferroni} and \code{MonteCarlo} are
                   given.}
  \item{pargs}{ control parameters for the computation of multivariate
                normal probabilities, see \code{\link[mvtnorm]{GenzBretz}}.}
  \item{nmax}{ an integer of length two defining the number of bins each variable
               (in the response \code{yx} and the partitioning variables
               \code{z})) and is divided into prior to tree building. The default \code{Inf}
               does not apply any binning. Highly experimental, use at your own
               risk.}
  \item{alpha}{ a double, the significance level for variable selection.}
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to implement a split. }
  \item{logmincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to implement a split on
                       the log-scale. }
  \item{minsplit}{ the minimum sum of weights in a node in order to be considered
                   for splitting. }
  \item{minbucket}{ the minimum sum of weights in a terminal node. }
  \item{minprob}{ proportion of observations needed to establish a terminal node.}
  \item{stump}{ a logical determining whether a stump (a tree with a maximum of three
                nodes only) is to be computed. }
  \item{lookahead}{ a logical determining whether a split is implemented only
                    after checking if tests in both daughter nodes can be performed.}
  \item{MIA}{ a logical determining the treatment of \code{NA} as a category in split,
              see Twala et al. (2008).}
  \item{nresample}{ number of permutations for \code{testtype = "MonteCarlo"}.}
  \item{tol}{tolerance for zero variances.}
  \item{maxsurrogate}{ number of surrogate splits to evaluate.}
  \item{numsurrogate}{ a logical for backward-compatibility with party. If
                       \code{TRUE}, only at least ordered variables are considered for surrogate splits.}
  \item{mtry}{ number of input variables randomly sampled as candidates
               at each node for random forest like algorithms. The default
               \code{mtry = Inf} means that no random selection takes place.
               If \code{\link{ctree_control}} is used in \code{\link{cforest}}
               this argument is ignored.}
  \item{maxdepth}{ maximum depth of the tree. The default \code{maxdepth = Inf}
                   means that no restrictions are applied to tree sizes.}
  \item{multiway}{ a logical indicating if multiway splits for all factor levels
                   are implemented for unordered factors.}
  \item{splittry}{ number of variables that are inspected for admissible splits
                   if the best split doesn't meet the sample size constraints.}
  \item{intersplit}{ a logical indicating if splits in numeric variables
                     are simply \code{x <= a} (the default) or interpolated
                     \code{x <= (a + b) / 2}. The latter feature is experimental, see
                     Galili and Meilijson (2016).}
  \item{majority}{ if \code{FALSE}, observations which can't be classified to a
                   daughter node because of missing information are randomly
                   assigned (following the node distribution). If \code{TRUE},
                   they go with the majority (the default in \code{\link[party]{ctree}}).}
  \item{caseweights}{ a logical interpreting \code{weights} as case weights.}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
                  \code{function(X, FUN, \dots)}. It is used for computing the variable selection criterion.
                  The default is to use the basic \code{lapply}
                  function unless the \code{cores} argument is specified (see below).
               If \code{\link{ctree_control}} is used in \code{\link{cforest}}
               this argument is ignored.}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
               \code{\link[parallel]{mclapply}} with the desired number of \code{cores}.
               If \code{\link{ctree_control}} is used in \code{\link{cforest}}
               this argument is ignored.}
  \item{saveinfo}{logical. Store information about variable selection
                  procedure in \code{info} slot of each \code{partynode}.}
  \item{update}{logical. If \code{TRUE}, the data transformation is updated
                in every node. The default always was and still is not to
                update unless \code{ytrafo} is a function.}
  \item{splitflavour}{use exhaustive search over splits instead of maximally
                      selected statistics (\code{ctree}). This feature may change.}

}
\details{

  The arguments \code{teststat}, \code{testtype} and \code{mincriterion}
  determine how the global null hypothesis of independence between all input
  variables and the response is tested (see \code{\link{ctree}}).
  The variable with most extreme p-value or test statistic is selected
  for splitting. If this isn't possible due to sample size constraints
  explained in the next paragraph, up to \code{splittry} other variables
  are inspected for possible splits.

  A split is established when all of the following criteria are met:
  1) the sum of the weights in the current node
  is larger than \code{minsplit}, 2) a fraction of the sum of weights of more than
  \code{minprob} will be contained in all daughter nodes, 3) the sum of
  the weights in all daughter nodes exceeds \code{minbucket}, and 4)
  the depth of the tree is smaller than \code{maxdepth}.
  This avoids pathological splits deep down the tree.
  When \code{stump = TRUE}, a tree with at most two terminal nodes is computed.

  The argument \code{mtry > 0} means that a random forest like `variable
  selection', i.e., a random selection of \code{mtry} input variables, is
  performed in each node.

  In each inner node, \code{maxsurrogate} surrogate splits are computed
  (regardless of any missing values in the learning sample). Factors
  in test samples whose levels were empty in the learning sample
  are treated as missing when computing predictions (in contrast
  to \code{\link[party]{ctree}}. Note also the different behaviour of
  \code{majority} in the two implementations.

}
\value{
  A list.
}
\references{

  B. E. T. H. Twala, M. C. Jones, and D. J. Hand (2008),
  Good Methods for Coping with Missing Data in Decision Trees,
  \emph{Pattern Recognition Letters}, \bold{29}(7), 950--956.

  Tal Galili, Isaac Meilijson (2016), Splitting Matters: How
  Monotone Transformation of Predictor Variables May Improve the
  Predictions of Decision Tree Models, \url{https://arxiv.org/abs/1611.04561}.


}
\keyword{misc}
