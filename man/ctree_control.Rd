\name{ctree_control}
\alias{ctree_control}
\title{ Control for Conditional Inference Trees }
\description{

  Various parameters that control aspects of the `ctree' fit.

}
\usage{
ctree_control(teststat = c("quad", "max"),
    testtype = c("Bonferroni", "Univariate", "Teststatistic"),
    mincriterion = 0.95, minsplit = 20L, minbucket = 7L, 
    minprob = 0.01, stump = FALSE, maxsurrogate = 0L, mtry = Inf, 
    maxdepth = Inf, multiway = FALSE, splittry = 2L, majority = FALSE)
}
\arguments{
  \item{teststat}{ a character specifying the type of the test statistic
                       to be applied. }
  \item{testtype}{ a character specifying how to compute the distribution of
                   the test statistic. }
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to implement a split. }
  \item{minsplit}{ the minimum sum of weights in a node in order to be considered
                   for splitting. }
  \item{minbucket}{ the minimum sum of weights in a terminal node. }
  \item{minprob}{ proportion of observations needed to establish a terminal node.}
  \item{stump}{ a logical determining whether a stump (a tree with three
                nodes only) is to be computed. }
  \item{maxsurrogate}{ number of surrogate splits to evaluate. Note the
                       currently only surrogate splits in ordered
                       covariables are implemented. }
  \item{mtry}{ number of input variables randomly sampled as candidates
               at each node for random forest like algorithms. The default
               \code{mtry = Inf} means that no random selection takes place.}
  \item{maxdepth}{ maximum depth of the tree. The default \code{maxdepth = Inf}
                   means that no restrictions are applied to tree sizes.}
  \item{multiway}{ a logical indicating if multiway splits for all factor levels
                   are implemented for unordered factors.}
  \item{splittry}{ number of variables that are inspected for admissible splits
                   if the best split doesn't meet the sample size constraints.}
  \item{majority}{ if \code{FALSE}, observations which can't be classified to a
                   daughter node because of missing information are randomly
                   assigned (following the node distribution). If \code{FALSE},
                   they go with the majority (the default in \code{\link[party]{ctree}}).}
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
\keyword{misc}
