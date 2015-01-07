\name{mob_control}
\alias{mob_control}

\title{Control Parameters for Model-Based Partitioning}

\description{
  Various parameters that control aspects the fitting algorithm
  for recursively partitioned \code{\link{mob}} models.
}

\usage{
mob_control(alpha = 0.05, bonferroni = TRUE, minsize = NULL, maxdepth = Inf,
  mtry = Inf, trim = 0.1, breakties = FALSE, parm = NULL, dfsplit = TRUE, prune = NULL,
  restart = TRUE, verbose = FALSE, caseweights = TRUE, ytype = "vector", xtype = "matrix",
  terminal = "object", inner = terminal, model = TRUE, numsplit = "left",
  catsplit = "binary", vcov = "opg", ordinal = "chisq", nrep = 10000,
  minsplit = minsize, minbucket = minsize, applyfun = NULL, cores = NULL)
}

\arguments{
  \item{alpha}{numeric significance level. A node is splitted when
    the (possibly Bonferroni-corrected) \eqn{p} value for any parameter
    stability test in that node falls below \code{alpha} (and the stopping
    criteria \code{minsize} and \code{maxdepth} are not fulfilled).}
  \item{bonferroni}{logical. Should \eqn{p} values be Bonferroni corrected?}
  \item{minsize, minsplit, minbucket}{integer. The minimum number of observations in a node.
    If \code{NULL}, the default is to use 10 times the number of parameters
    to be estimated (divided by the number of responses per observation
    if that is greater than 1). \code{minsize} is the recommended name and
    \code{minsplit}/\code{minbucket} are only included for backward compatibility with previous
    versions of \code{mob} and compatibility with \code{ctree}, respectively.}
  \item{maxdepth}{integer. The maximum depth of the tree.}
  \item{mtry}{integer. The number of partitioning variables randomly sampled
    as candidates in each node for forest-style algorithms. If \code{mtry}
    is greater than the number of partitioning variables, no random selection
    is performed. (Thus, by default all available partitioning variables are considered.)}
  \item{trim}{numeric. This specifies the trimming in the parameter instability
    test for the numerical variables. If smaller than 1, it is interpreted
    as the fraction relative to the current node size.}
  \item{breakties}{logical. Should ties in numeric variables be broken
    randomly for computing the associated parameter instability test?}
  \item{parm}{numeric or character. Number or name of model parameters
    included in the parameter instability tests (by default all parameters
    are included).}
  \item{dfsplit}{logical or numeric. \code{as.integer(dfsplit)} is the degrees of freedom
    per selected split employed when computing information criteria etc.}
  \item{prune}{character, numeric, or function for specifying post-pruning rule.
    If \code{prune} is \code{NULL} (the default), no post-pruning is performed.
    For likelihood-based \code{mob()} trees, \code{prune} can be set to
    \code{"AIC"} or \code{"BIC"} for post-pruning based on the corresponding
    information criteria. More general rules (also in scenarios that are
    not likelihood-based), can be specified by function arguments to
    \code{prune}, for details see below.}
  \item{restart}{logical. When determining the optimal split point in a numerical
    variable: Should model estimation be restarted with \code{NULL} starting
    values for each split? The default is \code{TRUE}. If \code{FALSE}, then
    the parameter estimates from the previous split point are used as starting
    values for the next split point (because in practice the difference are
    often not huge). (Note that in that case a \code{for} loop is used
    instead of the \code{applyfun} for fitting models across sample splits.)}
  \item{verbose}{logical. Should information about the fitting process
    of \code{\link{mob}} (such as test statistics, \eqn{p} values, selected
    splitting variables and split points) be printed to the screen?}
  \item{caseweights}{logical. Should weights be interpreted as case weights?
    If \code{TRUE}, the number of observations is \code{sum(weights)},
    otherwise it is \code{sum(weights > 0)}.}
  \item{ytype, xtype}{character. Specification of how \code{mob} should
    preprocess \code{y} and \code{x} variables. Possible choice are:
    \code{"vector"} (for \code{y} only), i.e., only one variable;
    \code{"matrix"}, i.e., the model matrix of all variables;
    \code{"data.frame"}, i.e., a data frame of all variables.}
  \item{terminal, inner}{character. Specification of which additional
    information (\code{"estfun"}, \code{"object"}, or both) should be
    stored in each node. If \code{NULL}, no additional information is
    stored.}
  \item{model}{logical. Should the full model frame be stored in the
    resulting object?}
  \item{numsplit}{character indicating how splits for numeric variables
    should be justified. Because any splitpoint in the interval between
    the last observation from the left child segment and the first observation
    from the right child segment leads to the same observed split, two
    options are available in \code{mob_control}:
    Either, the split is \code{"left"}-justified (the default for
    backward compatibility) or \code{"center"}-justified using the
    midpoint of the possible interval.}
  \item{catsplit}{character indicating how (unordered) categorical variables
    should be splitted. By default the best \code{"binary"} split is
    searched (by minimizing the objective function). Alternatively,
    if set to \code{"multiway"}, the node is simply splitted into all
    levels of the categorical variable.}
  \item{vcov}{character indicating which type of covariance matrix
    estimator should be employed in the parameter instability tests.
    The default is the outer product of gradients (\code{"opg"}).
    Alternatively, \code{vcov = "info"} employs the information matrix
    and \code{vcov = "sandwich"} the sandwich matrix (both of which are
    only sensible for maximum likelihood estimation).}
  \item{ordinal}{character indicating which type of parameter instability
    test should be employed for ordinal partitioning variables (i.e.,
    ordered factors). This can be \code{"chisq"}, \code{"max"}, or \code{"L2"}.
    If \code{"chisq"} then the variable is treated as unordered and a
    chi-squared test is performed. If \code{"L2"}, then a maxLM-type
    test as for numeric variables is carried out but correcting for ties.
    This requires simulation of p-values via \code{\link[strucchange]{catL2BB}}
    and requires some computation time. For \code{"max"} a weighted
    double maximum test is used that computes p-values via
    \code{\link[mvtnorm]{pmvnorm}}.}
  \item{nrep}{numeric. Number of replications in the simulation of
    p-values for the ordinal \code{"L2"} statistic (if used).}
  \item{applyfun}{an optional \code{\link[base]{lapply}}-style function with arguments
    \code{function(X, FUN, \dots)}. It is used for refitting the model across
    potential sample splits. The default is to use the basic \code{lapply}
    function unless the \code{cores} argument is specified (see below).}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
    \code{\link[parallel]{mclapply}} with the desired number of \code{cores}.}
}

\details{
  See \code{\link{mob}} for more details and references.
  
  For post-pruning, \code{prune} can be set to a \code{function(objfun, df, nobs)}
  which either returns \code{TRUE} to signal that a current node can be pruned
  or \code{FALSE}. All supplied arguments are of length two: \code{objfun} is the sum of objective
  function values in the current node and its child nodes, respectively.
  \code{df} is the degrees of freedom in the current node and its child nodes,
  respectively. \code{nobs} is vector with the number of observations in the
  current node and the total number of observations in the dataset, respectively.  

  If the objective function employed in the \code{mob()} call is the negative
  log-likelihood, then a suitable function is set up on the fly by comparing  
  \code{(2 * objfun + penalty * df)} in the current and the daughter nodes.
  The penalty can then be set via a numeric or character value for \code{prune}:
  AIC is used if \code{prune = "AIC"} or \code{prune = 2} and
  BIC if \code{prune = "BIC"} or \code{prune = log(n)}.
}

\seealso{\code{\link{mob}}}

\value{
  A list of class \code{mob_control} containing the control parameters.
}

\keyword{misc}
