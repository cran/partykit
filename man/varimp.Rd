\name{varimp}
\alias{varimp}
\alias{varimp.constparty}
\alias{varimp.cforest}
\title{ Variable Importance }
\description{
    Standard and conditional variable importance for `cforest', following the permutation
    principle of the `mean decrease in accuracy' importance in `randomForest'.
}
\usage{
\method{varimp}{constparty}(object, nperm = 1L, 
       risk = c("loglik", "misclassification"), conditions = NULL,
       mincriterion = 0, ...)
\method{varimp}{cforest}(object, nperm = 1L, 
       OOB = TRUE, risk = c("loglik", "misclassification"),
       conditional = FALSE, threshold = .2, ...)
}
\arguments{
  \item{object}{ an object as returned by \code{cforest}.}
  \item{mincriterion}{ the value of the test statistic or 1 - p-value that
                       must be exceeded in order to include a split in the 
                       computation of the importance. The default \code{mincriterion = 0}
                       guarantees that all splits are included.}
  \item{conditional}{ a logical determining whether unconditional or conditional 
                      computation of the importance is performed. }
  \item{threshold}{ the value of the test statistic or 1 - p-value of the association 
                    between the variable of interest and a covariate that must be 
                    exceeded inorder to include the covariate in the conditioning 
                    scheme for the variable of interest (only relevant if 
                    \code{conditional = TRUE}). }
  \item{nperm}{ the number of permutations performed.}
  \item{OOB}{ a logical determining whether the importance is computed from the out-of-bag 
              sample or the learning sample (not suggested).}
  \item{risk}{ a character determining the risk to be evaluated.}
  \item{conditions}{ a list of conditions. }
  \item{\dots}{additional arguments, not used.}
}
\details{

  NEEDS UPDATE

  Function \code{varimp} can be used to compute variable importance measures
  similar to those computed by \code{\link[randomForest]{importance}}. Besides the
  standard version, a conditional version is available, that adjusts for correlations between
  predictor variables. 
  
  If \code{conditional = TRUE}, the importance of each variable is computed by permuting 
  within a grid defined by the covariates that are associated  (with 1 - p-value 
  greater than \code{threshold}) to the variable of interest.
  The resulting variable importance score is conditional in the sense of beta coefficients in   
  regression models, but represents the effect of a variable in both main effects and interactions.
  See Strobl et al. (2008) for details.

  Note, however, that all random forest results are subject to random variation. Thus, before
  interpreting the importance ranking, check whether the same ranking is achieved with a
  different random seed -- or otherwise increase the number of trees \code{ntree} in 
  \code{\link{ctree_control}}.

  Note that in the presence of missings in the predictor variables the procedure
  described in Hapfelmeier et al. (2012) is performed.

}
\value{
  A vector of `mean decrease in accuracy' importance scores.
}
\references{ 

    Leo Breiman (2001). Random Forests. \emph{Machine Learning}, 45(1), 5--32.

    Alexander Hapfelmeier, Torsten Hothorn, Kurt Ulm, and Carolin Strobl (2012).
    A New Variable Importance Measure for Random Forests with Missing Data.
    \emph{Statistics and Computing}, \url{http://dx.doi.org/10.1007/s11222-012-9349-1}

    Torsten Hothorn, Kurt Hornik, and Achim Zeileis (2006b). Unbiased
    Recursive Partitioning: A Conditional Inference Framework.
    \emph{Journal of Computational and Graphical Statistics}, \bold{15} (3),
    651-674.  Preprint available from 
    \url{http://statmath.wu-wien.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf}

    Carolin Strobl, Anne-Laure Boulesteix, Thomas Kneib, Thomas Augustin, and Achim Zeileis (2008).
    Conditional Variable Importance for Random Forests. \emph{BMC Bioinformatics}, \bold{9}, 307. 
    \url{http://www.biomedcentral.com/1471-2105/9/307}
}
\examples{
    
   set.seed(290875)
   data("readingSkills", package = "party")
   readingSkills.cf <- cforest(score ~ ., data = readingSkills, 
                               mtry = 2, ntree = 50)

   # standard importance
   varimp(readingSkills.cf)

   # conditional importance, may take a while...
   varimp(readingSkills.cf, conditional = TRUE)

}
\keyword{tree}
