\name{ctree}
\alias{ctree}
\title{ Conditional Inference Trees }
\description{
  Recursive partitioning for continuous, censored, ordered, nominal and
  multivariate response variables in a conditional inference framework. 
}
\usage{
ctree(formula, data, weights, subset, na.action = na.pass,
    control = ctree_control(\dots), \dots)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ a data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Only non-negative integer valued weights are
                  allowed.}
  \item{na.action}{a function which indicates what should happen when the data
                   contain missing value.}
  \item{control}{a list with control parameters, see
                 \code{\link{ctree_control}}.}
  \item{\dots}{arguments passed to \code{\link{ctree_control}}.}
}
\details{

  Function \code{partykit::ctree} is a reimplementation of (most of)
  \code{party::ctree} employing the new \code{\link{party}} infrastructure
  of the \pkg{partykit} infrastructure. Although the new code was already 
  extensively tested, it is not yet as mature as the old code. If you notice
  differences in the structure/predictions of the resulting trees, please
  contact the package maintainers. See also below for some remarks about
  the internals of the different implementations and how they should be
  merged in future releases.
  
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
  the right sized tree is grown and no form of pruning or cross-validation
  or whatsoever is needed. The selection of the input variable to split in
  is based on the univariate p-values avoiding a variable selection bias
  towards input variables with many possible cutpoints.

  Predictions can be computed using \code{\link{predict}}, which returns predicted means,
  predicted classes or median predicted survival times and 
  more information about the conditional
  distribution of the response, i.e., class probabilities
  or predicted Kaplan-Meier curves. For observations
  with zero weights, predictions are computed from the fitted tree 
  when \code{newdata = NULL}.

  For a general description of the methodology see Hothorn, Hornik and
  Zeileis (2006) and Hothorn, Hornik, van de Wiel and Zeileis (2006).

  Implementation details and roadmap: As pointed above, the function \code{ctree}
  is a reimplementation of \code{\link[party]{ctree}}. Not only
  the R code changed but also the underlying C code which at the moment
  does not support the \code{xtrafo} and \code{ytrafo} arguments due to
  efficiency considerations. The roadmap for future releases is the following:
  (1) Make \pkg{party} depend on \pkg{partykit}. (2) Merge the R code into
  a single \code{ctree} function but keeping the two underlying C functions
  separate. (3) The new interface will always return an object of class
  \code{\link{party}} but call the new and more efficient C code only if
  \code{xtrafo} was not used (while \code{ytrafo} should be integrated into
  the new C code).

}
\value{
  An object of class \code{\link{party}}.
}
\references{ 

   Helmut Strasser and Christian Weber (1999). On the asymptotic theory of permutation
   statistics. \emph{Mathematical Methods of Statistics}, \bold{8}, 220--250.

   Torsten Hothorn, Kurt Hornik, Mark A. van de Wiel and Achim Zeileis (2006).
   A Lego System for Conditional Inference. \emph{The American Statistician},
   \bold{60}(3), 257--263.

   Torsten Hothorn, Kurt Hornik and Achim Zeileis (2006). Unbiased Recursive
   Partitioning: A Conditional Inference Framework. \emph{Journal of
   Computational and Graphical Statistics}, \bold{15}(3), 651--674. 
   Preprint available from 
   \url{http://eeecon.uibk.ac.at/~zeileis/papers/Hothorn+Hornik+Zeileis-2006.pdf}

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
    if (require("ipred")) {
        data("GBSG2", package = "ipred")
        GBSG2ct <- ctree(Surv(time, cens) ~ .,data = GBSG2)
        predict(GBSG2ct, newdata = GBSG2[1:2,], type = "response")        
    }

    ### multivariate responses
    airq2 <- ctree(Ozone + Temp ~ ., data = airq)
    airq2
    plot(airq2)

}
\keyword{tree}
