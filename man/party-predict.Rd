\name{party-predict}
\alias{party-predict}
\alias{predict.party}
\alias{predict_party}
\alias{predict_party.default}
\alias{predict_party.constparty}
\alias{predict_party.simpleparty}
\title{ Tree Predictions }
\description{

  Compute predictions from \code{party} objects.

}
\usage{
\method{predict}{party}(object, newdata = NULL, perm = NULL, \dots)
predict_party(party, id, newdata = NULL, \dots)
\method{predict_party}{default}(party, id, newdata = NULL, FUN = NULL, \dots) 
\method{predict_party}{constparty}(party, id, newdata = NULL,
    type = c("response", "prob", "quantile", "density", "node"),
    at = if (type == "quantile") c(0.1, 0.5, 0.9),
    FUN = NULL, simplify = TRUE, \dots)
\method{predict_party}{simpleparty}(party, id, newdata = NULL,
    type = c("response", "prob", "node"), \dots)
}
\arguments{
  \item{object}{ objects of class \code{\link{party}}. }
  \item{newdata}{ an optional data frame in which to look for variables with
          which to predict, if omitted, the fitted values are used.}
  \item{perm}{an optional character vector of variable names. Splits of
              nodes with a primary split in any of these variables will
              be permuted (after dealing with surrogates). Note that
              surrogate split in the \code{perm} variables will no be permuted.}
  \item{party}{ objects of class \code{\link{party}}. }
  \item{id}{ a vector of terminal node identifiers. }
  \item{type}{ a character string denoting the type of predicted value
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
          \code{"node"} returns an integer vector of terminal node
          identifiers.}
  \item{FUN}{ a function to extract (\code{default} method) or compute 
              (\code{constparty} method) summary statistics. For the \code{default} method,
              this is a function of a terminal node only, for the \code{constparty} method, 
              predictions for each node have to be computed based on arguments \code{(y, w)}
              where \code{y} is the response and \code{w} are case weights.}
  \item{at}{ if the return value is a function (as the empirical cumulative distribution
             function or the empirical quantile function), this function is evaluated
             at values \code{at} and these numeric values are returned. If \code{at} is
             \code{NULL}, the functions themselves are returned in a list.}
  \item{simplify}{ a logical indicating whether the resulting list 
                   of predictions should be converted to a suitable
                   vector or matrix (if possible).}
  \item{\dots}{ additional arguments. }
}
\details{

  The \code{\link{predict}} method for \code{\link{party}} objects
  computes the identifiers of the predicted terminal nodes, either
  for new data in \code{newdata} or for the learning samples
  (only possible for objects of class \code{constparty}).
  These identifiers are delegated to the corresponding
  \code{predict_party} method which computes (via 
  \code{FUN} for class \code{constparty})
  or extracts (class \code{simpleparty}) the actual predictions.

}
\value{
  A list of predictions, possibly simplified to a numeric vector, 
  numeric matrix or factor.
}
\examples{

  ## fit tree using rpart
  library("rpart")
  rp <- rpart(skips ~ Opening + Solder + Mask + PadType + Panel,
              data = solder, method = 'anova')

  ## coerce to `constparty'
  pr <- as.party(rp)

  ## mean predictions
  predict(pr, newdata = solder[c(3, 541, 640),])

  ## ecdf
  predict(pr, newdata = solder[c(3, 541, 640),], type = "prob")

  ## terminal node identifiers
  predict(pr, newdata = solder[c(3, 541, 640),], type = "node")

  ## median predictions
  predict(pr, newdata = solder[c(3, 541, 640),], 
          FUN = function(y, w = 1) median(y))

}
\keyword{tree}
