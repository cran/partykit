\name{extree_fit}
\alias{extree_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit Extensible Trees.
}
\description{
  Basic infrastructure for fitting extensible trees.
}
\usage{
extree_fit(data, trafo, converged, selectfun = ctrl$selectfun, splitfun = ctrl$splitfun, 
           svselectfun = ctrl$svselectfun, svsplitfun = ctrl$svsplitfun, partyvars, 
           subset, weights, ctrl, doFit = TRUE)
}
\arguments{
  \item{data}{an object of class \code{extree_data}, see
\code{\link{extree_data}}.
}
  \item{trafo}{a function with arguments \code{subset}, \code{weights},
\code{info}, \code{estfun} and \code{object}.
}
  \item{converged}{a function with arguments \code{subset}, \code{weights}.
}
  \item{selectfun}{an optional function for selecting variables.
}
  \item{splitfun}{an optional function for selecting splits.
}
  \item{svselectfun}{an optional function for selecting surrogate variables.
}
  \item{svsplitfun}{an optional function for selecting surrogate splits.
}
  \item{partyvars}{a numeric vector assigning a weight to each partitioning
variable (\code{z} in \code{\link{extree_data}}.
}
  \item{subset}{a sorted integer vector describing a subset.
}
  \item{weights}{an optional vector of weights.
}
  \item{ctrl}{control arguments.
}
  \item{doFit}{a logical indicating if the tree shall be grown (\code{TRUE})
or not \code{FALSE}.
}
}
\details{
This internal functionality will be the basis of implementations of other
tree algorithms in future versions. Currently, only \code{ctree} relies on
this function.
}
\value{An object of class \code{partynode}.}
\keyword{tree}
