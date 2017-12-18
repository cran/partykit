\name{extree_data}
\alias{extree_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Data Preprocessing for Extensible Trees.
}
\description{
  A routine for preprocessing data before an extensible tree can be grown by
  \code{extree_fit}.
}
\usage{
extree_data(formula, data, subset, na.action = na.pass, weights, offset, 
            cluster, strata, scores = NULL, yx = c("none", "matrix"),
            ytype = c("vector", "data.frame", "matrix"), 
            nmax = c(yx = Inf, z = Inf), ...)
}
\arguments{
  \item{formula}{a formula describing the model of the form \code{y1 + y2 + ... ~
x1 + x2 + ... | z1 + z2 + ...}.
}
  \item{data}{an optional data.frame containing the variables in the model.
}
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.
}
  \item{na.action}{a function which indicates what should happen when the
data contain missing values.
}
  \item{weights}{an optional vector of weights.
}
  \item{offset}{an optional offset vector.
}
  \item{cluster}{an optional factor describing clusters. The interpretation
depends on the specific tree algorithm.
}
  \item{strata}{an optional factor describing strata. The interpretation
depends on the specific tree algorithm.
}
  \item{scores}{an optional named list of numeric scores to be assigned to
ordered factors in the \code{z} part of the formula.
}
  \item{yx}{a character indicating if design matrices shall be computed. 
}
  \item{ytype}{a character indicating how response variables shall be
stored.
}
  \item{nmax}{a numeric vector of length two with the maximal number of
bins in the response and \code{x}-part (first element) and the \code{z}
part. Use \code{Inf} to switch-off binning.
}
  \item{\dots}{additional arguments.
}
}
\details{
This internal functionality will be the basis of implementations of other
tree algorithms in future versions. Currently, only \code{ctree} relies on
this function.
}
\value{An object of class \code{extree_data}.
}
\examples{

  data("iris")

  ed <- extree_data(Species ~ Sepal.Width + Sepal.Length | Petal.Width + Petal.Length, 
                    data = iris, nmax = c("yx" = 25, "z" = 10), yx = "matrix")

  ### the model.frame
  mf <- model.frame(ed)
  all.equal(mf, iris[, names(mf)])

  ### binned y ~ x part
  model.frame(ed, yxonly = TRUE)

  ### binned Petal.Width
  ed[[4, type = "index"]]

  ### response
  ed$yx$y

  ### model matrix 
  ed$yx$x

}
\keyword{tree}
