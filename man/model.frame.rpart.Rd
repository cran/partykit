\name{model_frame_rpart}
\alias{model_frame_rpart}
\title{ Model Frame Method for rpart }
\description{

  A model.frame method for rpart objects.

}
\usage{
model_frame_rpart(formula, \dots)
}
\arguments{
  \item{formula}{ an object of class \code{\link{rpart}}.}
  \item{\dots}{ additional arguments.}
}
\details{
  A \code{\link{model.frame}} method for \code{\link{rpart}} objects.
  Because it is no longer possible to overwrite existing methods,
  the function name is a little different here.
}
\value{
  A model frame.
}
\keyword{tree}
