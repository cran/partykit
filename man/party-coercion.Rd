\name{party-coercion}
\alias{party-coercion}
\alias{as.party}
\alias{as.party.rpart}
\alias{as.party.Weka_tree}
\alias{as.party.XMLNode}
\alias{as.constparty}
\alias{as.simpleparty}
\alias{as.simpleparty.party}
\alias{as.simpleparty.simpleparty}
\alias{as.simpleparty.XMLNode}
\alias{as.simpleparty.constparty}
\alias{pmmlTreeModel}

\title{Coercion Functions}

\description{
Functions coercing various objects to objects of class party.
}

\usage{
as.party(obj, \dots)
\method{as.party}{rpart}(obj, data = TRUE, \dots)
\method{as.party}{Weka_tree}(obj, data = TRUE, \dots)
\method{as.party}{XMLNode}(obj, \dots)
pmmlTreeModel(file, \dots)
as.constparty(obj, \dots)
as.simpleparty(obj, \dots)
\method{as.simpleparty}{party}(obj, \dots)
\method{as.simpleparty}{simpleparty}(obj, \dots)
\method{as.simpleparty}{constparty}(obj, \dots)
\method{as.simpleparty}{XMLNode}(obj, \dots)
}

\arguments{
  \item{obj}{ an object of class \code{\link[rpart]{rpart}},
              \code{\link[RWeka:Weka_classifier_trees]{Weka_tree}}, \code{XMLnode} or
              objects inheriting from \code{party}.}
  \item{data}{logical. Should the model frame associated with the fitted
              \code{obj} be included in the \code{data} of the \code{party}?}
  \item{file}{ a file name of a XML file 
              containing a PMML description of a tree.}
  \item{\dots}{ additional arguments.}
}

\details{
Trees fitted using functions \code{\link[rpart]{rpart}} or 
\code{\link[RWeka:Weka_classifier_trees]{J48}} are coerced to \code{\link{party}}
objects. By default, objects of class \code{constparty}
are returned.

When information about the learning sample is available, 
\code{\link{party}} objects can be coerced to objects of 
class \code{constparty} or \code{simpleparty} 
(see \code{\link{party}} for details).
}

\value{
  All methods return objects of class \code{\link{party}}.
}

\examples{
## fit tree using rpart
library("rpart")
rp <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)

## coerce to `constparty'
as.party(rp)
}

\keyword{tree}
