\name{party}
\alias{party}
\alias{names.party}
\alias{names<-.party}
\alias{node_party}
\alias{is.constparty}
\alias{is.simpleparty}
\alias{data_party}
\alias{data_party.default}
\title{ Recursive Partytioning }
\description{

  A class for representing decision trees and corresponding
  accessor functions.

}
\usage{
party(node, data, fitted = NULL, terms = NULL, names = NULL, 
    info = NULL)
\method{names}{party}(x)
\method{names}{party}(x) <- value
data_party(party, id = 1L)
\method{data_party}{default}(party, id = 1L)
node_party(party)
is.constparty(party)
is.simpleparty(party)
}
\arguments{
  \item{node}{ an object of class \code{\link{partynode}}.}
  \item{data}{ a (potentially empty) \code{\link{data.frame}}.}
  \item{fitted}{ an optional \code{\link{data.frame}} with \code{nrow(data)} 
                 rows (only if \code{nrow(data) != 0} and 
                 containing at least the fitted terminal node identifiers
                 as element \code{(fitted)}. In addition, weights
                 may be contained as element \code{(weights)} and 
                 responses as \code{(response)}.}
  \item{terms}{ an optional \code{\link{terms}} object.  }
  \item{names}{ an optional vector of names to be assigned to each node of \code{node}. }
  \item{info}{ additional information. }
  \item{x}{ an object of class \code{party}.}
  \item{party}{ an object of class \code{party}.}
  \item{value}{a character vector of up to the same length as \code{x}, or
    \code{NULL}.}
  \item{id}{ a node identifier.}
}
\details{

  Objects of class \code{party} basically consist of a \code{\link{partynode}}
  object representing the tree structure in a recursive way and 
  data. The \code{data} argument takes a \code{data.frame} which, however,
  might have zero columns. Optionally, a \code{data.frame} with at least one
  variable \code{(fitted)} containing the terminal node numbers of 
  data used for fitting the tree may be specified along with a
  \code{\link{terms}} object or any additional (currently unstructured)
  information as \code{info}. Argument \code{names} defines names
  for all nodes in \code{node}.

  Method \code{names} can be used to extract or alter names for nodes.
  Function \code{node_party} returns the \code{node} element of a
  \code{party} object. Further methods for \code{party} objects
  are documented in \code{\link{party-methods}} and 
  \code{\link{party-predict}}. Trees of various flavors can be coerced
  to \code{party}, see \code{\link{party-coercion}}.
  
  Two classes inherit from class \code{party} and impose additional
  assumptions on the structure of this object:
  Class \code{constparty} requires that the \code{fitted} slot
  contains a partitioning of the learning sample as a factor \code{("fitted")}
  and the response values of all observations in the learning sample
  as \code{("response")}. This structure is most flexible and 
  allows for graphical display of the response values in terminal
  nodes as well as for computing predictions based on 
  arbitrary summary statistics.

  Class \code{simpleparty} assumes that certain pre-computed information
  about the distribution of the response variable is contained
  in the \code{info} slot nodes. At the moment, no formal
  class is used to describe this information.

}
\value{
  The constructor returns an object of class \code{party}:
  \item{node}{ an object of class \code{\link{partynode}}.}
  \item{data}{ a (potentially empty) \code{\link{data.frame}}.}
  \item{fitted}{ an optional \code{\link{data.frame}} with \code{nrow(data)} 
                 rows (only if \code{nrow(data) != 0} and 
                 containing at least the fitted terminal node identifiers
                 as element \code{(fitted)}. In addition, weights
                 may be contained as element \code{(weights)} and 
                 responses as \code{(response)}.}
  \item{terms}{ an optional \code{\link{terms}} object.  }
  \item{names}{ an optional vector of names to be assigned to each node of \code{node}. }
  \item{info}{ additional information. }
  \code{names} can be used to set and retrieve names of nodes and 
  \code{node_party} returns an object of class \code{\link{partynode}}.
  \code{data_party} returns a data frame with observations contained in node
  \code{id}.
}
\references{ 
Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.
}
\examples{
### data ###
## artificial WeatherPlay data
data("WeatherPlay", package = "partykit")
str(WeatherPlay)


### splits ###
## split in overcast, humidity, and windy
sp_o <- partysplit(1L, index = 1:3)
sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)

## query labels
character_split(sp_o)


### nodes ###
## set up partynode structure
pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))
pn


### tree ###
## party: associate recursive partynode structure with data
py <- party(pn, WeatherPlay)
py
plot(py)


### variations ###
## tree stump
n1 <- partynode(id = 1L, split = sp_o, kids = lapply(2L:4L, partynode))
print(n1, data = WeatherPlay)

## query fitted nodes and kids ids
fitted_node(n1, data = WeatherPlay)
kidids_node(n1, data = WeatherPlay)

## tree with full data sets
t1 <- party(n1, data = WeatherPlay)

## tree with empty data set
party(n1, data = WeatherPlay[0, ])

## constant-fit tree
t2 <- party(n1, 
  data = WeatherPlay,
  fitted = data.frame(
    "(fitted)" = fitted_node(n1, data = WeatherPlay),
    "(response)" = WeatherPlay$play,
    check.names = FALSE),
  terms = terms(play ~ ., data = WeatherPlay),
)
t2 <- as.constparty(t2)
t2
plot(t2)
}
\keyword{tree}
