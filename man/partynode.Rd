\name{partynode}
\alias{partynode}
\alias{kidids_node}
\alias{fitted_node}
\alias{id_node}
\alias{split_node}
\alias{surrogates_node}
\alias{kids_node}
\alias{info_node}
\alias{formatinfo_node}
\title{ Inner and Terminal Nodes }
\description{

  A class for representing inner and terminal 
  nodes in trees and functions for data partitioning.

}
\usage{
partynode(id, split = NULL, kids = NULL, surrogates = NULL, 
    info = NULL)
kidids_node(node, data, vmatch = 1:ncol(data), 
    obs = NULL, perm = NULL)
fitted_node(node, data, vmatch = 1:ncol(data), 
    obs = 1:nrow(data), perm = NULL)
id_node(node)
split_node(node)
surrogates_node(node)
kids_node(node)
info_node(node)
formatinfo_node(node, FUN = NULL, default = "", prefix = NULL, \dots)
}
\arguments{
  \item{id}{ integer, a unique identifier for a node. }
  \item{split}{ an object of class \code{\link{partysplit}}. }
  \item{kids}{ a list of \code{partynode} objects. }
  \item{surrogates}{ a list of \code{partysplit} objects.}
  \item{info}{ additional information. }
  \item{node}{ an object of class \code{partynode}.}
  \item{data}{ a \code{\link{list}} or \code{\link{data.frame}}.}
  \item{vmatch}{ a permutation of the variable numbers in \code{data}.}
  \item{obs}{ a logical or integer vector indicating a subset of the
              observations in \code{data}.}
  \item{perm}{ a vector of integers specifying the variables
               to be permuted prior before splitting (i.e., for
               computing permutation variable importances). The
               default \code{NULL} doesn't alter the data.}
  \item{FUN}{ function for formatting the \code{info}, for default see below.}
  \item{default}{ a character used if the \code{info} in \code{node}
                  is \code{NULL}.}
  \item{prefix}{ an optional prefix to be added to the returned character. }
  \item{\dots}{ further arguments passed to
                \code{\link[utils]{capture.output}}.}
}
\details{

  A node represents both inner and terminal nodes in a tree structure. 
  Each node has a unique identifier \code{id}. A node consisting only
  of such an identifier (and possibly additional information
  in \code{info}) is a terminal node.

  Inner nodes consist of a primary split 
  (an object of class \code{\link{partysplit}})
  and at least two kids (daughter nodes). 
  Kid nodes are objects of class \code{partynode}
  itself, so the tree structure is defined recursively.
  In addition, a list of \code{partysplit} objects 
  offering surrogate splits can be supplied. Like
  \code{\link{partysplit}} objects, \code{partynode}
  objects aren't connected to the actual data.

  Function \code{kidids_node()} determines how
  the observations in \code{data[obs,]} are partitioned
  into the kid nodes and returns the number of the list element
  in list \code{kids} each observations belongs to 
  (and not it's identifier).
  This is done by evaluating \code{split} (and possibly
  all surrogate splits) on \code{data} using 
  \code{\link{kidids_split}}.

  Function \code{fitted_node()} performs all 
  splits recursively and returns the identifier \code{id}
  of the terminal node each observation in
  \code{data[obs,]} belongs to. Arguments \code{vmatch}, 
  \code{obs} and \code{perm} are passed to \code{\link{kidids_split}}.

  Function \code{formatinfo_node()} extracts the the \code{info}
  from \code{node} and formats it to a \code{character} vector using the following
  strategy: If \code{is.null(info)}, the \code{default} is returned.  
  Otherwise, \code{FUN} is applied for formatting. The default function uses
  \code{as.character} for atomic objects and applies \code{\link[utils]{capture.output}}
  to \code{print(info)} for other objects. Optionally, a \code{prefix} can be added
  to the computed character string.

  All other functions are accessor functions for
  extracting information from objects of class \code{partynode}.

}
\value{
  The constructor \code{partynode()} returns an object of class \code{partynode}:
  \item{id}{ a unique integer identifier for a node. }
  \item{split}{ an object of class \code{\link{partysplit}}. }
  \item{kids}{ a list of \code{partynode} objects. }
  \item{surrogates}{ a list of \code{\link{partysplit}} objects.}
  \item{info}{ additional information. }

  \code{kidids_split()} returns an integer vector describing
  the partition of the observations into kid nodes by their position
  in list \code{kids}.

  \code{fitted_node()} returns the node identifiers (\code{id}) 
  of the terminal nodes each observation belongs to.

}
\references{ 
Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.
}
\examples{
data("iris", package = "datasets")

## a stump defined by a binary split in Sepal.Length
stump <- partynode(id = 1L, 
    split = partysplit(which(names(iris) == "Sepal.Length"),
	breaks = 5),
    kids = lapply(2:3, partynode))

## textual representation
print(stump, data = iris)  

## list element number and node id of the two terminal nodes
table(kidids_node(stump, iris), 
    fitted_node(stump, data = iris))

## assign terminal nodes with probability 0.5
## to observations with missing `Sepal.Length'
iris_NA <- iris
iris_NA[sample(1:nrow(iris), 50), "Sepal.Length"] <- NA
table(fitted_node(stump, data = iris_NA, 
    obs = !complete.cases(iris_NA)))

## a stump defined by a primary split in `Sepal.Length'
## and a surrogate split in `Sepal.Width' which
## determines terminal nodes for observations with
## missing `Sepal.Length'
stump <- partynode(id = 1L, 
    split = partysplit(which(names(iris) == "Sepal.Length"),
	breaks = 5),
    kids = lapply(2:3, partynode),
    surrogates = list(partysplit(
	which(names(iris) == "Sepal.Width"), breaks = 3)))
f <- fitted_node(stump, data = iris_NA, 
    obs = !complete.cases(iris_NA))
tapply(iris_NA$Sepal.Width[!complete.cases(iris_NA)], f, range)

}
\keyword{tree}
