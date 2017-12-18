\name{nodeids}
\alias{nodeids}
\alias{nodeids.party}
\alias{nodeids.partynode}
\alias{get_paths}
\title{ Extract Node Identifiers }
\description{

  Extract unique identifiers from inner and terminals
  nodes of a \code{partynode} object.

}
\usage{
nodeids(obj, \dots)
\method{nodeids}{partynode}(obj, from = NULL, terminal = FALSE, \dots)
\method{nodeids}{party}(obj, from = NULL, terminal = FALSE, \dots)
get_paths(obj, i)
}
\arguments{
  \item{obj}{ an object of class \code{\link{partynode}} 
             or \code{\link{party}}.}
  \item{from}{ an integer specifying node to start from.}
  \item{terminal}{ logical specifying if only node identifiers 
                   of terminal nodes are returned. }
  \item{i}{a vector of node identifiers.}
  \item{\dots}{ additional arguments.}
}
\details{
  The identifiers of each node are extracted from \code{nodeids}.
  \code{get_paths} returns the paths for extracting the corresponding
  nodes using list subsets.
}
\value{
  A vector of node identifiers.
}
\examples{

  ## a tree as flat list structure
  nodelist <- list(
      # root node
      list(id = 1L, split = partysplit(varid = 4L, breaks = 1.9),
          kids = 2:3),
      # V4 <= 1.9, terminal node
      list(id = 2L),
      # V4 > 1.9
      list(id = 3L, split = partysplit(varid = 1L, breaks = 1.7),
          kids = c(4L, 7L)),
      # V1 <= 1.7
      list(id = 4L, split = partysplit(varid = 4L, breaks = 4.8),
          kids = 5:6),
      # V4 <= 4.8, terminal node
      list(id = 5L),
      # V4 > 4.8, terminal node
      list(id = 6L),
      # V1 > 1.7, terminal node
      list(id = 7L)
  )

  ## convert to a recursive structure
  node <- as.partynode(nodelist)

  ## set up party object
  data("iris")
  tree <- party(node, data = iris, 
      fitted = data.frame("(fitted)" = 
                          fitted_node(node, data = iris),
                          check.names = FALSE))
  tree
  
  ### ids of all nodes
  nodeids(tree)

  ### ids of all terminal nodes
  nodeids(tree, terminal = TRUE)

  ### ids of terminal nodes in subtree with root [3]
  nodeids(tree, from = 3, terminal = TRUE)

  ### get paths and extract all terminal nodes
  tr <- unclass(node_party(tree))
  lapply(get_paths(tree, nodeids(tree, terminal = TRUE)), 
         function(path) tr[path])

}
\keyword{tree}
