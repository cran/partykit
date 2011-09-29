\name{nodeapply}
\alias{nodeapply}
\alias{nodeapply.party}
\alias{nodeapply.partynode}
\title{ Apply Functions Over Nodes }
\description{

  Returns a list of values obtained by applying a function to
  \code{party} or \code{partynode} objects.

}
\usage{
nodeapply(obj, ids = 1, FUN = NULL, \dots)
\method{nodeapply}{partynode}(obj, ids = 1, FUN = NULL, \dots)
\method{nodeapply}{party}(obj, ids = 1, FUN = NULL, by_node = TRUE, \dots)
}
\arguments{
  \item{obj}{ an object of class \code{\link{partynode}} 
             or \code{\link{party}}.}
  \item{ids}{ integer vector of node identifiers to apply over.}
  \item{FUN}{ a function to be applied to nodes. By default, the 
              node itself is returned.}
  \item{by_node}{ a logical indicating if \code{FUN} is applied
                  to subsets of \code{\link{party}} objects or 
                  \code{\link{partynode}} objects (default). }
  \item{\dots}{ additional arguments.}
}
\details{

  Function \code{FUN} is applied to all nodes with node identifiers in
  \code{ids} for a \code{partynode} object.
  The method for \code{party} by default calls the \code{nodeapply}
  method on it's \code{node} slot. If \code{by_node} is \code{FALSE},
  it is applied to a \code{party} object with root node \code{ids}.

}
\value{
  A list of results of length \code{length(ids)}.
}
\examples{

    ## a tree as flat list structure
    nodelist <- list(
        # root node
        list(id = 1L, split = partysplit(varid = 4L, breaks = 1.9),
            kids = 2:3),
        # V4 <= 1.9, terminal node
        list(id = 2L, info = "terminal A"),
        # V4 > 1.9
        list(id = 3L, split = partysplit(varid = 5L, breaks = 1.7),
            kids = c(4L, 7L)),
        # V5 <= 1.7
        list(id = 4L, split = partysplit(varid = 4L, breaks = 4.8),
            kids = 5:6),
        # V4 <= 4.8, terminal node
        list(id = 5L, info = "terminal B"),
        # V4 > 4.8, terminal node
        list(id = 6L, info = "terminal C"),
        # V5 > 1.7, terminal node
        list(id = 7L, info = "terminal D")
    )

    ## convert to a recursive structure
    node <- as.partynode(nodelist)

    ## return root node 
    nodeapply(node)

    ## return info slots of terminal nodes
    nodeapply(node, ids = nodeids(node, terminal = TRUE), 
        FUN = function(x) info_node(x))

    ## fit tree using rpart
    library("rpart")
    rp <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)

    ## coerce to `constparty'
    rpk <- as.party(rp)

    ## extract nodeids
    nodeids(rpk)
    unlist(nodeapply(node_party(rpk), ids = nodeids(rpk), 
        FUN = id_node))
    unlist(nodeapply(rpk, ids = nodeids(rpk), FUN = id_node))

    ## but root nodes of party objects always have id = 1
    unlist(nodeapply(rpk, ids = nodeids(rpk), FUN = function(x) 
        id_node(node_party(x)), by_node = FALSE))
}
\keyword{tree}
