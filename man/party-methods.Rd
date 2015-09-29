\name{party-methods}
\alias{party-methods}
\alias{length.party}
\alias{print.party}
\alias{print.simpleparty}
\alias{print.constparty}
\alias{[.party}
\alias{[[.party}
\alias{depth.party}
\alias{width.party}
\alias{getCall.party}
\alias{nodeprune}
\alias{nodeprune.party}
\title{ Methods for Party Objects }
\description{

  Methods for computing on \code{party} objects.

}
\usage{
\method{print}{party}(x, 
    terminal_panel = function(node)
      formatinfo_node(node, default = "*", prefix = ": "),
    tp_args = list(),
    inner_panel = function(node) "", ip_args = list(),
    header_panel = function(party) "",
    footer_panel = function(party) "",
    digits = getOption("digits") - 2, \dots)
\method{print}{simpleparty}(x, digits = getOption("digits") - 4,
    header = NULL, footer = TRUE, \dots)
\method{print}{constparty}(x, FUN = NULL, digits = getOption("digits") - 4,
    header = NULL, footer = TRUE, \dots)
\method{length}{party}(x)
\method{[}{party}(x, i, \dots)
\method{[[}{party}(x, i, \dots)
\method{depth}{party}(x, root = FALSE, \dots)
\method{width}{party}(x, \dots)
\method{nodeprune}{party}(x, ids, ...)
}
\arguments{
  \item{x}{ an object of class \code{\link{party}}.}
  \item{i}{ an integer specifying the root of the subtree to extract.}
  \item{terminal_panel}{ a panel function for printing terminal nodes.}
  \item{tp_args}{ a list containing arguments to \code{terminal_panel}.}
  \item{inner_panel}{ a panel function for printing inner nodes.}
  \item{ip_args}{ a list containing arguments to \code{inner_panel}.}
  \item{header_panel}{ a panel function for printing the header.}
  \item{footer_panel}{ a panel function for printing the footer.}
  \item{digits}{ number of digits to be printed.}
  \item{header}{ header to be printed.}
  \item{footer}{ footer to be printed.}
  \item{FUN}{ a function to be applied to nodes.}
  \item{root}{ a logical. Should the root count be counted in \code{depth}? }
  \item{ids}{ a vector of node ids (or their names) to be pruned-off.}
  \item{\dots}{ additional arguments.}
}
\details{

  \code{length} gives the number of nodes in the tree (in contrast
  to the \code{length} method for \code{\link{partynode}} objects which
  returns the number of kid nodes in the root),
  \code{depth} the depth of the tree and \code{width} 
  the number of terminal nodes. The subset methods extract subtrees
  and the \code{print} method generates a textual representation
  of the tree. \code{nodeprune} prunes-off nodes and makes sure that
  the node ids of the resulting tree are in pre-order starting with root node id 1.
  For \code{constparty} objects, the \code{fitted} slot is also changed.

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
        list(id = 3L, split = partysplit(varid = 5L, breaks = 1.7), 
            kids = c(4L, 7L)),
        # V5 <= 1.7
        list(id = 4L, split = partysplit(varid = 4L, breaks = 4.8), 
            kids = 5:6),
        # V4 <= 4.8, terminal node
        list(id = 5L),
        # V4 > 4.8, terminal node
        list(id = 6L),
        # V5 > 1.7, terminal node
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
    names(tree) <- paste("Node", nodeids(tree), sep = " ")

    ## number of kids in root node
    length(tree)

    ## depth of tree
    depth(tree)

    ## number of terminal nodes
    width(tree)

    ## node number four
    tree["Node 4"]
    tree[["Node 4"]]

}
\keyword{tree}
