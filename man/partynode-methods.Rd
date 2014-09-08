\name{partynode-methods}
\alias{partynode-methods}
\alias{is.partynode}
\alias{as.partynode}
\alias{as.partynode.partynode}
\alias{as.partynode.list}
\alias{as.list.partynode}
\alias{length.partynode}
\alias{[.partynode}
\alias{[[.partynode}
\alias{is.terminal}
\alias{is.terminal.partynode}
\alias{depth.partynode}
\alias{width}
\alias{width.partynode}
\alias{print.partynode}
\alias{nodeprune.partynode}
\title{ Methods for Node Objects}
\description{

  Methods for computing on \code{partynode} objects.

}
\usage{
is.partynode(x)
as.partynode(x, \dots)
\method{as.partynode}{partynode}(x, from = NULL, recursive = TRUE, \dots)
\method{as.partynode}{list}(x, \dots)
\method{as.list}{partynode}(x, \dots)
\method{length}{partynode}(x)
\method{[}{partynode}(x, i, \dots)
\method{[[}{partynode}(x, i, \dots)
is.terminal(x, \dots)
\method{is.terminal}{partynode}(x, \dots)
\method{depth}{partynode}(x, root = FALSE, \dots)
width(x, \dots)
\method{width}{partynode}(x, \dots)
\method{print}{partynode}(x, data = NULL, names = NULL,
    inner_panel = function(node) "", 
    terminal_panel = function(node) " *",
    prefix = "", first = TRUE, digits = getOption("digits") - 2, 
    \dots)
\method{nodeprune}{partynode}(x, ids, ...)
}
\arguments{
  \item{x}{ an object of class \code{partynode} or \code{list}.}
  \item{from}{ an integer giving the identifier of the root node.}
  \item{recursive}{ a logical, if \code{FALSE}, only the id of the root node
                    is checked against \code{from}. If 
                    \code{TRUE}, the ids of all nodes are checked.}
  \item{i}{ an integer specifying the kid to extract.}
  \item{root}{ a logical. Should the root count be counted in \code{depth}? }
  \item{data}{ an optional \code{data.frame}.}
  \item{names}{ a vector of names for nodes.}
  \item{terminal_panel}{ a panel function for printing terminal nodes.}
  \item{inner_panel}{ a panel function for printing inner nodes.}
  \item{prefix}{ lines start with this symbol.}
  \item{first}{ a logical.}
  \item{digits}{ number of digits to be printed.}
  \item{ids}{ a vector of node ids to be pruned-off.}
  \item{\dots}{ additional arguments.}
}
\details{

  \code{is.partynode} checks if the argument is a valid \code{partynode}
  object. \code{is.terminal} is \code{TRUE} for terminal nodes
  and \code{FALSE} for inner nodes. The subset methods
  return the \code{partynode} object corresponding to the \code{i}th
  kid.

  The \code{as.partynode} and \code{as.list} methods can be used
  to convert flat list structures into recursive \code{partynode}
  objects and vice versa. \code{as.partynode} applied to
  \code{partynode} objects renumbers the recursive nodes
  starting with root node identifier \code{from}.

  \code{length} gives the number of kid nodes of the root node,
  \code{depth} the depth of the tree and \code{width} 
  the number of terminal nodes.

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

## print raw recursive structure without data
print(node)

## print tree along with the associated iris data
data("iris", package = "datasets")
print(node, data = iris)

## print subtree
print(node[2], data = iris)

## print subtree, with root node number one
print(as.partynode(node[2], from = 1), data = iris)

## number of kids in root node
length(node)

## depth of tree
depth(node)

## number of terminal nodes
width(node)

## convert back to flat structure
as.list(node)
}
\keyword{tree}
