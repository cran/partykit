\name{partysplit}
\alias{partysplit}
\alias{kidids_split}
\alias{character_split}
\alias{varid_split}
\alias{breaks_split}
\alias{index_split}
\alias{right_split}
\alias{prob_split}
\alias{info_split}
\title{ Binary and Multiway Splits }
\description{

  A class for representing multiway splits and functions 
  for computing on splits.

}
\usage{
partysplit(varid, breaks = NULL, index = NULL, right = TRUE, 
    prob = NULL, info = NULL)
kidids_split(split, data, vmatch = 1:length(data), obs = NULL)
character_split(split, data = NULL, 
    digits = getOption("digits") - 2)
varid_split(split)
breaks_split(split)
index_split(split)
right_split(split)
prob_split(split)
info_split(split)
}
\arguments{
  \item{varid}{ an integer specifying the variable to split in, i.e., 
                a column number in \code{data}. }
  \item{breaks}{ a numeric vector of split points. }
  \item{index}{ an integer vector containing a contiguous sequence 
                from one to the number of kid nodes. May contain \code{NA}s.}
  \item{right}{ a logical, indicating if the intervals defined by 
                \code{breaks} should be closed on the
                right (and open on the left) or vice versa.}
  \item{prob}{ a numeric vector representing a probability distribution 
               over kid nodes. }
  \item{info}{ additional information. }
  \item{split}{ an object of class \code{partysplit}.}
  \item{data}{ a \code{\link{list}} or \code{\link{data.frame}}.}
  \item{vmatch}{ a permutation of the variable numbers in \code{data}.}
  \item{obs}{ a logical or integer vector indicating a subset of the 
              observations in \code{data}.}
  \item{digits}{ minimal number of significant digits.}
}
\details{

  A split is basically a function that maps data, 
  more specifically a partitioning variable,
  to a set of integers indicating the kid nodes to send observations to. 
  Objects of class \code{partysplit} describe such a function and can
  be set-up via the \code{partysplit()} constructor.
  The variables are available in a \code{list} or \code{data.frame} 
  (here called \code{data}) and \code{varid} specifies the 
  partitioning variable, i.e., the variable or list element to split in. 
  The constructor \code{partysplit()} doesn't have access 
  to the actual data, i.e., doesn't \emph{estimate} splits.

  \code{kidids_split(split, data)} actually partitions the data
  \code{data[obs,varid_split(split)]} and assigns an integer (giving the
  kid node number) to each observation. If \code{vmatch} is given,
  the variable \code{vmatch[varid_split(split)]} is used.
  
  \code{character_split()} returns a character representation
  of its \code{split} argument. The remaining functions 
  defined here are accessor functions for \code{partysplit} 
  objects.

  The numeric vector \code{breaks} defines how the range of 
  the partitioning variable (after coercing to a numeric via 
  \code{\link{as.numeric}}) is divided into intervals 
  (like in \code{\link{cut}}) and may be 
  \code{NULL}. These intervals are represented by the
  numbers one to \code{length(breaks) + 1}.

  \code{index} assigns these \code{length(breaks) + 1} 
  intervals to one of at least two kid nodes. Thus, \code{index}
  is a vector of integers where each element corresponds
  to one element in a list \code{kids} containing \code{\link{partynode}}
  objects, see \code{\link{partynode}} for details. The vector
  \code{index} may contain \code{NA}s, in that case, the corresponding
  values of the splitting variable are treated as missings (for
  example factor levels that are not present in the learning sample).
  Either \code{breaks} or \code{index} must be given.
  When \code{breaks} is \code{NULL}, it is assumed that
  the partitioning variable itself has storage mode \code{integer}
  (e.g., is a \code{\link{factor}}).

  \code{prob} defines a probability distribution over
  all kid nodes which is used for random splitting
  when a deterministic split isn't possible (due to missing 
  values, for example).

  \code{info} takes arbitrary user-specified information.

}
\value{
  The constructor \code{partysplit()} returns an object of class \code{partysplit}:
  \item{varid}{ an integer specifying the variable to split in, i.e., 
                a column number in \code{data}, }
  \item{breaks}{ a numeric vector of split points, }
  \item{index}{ an integer vector containing a contiguous sequence 
                from one to the number of kid nodes,}
  \item{right}{ a logical, indicating if the intervals defined by 
                \code{breaks} should be closed on the
                right (and open on the left) or vice versa}
  \item{prob}{ a numeric vector representing a probability distribution 
               over kid nodes, }
  \item{info}{ additional information. }

  \code{kidids_split()} returns an integer vector describing
  the partition of the observations into kid nodes.

  \code{character_split()} gives a character representation of the
  split and the remaining functions return the corresponding slots
  of \code{partysplit} objects.
}
\seealso{\code{\link{cut}}}
\references{ 
Hothorn T, Zeileis A (2015).
  partykit: A Modular Toolkit for Recursive Partytioning in R.
  \emph{Journal of Machine Learning Research}, \bold{16}, 3905--3909.
}
\examples{
data("iris", package = "datasets")

## binary split in numeric variable `Sepal.Length'
sl5 <- partysplit(which(names(iris) == "Sepal.Length"),
    breaks = 5)
character_split(sl5, data = iris)
table(kidids_split(sl5, data = iris), iris$Sepal.Length <= 5)

## multiway split in numeric variable `Sepal.Width', 
## higher values go to the first kid, smallest values 
## to the last kid
sw23 <- partysplit(which(names(iris) == "Sepal.Width"),    
    breaks = c(3, 3.5), index = 3:1)	
character_split(sw23, data = iris)    
table(kidids_split(sw23, data = iris), 
    cut(iris$Sepal.Width, breaks = c(-Inf, 2, 3, Inf)))   

## binary split in factor `Species'
sp <- partysplit(which(names(iris) == "Species"),
    index = c(1L, 1L, 2L))
character_split(sp, data = iris)
table(kidids_split(sp, data = iris), iris$Species)

## multiway split in factor `Species'
sp <- partysplit(which(names(iris) == "Species"), index = 1:3)
character_split(sp, data = iris)
table(kidids_split(sp, data = iris), iris$Species)

## multiway split in numeric variable `Sepal.Width'
sp <- partysplit(which(names(iris) == "Sepal.Width"), 
    breaks = quantile(iris$Sepal.Width))
character_split(sp, data = iris)

}
\keyword{tree}
