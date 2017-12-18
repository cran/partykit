\name{party-plot}
\alias{party-plot}
\alias{plot.party}
\alias{plot.constparty}
\alias{plot.simpleparty}
\title{ Visualization of Trees }

\description{
  \code{plot} method for \code{party} objects with
  extended facilities for plugging in panel functions.
}

\usage{
\method{plot}{party}(x, main = NULL,
    terminal_panel = node_terminal, tp_args = list(),
    inner_panel = node_inner, ip_args = list(),
    edge_panel = edge_simple, ep_args = list(),
    drop_terminal = FALSE, tnex = 1,
    newpage = TRUE, pop = TRUE, gp = gpar(),
    margins = NULL, \dots)
\method{plot}{constparty}(x, main = NULL,
    terminal_panel = NULL, tp_args = list(),
    inner_panel = node_inner, ip_args = list(),
    edge_panel = edge_simple, ep_args = list(),
    type = c("extended", "simple"), drop_terminal = NULL, 
    tnex = NULL, newpage = TRUE, pop = TRUE, gp = gpar(), 
    \dots)
\method{plot}{simpleparty}(x, digits = getOption("digits") - 4, tp_args = NULL, \dots)
}
\arguments{
  \item{x}{ an object of class \code{party} or \code{constparty}.}
  \item{main}{ an optional title for the plot.}
  \item{type}{ a character specifying the complexity of the plot:
               \code{extended} tries to visualize the distribution of the
               response variable in each terminal node whereas \code{simple} 
               only gives some summary information.}
  \item{terminal_panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the terminal nodes.
                         Alternatively, a panel generating function of class
			 \code{"grapcon_generator"} that is called with arguments
			 \code{x} and \code{tp_args} to set up a panel function.
			 By default, an appropriate panel function is chosen 
			 depending on the scale of the dependent variable.}
  \item{tp_args}{ a list of arguments passed to \code{terminal_panel} if this
                  is a \code{"grapcon_generator"} object.}
  \item{inner_panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the inner nodes.
                         Alternatively, a panel generating function of class
			 \code{"grapcon_generator"} that is called with arguments
			 \code{x} and \code{ip_args} to set up a panel function.}
  \item{ip_args}{ a list of arguments passed to \code{inner_panel} if this
                  is a \code{"grapcon_generator"} object.}
  \item{edge_panel}{ an optional panel function of the form 
                    \code{function(split, ordered = FALSE, left = TRUE)}
                    plotting the edges.
                    Alternatively, a panel generating function of class
		    \code{"grapcon_generator"} that is called with arguments
		    \code{x} and \code{ip_args} to set up a panel function.}		       
  \item{ep_args}{ a list of arguments passed to \code{edge_panel} if this
                  is a \code{"grapcon_generator"} object.}
  \item{drop_terminal}{ a logical indicating whether all terminal nodes
                        should be plotted at the bottom.}
  \item{tnex}{a numeric value giving the terminal node extension in relation
    to the inner nodes.}
  \item{newpage}{ a logical indicating whether \code{grid.newpage()} should be called. }
  \item{pop}{ a logical whether the viewport tree should be popped before
              return. }
  \item{gp}{graphical parameters.}
  \item{margins}{numeric vector of margin sizes.}
  \item{digits}{number of digits to be printed.}
  \item{\dots}{ additional arguments passed to callies.}
}

\details{
  This \code{plot} method for \code{party} objects provides an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. Panel functions for plotting
  inner nodes, edges and terminal nodes are available for the most important
  cases and can serve as the basis for user-supplied extensions, see
  \code{\link{node_inner}}.

  More details on the ideas and concepts of panel-generating functions and
  \code{"grapcon_generator"} objects in general can be found in Meyer, Zeileis
  and Hornik (2005).
}

\references{
Meyer D, Zeileis A, Hornik K (2006).
  The Strucplot Framework: Visualizing Multi-Way Contingency Tables with vcd.
  \emph{Journal of Statistical Software}, \bold{17}(3), 1--48.
  \url{http://www.jstatsoft.org/v17/i03/}
}


\seealso{\code{\link{node_inner}}, \code{\link{node_terminal}}, \code{\link{edge_simple}}, 
  \code{\link{node_barplot}}, \code{\link{node_boxplot}}.}

\keyword{hplot}
