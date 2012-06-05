\name{panelfunctions}
\alias{panelfunctions}
\alias{node_inner}
\alias{node_terminal}
\alias{edge_simple}
\alias{node_barplot}
\alias{node_boxplot}
\alias{node_surv}

\title{ Panel-Generators for Visualization of Party Trees }
\description{
    The plot method for \code{party} and \code{constparty} objects are rather
    flexible and can be extended by panel functions. Some pre-defined
    panel-generating functions of class \code{grapcon_generator}
    for the most important cases are documented here.
}
\usage{
node_inner(obj, id = TRUE, abbreviate = FALSE, fill = "white", 
    gp = gpar())
node_terminal(obj, digits = 3, abbreviate = FALSE,
    fill = c("lightgray", "white"), id = TRUE, gp = gpar(), FUN = NULL)
edge_simple(obj, digits = 3, abbreviate = FALSE)
node_boxplot(obj, col = "black", fill = "lightgray", width = 0.5,
    yscale = NULL, ylines = 3, cex = 0.5, id = TRUE, gp = gpar())
node_barplot(obj, col = "black", fill = NULL, beside = NULL,
    ymax = NULL, ylines = NULL, widths = 1, gap = NULL,
    reverse = NULL, id = TRUE, gp = gpar())
node_surv(obj, col = "black", ylines = 2, id = TRUE, gp = gpar(), \dots)
}
\arguments{
  \item{obj}{ an object of class \code{party}.}
  \item{digits}{ integer, used for formating numbers. }
  \item{abbreviate}{ logical indicating whether strings should be 
    abbreviated. }
  \item{col}{ a color for points and lines. }
  \item{fill}{ a color to filling rectangles. }
  \item{id}{ logical. Should node IDs be plotted?}
  \item{ylines}{ number of lines for spaces in y-direction. }
  \item{widths}{ widths in barplots. }
  \item{width}{ width in boxplots. }
  \item{gap}{ gap between bars in a barplot (\code{node_barplot}). }
  \item{yscale}{ limits in y-direction}
  \item{ymax}{ upper limit in y-direction}
  \item{cex}{character extension of points in scatter plots.}
  \item{beside}{ logical indicating if barplots should be side by side or stacked. }
  \item{reverse}{logical indicating whether the order of levels should be reversed
    for barplots.}
  \item{gp}{graphical parameters.}
  \item{FUN}{function for formatting the \code{info}, passed to \code{\link{formatinfo_node}}.}
  \item{\dots}{ additional arguments passed to callies 
                (for example to \code{\link[survival]{survfit}}).}
}

\details{
  The \code{plot} methods for \code{party} and \code{constparty} objects provide an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. The panel functions to be used
  should depend only on the node being visualized, however, for setting up
  an appropriate panel function, information from the whole tree is typically
  required. Hence, \pkg{party} adopts the framework of \code{grapcon_generator}
  (graphical appearance control) from the \pkg{vcd} package (Meyer, Zeileis and
  Hornik, 2005) and provides several panel-generating functions. For convenience,
  the panel-generating functions \code{node_inner} and \code{edge_simple} 
  return panel functions to draw inner nodes and left and right edges. 
  For drawing terminal nodes, the functions returned by the other panel 
  functions can be used. The panel generating function \code{node_terminal} 
  is a terse text-based representation of terminal nodes.

  Graphical representations of terminal nodes are available and depend on
  the kind of model and the measurement scale of the variables modeled.

  For univariate regressions (typically fitted by \code{}),
  \code{node_surv}
  returns a functions that plots Kaplan-Meier curves in each 
  terminal node; \code{node_barplot}, \code{node_boxplot}, \code{node_hist} and
  \code{node_density} can be used to plot bar plots, box plots, histograms and 
  estimated densities into the terminal nodes.

  For multivariate regressions (typically fitted by \code{mob}),
  \code{node_bivplot} returns a panel function that creates bivariate plots
  of the response against all regressors in the model. Depending on the scale
  of the variables involved, scatter plots, box plots, spinograms (or CD plots)
  and spine plots are created. For the latter two \code{\link[vcd]{spine}} and
  \code{\link[vcd]{cd_plot}} from the \pkg{vcd} package are re-used.
}

\references{
  David Meyer, Achim Zeileis, and Kurt Hornik (2006).
  The Strucplot Framework: Visualizing Multi-Way Contingency Tables with vcd.
  \emph{Journal of Statistical Software}, \bold{17}(3).
  \url{http://www.jstatsoft.org/v17/i03/}
}
\keyword{hplot}
