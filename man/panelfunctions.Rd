\name{panelfunctions}
\alias{panelfunctions}
\alias{node_inner}
\alias{node_terminal}
\alias{edge_simple}
\alias{node_barplot}
\alias{node_bivplot}
\alias{node_boxplot}
\alias{node_surv}
\alias{node_ecdf}
\alias{node_mvar}

\title{ Panel-Generators for Visualization of Party Trees }

\description{
    The plot method for \code{party} and \code{constparty} objects are rather
    flexible and can be extended by panel functions. Some pre-defined
    panel-generating functions of class \code{grapcon_generator}
    for the most important cases are documented here.
}

\usage{
node_inner(obj, id = TRUE, pval = TRUE, abbreviate = FALSE, fill = "white", 
  gp = gpar())

node_terminal(obj, digits = 3, abbreviate = FALSE,
  fill = c("lightgray", "white"), id = TRUE,
  just = c("center", "top"), top = 0.85,
  align = c("center", "left", "right"), gp = NULL, FUN = NULL,
  height = NULL, width = NULL)

edge_simple(obj, digits = 3, abbreviate = FALSE, justmin = Inf,
  just = c("alternate", "increasing", "decreasing", "equal"),
  fill = "white")

node_boxplot(obj, col = "black", fill = "lightgray", bg = "white", width = 0.5,
  yscale = NULL, ylines = 3, cex = 0.5, id = TRUE, mainlab = NULL, gp = gpar())

node_barplot(obj, col = "black", fill = NULL, bg = "white",
  beside = NULL, ymax = NULL, ylines = NULL, widths = 1, gap = NULL,
  reverse = NULL, rot = 0, just = c("center", "top"), id = TRUE,
  mainlab = NULL, text = c("none", "horizontal", "vertical"), gp = gpar())

node_surv(obj, col = "black", bg = "white", yscale = c(0, 1), ylines = 2,
  id = TRUE, mainlab = NULL, gp = gpar(), \dots)

node_ecdf(obj, col = "black", bg = "white", ylines = 2,
  id = TRUE, mainlab = NULL, gp = gpar(), \dots)

node_bivplot(mobobj, which = NULL, id = TRUE, pop = TRUE,
  pointcol = "black", pointcex = 0.5,
  boxcol = "black", boxwidth = 0.5, boxfill = "lightgray",
  bg = "white", fitmean = TRUE, linecol = "red",
  cdplot = FALSE, fivenum = TRUE, breaks = NULL,
  ylines = NULL, xlab = FALSE, ylab = FALSE, margins = rep(1.5, 4),
  mainlab = NULL, \dots)

node_mvar(obj, which = NULL, id = TRUE, pop = TRUE, ylines = NULL,
  mainlab = NULL, varlab = TRUE, bg = "white", ...)
}

\arguments{
  \item{obj}{ an object of class \code{party}.}
  \item{digits}{ integer, used for formating numbers. }
  \item{abbreviate}{ logical indicating whether strings should be 
    abbreviated. }
  \item{col, pointcol, boxcol, linecol}{ a color for points and lines. }
  \item{fill, boxfill, bg}{ a color to filling rectangles and backgrounds. }
  \item{id}{ logical. Should node IDs be plotted?}
  \item{pval}{ logical. Should node p values be plotted (if they are
    available)?}
  \item{just}{justification of terminal panel viewport (\code{node_terminal}),
    or labels (\code{edge_simple}, \code{node_barplot}).}
  \item{justmin}{minimum average edge label length to employ justification
    via \code{just} in \code{edge_panel}, otherwise \code{just = "equal"}
    is used. Thus, by default \code{"equal"} justification is always used
    but other justifications could be employed for finite \code{justmin}.}
  \item{top}{in case of top justification, the npc coordinate at which the
    viewport is justified.}
  \item{align}{alignment of text within terminal panel viewport.}
  \item{ylines}{ number of lines for spaces in y-direction. }
  \item{widths}{ widths in barplots. }
  \item{boxwidth}{ width in boxplots (called \code{width} in \code{node_boxplot}). }
  \item{gap}{ gap between bars in a barplot (\code{node_barplot}). }
  \item{yscale}{ limits in y-direction}
  \item{ymax}{ upper limit in y-direction}
  \item{cex, pointcex}{character extension of points in scatter plots.}
  \item{beside}{ logical indicating if barplots should be side by side or stacked. }
  \item{reverse}{logical indicating whether the order of levels should be reversed
    for barplots. }
  \item{rot}{ arguments passed to \code{\link[grid]{grid.text}} for the
    x-axis labeling. }
  \item{gp}{graphical parameters.}
  \item{FUN}{function for formatting the \code{info}, passed to \code{\link{formatinfo_node}}.}
  \item{height, width}{ numeric, number of lines/columns for printing text. }
  \item{mobobj}{an object of class \code{modelparty} as computed by \code{\link{mob}}.}
  \item{which}{numeric or character. Optional selection of subset of regressor
    variables. By default one panel for each regressor variable is drawn.}
  \item{pop}{logical. Should the viewports in the individual nodes be popped after
    drawing?}
  \item{fitmean}{logical. Should the fitted mean function be visualized?}
  \item{cdplot}{logical. Should a CD plot (or a spineplot) be drawn when
    the response variable is categorical?}
  \item{fivenum}{logical. Should the five-number summary be used for splitting
    the x-axis in spineplots?}
  \item{breaks}{numeric. Optional numeric vector with breaks for the x-axis in
    splineplots.}
  \item{xlab, ylab}{character. Optional annotation for x-axis and y-axis.}
  \item{margins}{numeric. Margins around drawing area in viewport.}
  \item{mainlab}{character or function. An optional title for the plot. Either
    a character or a \code{function(id, nobs)}.}
  \item{varlab}{logical. Should the individual variable labels be attached to
    the \code{mainlab} for multivariate responses?}
  \item{text}{logical or character. Should percentage labels be drawn
    for each bar? The default is \code{"none"} or equivalently \code{FALSE}.
    Can be set to \code{TRUE} (or \code{"horizontal"}) or alternatively
    \code{"vertical"}.}
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
  terminal node; \code{node_barplot}, \code{node_boxplot}, \code{node_hist}, \code{node_ecdf} and
  \code{node_density} can be used to plot bar plots, box plots, histograms, 
  empirical cumulative distribution functions and 
  estimated densities into the terminal nodes.

  For multivariate regressions (typically fitted by \code{mob}),
  \code{node_bivplot} returns a panel function that creates bivariate plots
  of the response against all regressors in the model. Depending on the scale
  of the variables involved, scatter plots, box plots, spinograms (or CD plots)
  and spine plots are created. For the latter two \code{\link[vcd]{spine}} and
  \code{\link[vcd]{cd_plot}} from the \pkg{vcd} package are re-used.

  For multivariate responses in \code{\link{ctree}}, the panel function \code{node_mvar}
  generates one plot for each response.
}

\references{
Meyer D, Zeileis A, Hornik K (2006).
  The Strucplot Framework: Visualizing Multi-Way Contingency Tables with vcd.
  \emph{Journal of Statistical Software}, \bold{17}(3), 1--48.
  \url{http://www.jstatsoft.org/v17/i03/}
}

\keyword{hplot}
