\name{WeatherPlay}
\alias{WeatherPlay}

\title{Weather Conditions and Playing a Game}

\description{
Artificial data set concerning the conditions suitable for
playing some unspecified game.
}

\usage{data("WeatherPlay")}

\format{
A data frame containing 14 observations on 5 variables.
\describe{
  \item{outlook}{factor.}
  \item{temperature}{numeric.}
  \item{humidity}{numeric.}
  \item{windy}{factor.}
  \item{play}{factor.}
}
}

\source{
Table 1.3 in Witten and Frank (2011).
}

\references{
Witten IH, Frank E (2011).
  \emph{Data Mining: Practical Machine Learning Tools and Techniques}.
  3rd Edition, Morgan Kaufmann, San Francisco. 
}

\seealso{\code{\link{party}}, \code{\link{partynode}}, \code{\link{partysplit}}}

\examples{
## load weather data
data("WeatherPlay", package = "partykit")
WeatherPlay

## construct simple tree
pn <- partynode(1L,
  split = partysplit(1L, index = 1:3),
  kids = list(
    partynode(2L,
      split = partysplit(3L, breaks = 75),
      kids = list(
        partynode(3L, info = "yes"),
        partynode(4L, info = "no"))),
    partynode(5L, info = "yes"),
    partynode(6L,
      split = partysplit(4L, index = 1:2),
      kids = list(
        partynode(7L, info = "yes"),
        partynode(8L, info = "no")))))
pn

## couple with data
py <- party(pn, WeatherPlay)

## print/plot/predict
print(py)
plot(py)
predict(py, newdata = WeatherPlay)

## customize printing
print(py,
  terminal_panel = function(node) paste(": play=", info_node(node), sep = ""))
}

\keyword{datasets}
