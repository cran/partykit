\name{HuntingSpiders}
\alias{HuntingSpiders}

\title{Abundance of Hunting Spiders}

\description{
Abundances for 12 species of hunting spiders along with
environmental predictors, all rated on a 0--9 scale.
}

\usage{data("HuntingSpiders")}

\format{
A data frame containing 28 observations on 18 variables
(12 species abundances and 6 environmental predictors).
\describe{
  \item{arct.lute}{numeric. Abundance of species \emph{Arctosa lutetiana} (on a scale 0--9).}
  \item{pard.lugu}{numeric. Abundance of species \emph{Pardosa lugubris} (on a scale 0--9).}
  \item{zora.spin}{numeric. Abundance of species \emph{Zora spinimana} (on a scale 0--9).}
  \item{pard.nigr}{numeric. Abundance of species \emph{Pardosa nigriceps} (on a scale 0--9).}
  \item{pard.pull}{numeric. Abundance of species \emph{Pardosa pullata} (on a scale 0--9).}
  \item{aulo.albi}{numeric. Abundance of species \emph{Aulonia albimana} (on a scale 0--9).}
  \item{troc.terr}{numeric. Abundance of species \emph{Trochosa terricola} (on a scale 0--9).}
  \item{alop.cune}{numeric. Abundance of species \emph{Alopecosa cuneata} (on a scale 0--9).}
  \item{pard.mont}{numeric. Abundance of species \emph{Pardosa monticola} (on a scale 0--9).}
  \item{alop.acce}{numeric. Abundance of species \emph{Alopecosa accentuata} (on a scale 0--9).}
  \item{alop.fabr}{numeric. Abundance of species \emph{Alopecosa fabrilis} (on a scale 0--9).}
  \item{arct.peri}{numeric. Abundance of species \emph{Arctosa perita} (on a scale 0--9).}
  \item{water}{numeric. Environmental predictor on a scale 0--9.}
  \item{sand}{numeric. Environmental predictor on a scale 0--9.}
  \item{moss}{numeric. Environmental predictor on a scale 0--9.}
  \item{reft}{numeric. Environmental predictor on a scale 0--9.}
  \item{twigs}{numeric. Environmental predictor on a scale 0--9.}
  \item{herbs}{numeric. Environmental predictor on a scale 0--9.}
}
}

\details{
  The data were originally analyzed by Van der Aart and Smeenk-Enserink (1975).
  De'ath (2002) transformed all variables to the 0--9 scale and employed
  multivariate regression trees.
}

\source{
  Package \pkg{mvpart} (currently archived, see
  \url{http://CRAN.R-project.org/package=mvpart}).
}

\references{
Van der Aart PJ, Smeenk-Enserink N (1975).
  Correlations between Distributions of Hunting Spiders (Lycosidae, Ctenidae)
  and Environmental Characteristics in a Dune Area.
  \emph{Netherlands Journal of Zoology}, \bold{25}, 1--45.

De'ath G (2002).
  Multivariate Regression Trees: A New Technique for Modelling Species-Environment Relationships.
  \emph{Ecology}, \bold{83}(4), 1103--1117.
}

\examples{
## load data
data("HuntingSpiders", package = "partykit")

## fit multivariate tree for 12-dimensional species abundance
sptree <- ctree(arct.lute + pard.lugu + zora.spin + pard.nigr + pard.pull +
  aulo.albi + troc.terr + alop.cune + pard.mont + alop.acce + alop.fabr +
  arct.peri ~ herbs + reft + moss + sand + twigs + water, data = HuntingSpiders,
  teststat = "max", minsplit = 5)
plot(sptree, terminal_panel = node_barplot)
}

\keyword{datasets}
