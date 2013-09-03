### R code from vignette source 'partykit.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width = 70)
library("partykit")
set.seed(290875)
data("iris")


###################################################
### code chunk number 2: weather-data
###################################################
data("WeatherPlay", package = "partykit")
WeatherPlay


###################################################
### code chunk number 3: weather-plot0
###################################################
py <- party(
  partynode(1L,
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
          partynode(8L, info = "no"))))),
  WeatherPlay)
plot(py)


###################################################
### code chunk number 4: weather-splits
###################################################
sp_o <- partysplit(1L, index = 1:3)
sp_h <- partysplit(3L, breaks = 75)
sp_w <- partysplit(4L, index = 1:2)


###################################################
### code chunk number 5: weather-nodes
###################################################
pn <- partynode(1L, split = sp_o, kids = list(
  partynode(2L, split = sp_h, kids = list(
    partynode(3L, info = "yes"),
    partynode(4L, info = "no"))),
  partynode(5L, info = "yes"),
  partynode(6L, split = sp_w, kids = list(
    partynode(7L, info = "yes"),
    partynode(8L, info = "no")))))


###################################################
### code chunk number 6: weather-nodes-print
###################################################
pn


###################################################
### code chunk number 7: weather-party
###################################################
py <- party(pn, WeatherPlay)
py


###################################################
### code chunk number 8: weahter-plot (eval = FALSE)
###################################################
## plot(py)


###################################################
### code chunk number 9: weahter-predict (eval = FALSE)
###################################################
## predict(py, newdata = WeatherPlay)


###################################################
### code chunk number 10: weather-party-print
###################################################
print(py,
  terminal_panel = function(node) paste(": play=", node$info, sep = ""))


###################################################
### code chunk number 11: partysplit-1
###################################################
## binary split in numeric variable `Sepal.Length'
sl5 <- partysplit(which(names(iris) == "Sepal.Length"), breaks = 5)
class(sl5)


###################################################
### code chunk number 12: partysplit-2
###################################################
unclass(sl5)


###################################################
### code chunk number 13: partysplit-3
###################################################
character_split(sl5, data = iris)


###################################################
### code chunk number 14: partysplit-4
###################################################
kidids_split(sl5, data = iris)


###################################################
### code chunk number 15: partysplit-5
###################################################
(!with(iris, Sepal.Length <= 5)) + 1


###################################################
### code chunk number 16: partysplit-6
###################################################
## binary split in factor `Species'
sp <- partysplit(which(names(iris) == "Species"), index = c(1L, 1L, 2L))
character_split(sp, data = iris)
table(kidids_split(sp, data = iris), iris$Species)


###################################################
### code chunk number 17: partysplit-6
###################################################
unclass(sp)


###################################################
### code chunk number 18: partysplit-7
###################################################
## multiway split in numeric variable `Sepal.Width',    
## higher values go to the first kid, smallest values
## to the last kid
sw23 <- partysplit(which(names(iris) == "Sepal.Width"),
  breaks = c(3, 3.5), index = 3:1)
character_split(sw23, data = iris)
table(kidids_split(sw23, data = iris),
  cut(iris$Sepal.Width, breaks = c(-Inf, 2, 3, Inf)))


###################################################
### code chunk number 19: partysplit-8
###################################################
sw23 <- partysplit(which(names(iris) == "Sepal.Width"),
  breaks = c(3, 3.5), index = c(1L, 3L, 2L))
character_split(sw23, data = iris)


###################################################
### code chunk number 20: partynode-1
###################################################
n1 <- partynode(id = 1L)
is.terminal(n1)
print(n1)


###################################################
### code chunk number 21: partynode-2
###################################################
n1 <- partynode(id = 1L, split = sl5, kids = sapply(2:3, partynode))
print(n1, data = iris)


###################################################
### code chunk number 22: partynode-3
###################################################
fitted_node(n1, data = iris)


###################################################
### code chunk number 23: partynode-4
###################################################
kidids_node(n1, data = iris)


###################################################
### code chunk number 24: partynode-5
###################################################
n1[2]


###################################################
### code chunk number 25: party-1
###################################################
t1 <- party(n1, 
  data = iris,
  fitted = data.frame(
    "(fitted)" = fitted_node(n1, data = iris),
    "(response)" = iris$Species,
    check.names = FALSE)
)
t1


###################################################
### code chunk number 26: mytree-1
###################################################
findsplit <- function(response, data, weights) {

  ### extract response values from data
  y <- data[[response]]

  logpmin <- 0
  xselect <- NULL

  ### cycle through all features
  for (i in which(names(data) != response)) {

    ### expand data
    x <- data[[i]]
    xt <- rep(x, weights)
    yt <- rep(y, weights)

    ### potential split points (not too many)
    qx <- unique(quantile(xt, 
        	 prob = seq(from = 0.1, to = 0.9, by = 0.05)))

    ### assess all potential splits by their t-test
    ### log-p-value
    logp <- sapply(qx, function(q) {
      tt <- t.test(yt[xt <= q], yt[xt > q])
      pt(-abs(tt$statistic), tt$parameter, log = TRUE) + log(2)
    })

    ### if the best split in variable i significant AND
    ### better than what we already had store this information
    if (min(logp) < logpmin & min(logp) < log(0.05)) {
      logpmin <- min(logp)
      xselect <- i
      splitpoint <- qx[which.min(logp)]
    }
  }

  ### no significant split found, give up
  if (is.null(xselect)) return(NULL)

  ### return split as partysplit object
  return(partysplit(
      varid = as.integer(xselect),	 ### which variable?
      breaks = as.numeric(splitpoint),   ### which split point?
      info = list(pvalue = exp(logpmin)  ### save p-value in addition
  )))
}


###################################################
### code chunk number 27: mytree-2
###################################################
growtree <- function(id = 1L, response, data, weights) {

  ### for less than 30 obs. stop here
  if (sum(weights) < 30) return(partynode(id = id))

  ### find best split
  sp <- findsplit(response, data, weights)
  ### no split found, stop here
  if (is.null(sp)) return(partynode(id = id))

  ### actually split the data
  kidids <- kidids_split(sp, data = data)

  ### set-up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids))
  for (kidid in 1:max(kidids)) {
  ### select obs for current node
  w <- weights
  w[kidids != kidid] <- 0
  ### get next node id
  if (kidid > 1) {
    myid <- max(nodeids(kids[[kidid - 1]]))
  } else {
    myid <- id
  }
  ### start recursion on this daugther node
  kids[[kidid]] <- growtree(id = as.integer(myid + 1), response, data, w)
  }

  ### return nodes
  return(partynode(id = as.integer(id), split = sp, kids = kids))
}


###################################################
### code chunk number 28: mytree-3
###################################################
mytree <- function(formula, data, weights = NULL) {

  ### name of the response variable
  response <- all.vars(formula)[1]
  ### data without missing values, response comes last
  data <- data[complete.cases(data), c(all.vars(formula)[-1], response)]
  ### data is numeric
  stopifnot(all(sapply(data, is.numeric)))

  if (is.null(weights)) weights <- rep(1, nrow(data))
  ### weights are case weights, i.e., integers
  stopifnot(length(weights) == nrow(data) &
    max(abs(weights - floor(weights))) < .Machine$double.eps)

  ### grow tree
  nodes <- growtree(id = 1L, response, data, weights)

  ### compute terminal node number for each obs.
  fitted <- fitted_node(nodes, data = data)
  ### return rich object
  ret <- party(nodes, 
    data = data,
    fitted = data.frame(
      "(fitted)" = fitted,
      "(response)" = data[[response]],
      "(weights)" = weights,
      check.names = FALSE),
    terms = terms(formula))
  as.constparty(ret)
}


###################################################
### code chunk number 29: mytree-4
###################################################
aqt <- mytree(Ozone ~ Solar.R + Wind + Temp, data = airquality)
aqt


###################################################
### code chunk number 30: mytree-5
###################################################
plot(aqt)


###################################################
### code chunk number 31: mytree-6
###################################################
predict(aqt, newdata = airquality[1:10,])


###################################################
### code chunk number 32: mytree-7
###################################################
aqt4 <- aqt[4]
aqt4


###################################################
### code chunk number 33: mytree-8
###################################################
plot(aqt4)


###################################################
### code chunk number 34: mytree-10
###################################################
fun <- function(x) format.pval(info_split(split_node(x))$pvalue,
  digits = 3, eps = 0.001)
nid <- nodeids(aqt)
iid <- nid[!(nid %in% nodeids(aqt, terminal = TRUE))]
unlist(nodeapply(aqt, ids = iid, FUN = fun))


