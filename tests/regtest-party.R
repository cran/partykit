## load package and fix seed
library("partykit")
library("survival")
set.seed(1)

## rpart: kyphosis data
library("rpart")
data("kyphosis", package = "rpart")
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
pfit <- as.party(fit)
all(predict(pfit, newdata = kyphosis, type = "node") == fit$where)

## J48: iris data
library("RWeka")
data("iris", package = "datasets")
itree <- J48(Species ~ ., data = iris)
pitree <- as.party(itree)
stopifnot(all(predict(pitree) == predict(pitree, newdata = iris[, 3:4])))

all.equal(predict(itree, type = "prob", newdata = iris),
          predict(pitree, type = "prob", newdata = iris))
all.equal(predict(itree,  newdata = iris),
          predict(pitree, newdata = iris))

## rpart/J48: GlaucomaM data
data("GlaucomaM", package = "TH.data")
w <- runif(nrow(GlaucomaM))
fit <- rpart(Class ~ ., data = GlaucomaM, weights = w)
pfit <- as.party(fit)
all(predict(pfit, type = "node") == fit$where)
tmp <- GlaucomaM[sample(1:nrow(GlaucomaM), 100),]
all.equal(predict(fit, type = "prob", newdata = tmp), predict(pfit, type = "prob", newdata = tmp))
all.equal(predict(fit, type = "class", newdata = tmp), predict(pfit, newdata = tmp))
itree <- J48(Class ~ ., data = GlaucomaM)
pitree <- as.party(itree)
all.equal(predict(itree, newdata = tmp, type = "prob"), predict(pitree, newdata = tmp, type = "prob"))

## rpart: airquality data
data("airquality")
aq <- subset(airquality, !is.na(Ozone))
w <- runif(nrow(aq), max = 3)
aqr <- rpart(Ozone ~ ., data = aq, weights = w)
aqp <- as.party(aqr)
tmp <- subset(airquality, is.na(Ozone))
all.equal(predict(aqr, newdata = tmp), predict(aqp, newdata = tmp))

## rpart: GBSG2 data
data("GBSG2", package = "TH.data")
library("survival")
fit <- rpart(Surv(time, cens) ~ ., data = GBSG2)
pfit <- as.party(fit)
pfit$fitted
predict(pfit, newdata = GBSG2[1:100,], type = "prob")
predict(pfit, newdata = GBSG2[1:100,], type = "response")

### multiple responses
f <- fitted(pfit)
f[["(response)"]] <- data.frame(srv = f[["(response)"]], hansi = runif(nrow(f)))
mp <- party(node_party(pfit), fitted = f, data = pfit$data)
class(mp) <- c("constparty", "party")
predict(mp, newdata = GBSG2[1:10,])

