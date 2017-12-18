
library("partykit")
stopifnot(require("party"))
set.seed(29)

### regression
airq <- airquality[complete.cases(airquality),]

mtry <- ncol(airq) - 1L
ntree <- 25

cf_partykit <- partykit::cforest(Ozone ~ ., data = airq,
    ntree = ntree, mtry = mtry)

w <- do.call("cbind", cf_partykit$weights)

cf_party <- party::cforest(Ozone ~ ., data = airq, 
    control = party::cforest_unbiased(ntree = ntree, mtry = mtry),
    weights = w)

p_partykit <- predict(cf_partykit)
p_party <- predict(cf_party)

stopifnot(max(abs(p_partykit - p_party)) < .Machine$double.eps)

prettytree(cf_party@ensemble[[1]], inames = names(airq)[-1])
party(cf_partykit$nodes[[1]], data = model.frame(cf_partykit))

v_party <- do.call("rbind", lapply(1:5, function(i) party::varimp(cf_party)))

v_partykit <- do.call("rbind", lapply(1:5, function(i) partykit::varimp(cf_partykit)))

summary(v_party)
summary(v_partykit)

party::varimp(cf_party, conditional = TRUE)
partykit::varimp(cf_partykit, conditional = TRUE)

### classification
set.seed(29)
mtry <- ncol(iris) - 1L
ntree <- 25

cf_partykit <- partykit::cforest(Species ~ ., data = iris,
    ntree = ntree, mtry = mtry)

w <- do.call("cbind", cf_partykit$weights)

cf_party <- party::cforest(Species ~ ., data = iris, 
    control = party::cforest_unbiased(ntree = ntree, mtry = mtry),
    weights = w)

p_partykit <- predict(cf_partykit, type = "prob")
p_party <- do.call("rbind", treeresponse(cf_party))

stopifnot(max(abs(unclass(p_partykit) - unclass(p_party))) < .Machine$double.eps)

prettytree(cf_party@ensemble[[1]], inames = names(iris)[-5])
party(cf_partykit$nodes[[1]], data = model.frame(cf_partykit))

v_party <- do.call("rbind", lapply(1:5, function(i) party::varimp(cf_party)))

v_partykit <- do.call("rbind", lapply(1:5, function(i)
    partykit::varimp(cf_partykit, risk = "mis")))

summary(v_party)
summary(v_partykit)

party::varimp(cf_party, conditional = TRUE)
partykit::varimp(cf_partykit, risk = "misclass", conditional = TRUE)
