suppressWarnings(RNGversion("3.5.2"))

## package
library("partykit")

## iris data
data("iris", package = "datasets")
irisct <- ctree(Species ~ ., data = iris)
print(irisct)
table(fit = predict(irisct), true = iris$Species)

## airquality data
data("airquality", package = "datasets")
airq <- subset(airquality, !is.na(Ozone))
airqct <- ctree(Ozone ~ ., data = airq)
print(airqct)
sum((airq$Ozone - predict(airqct))^2)

### split in one variable only: Temp is selected freely in the root node
### but none of the other variables is allowed deeper in the tree
airqct1 <- ctree(Ozone ~ ., data = airq, 
                 control = ctree_control(maxvar = 1L))
psplitids <- unique(do.call("c", 
        nodeapply(node_party(airqct1), 
                  ids = nodeids(node_party(airqct1)),
                  FUN = function(x) split_node(x)$varid)))
stopifnot(length(psplitids) == 1L)
