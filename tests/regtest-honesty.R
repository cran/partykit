
library("partykit")

n <- 100
x <- runif(n)
y <- rnorm(n, mean = sin(x), sd = .1)
s <- gl(4, n / 4)

set.seed(29)
### estimate with honesty
cf_ss <- cforest(y ~ x, strata = s, ntree = 5, mtry = 1,
                 perturb = list(replace = FALSE, fraction = c(.5, .5)))
### sample used for tree induction
stopifnot(sum(tapply(cf_ss$weights[[1]], s, sum)) == n / 2)
### sample used for parameter estimation
stopifnot(sum(tapply(cf_ss$honest_weights[[1]], s, sum)) == n / 2)

p <- predict(cf_ss)

set.seed(29)
### w/o honesty
cf_ss2 <- cforest(y ~ x, strata = s, ntree = 5, mtry = 1,
                  perturb = list(replace = FALSE, fraction = .5))

stopifnot(sum(tapply(cf_ss2$weights[[1]], s, sum)) == n / 2)

stopifnot(all.equal(cf_ss$nodes, cf_ss2$nodes))
stopifnot(all.equal(cf_ss$weights, cf_ss2$weights))
stopifnot(all.equal(predict(cf_ss, type = "node"), 
                    predict(cf_ss2, type = "node")))

tmp <- cf_ss2
tmp$weights <- lapply(tmp$weights, function(x) 1L - x)
pp <- predict(tmp)

stopifnot(all.equal(p, pp))

### bootstrap ignores honesty
cf_bs <- cforest(y ~ x, strata = s, ntree = 5,
              perturb = list(replace = TRUE, fraction = c(.5, .5)))

stopifnot(sum(tapply(cf_bs$weights[[1]], s, sum)) == n)
stopifnot(is.null(cf_bs$honest_weights))
