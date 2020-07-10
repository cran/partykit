## ---- Data -------------------------------------------------------------------

## Titanic survival data from base R
data("Titanic", package = "datasets")

## turn four-way contingency table into long data frame
ttnc <- as.data.frame(Titanic)
ttnc <- ttnc[rep(1:nrow(ttnc), ttnc$Freq), 1:4]
names(ttnc)[2] <- "Gender"

## wage determinants data from "Applied Econometrics with R"
data("CPS1985", package = "AER")


## ---- Motivational tree examples ---------------------------------------------

## Titanic survival: CTree
library("partykit")
ct_ttnc <- ctree(Survived ~ Gender + Age + Class, data = ttnc, alpha = 0.01)
plot(ct_ttnc)

## Wage determinants: CTree
ct_cps <- ctree(log(wage) ~ education + experience + age + ethnicity + gender + union, data = CPS1985, alpha = 0.01)
plot(ct_cps)

## Wage determinants: MOB (based on lm)
mob_cps <- lmtree(log(wage) ~ education | experience + age + ethnicity + gender + union, data = CPS1985)
plot(mob_cps)


## ---- Determine first split using classical statistical inference ------------

## Titanic survival
## - Gender
plot(Survived ~ Gender, data = ttnc)
tab <- xtabs(~ Survived + Gender, data = ttnc)
chisq.test(tab)
## - Age
plot(Survived ~ Age, data = ttnc)
tab <- xtabs(~ Survived + Age, data = ttnc)
chisq.test(tab)
## - Class
plot(Survived ~ Class, data = ttnc)
tab <- xtabs(~ Survived + Class, data = ttnc)
chisq.test(tab)

## Wage determinants
## - education
plot(log(wage) ~ education, data = CPS1985)
cor.test(~ log(wage) + education, data = CPS1985)
## - gender
plot(log(wage) ~ gender, data = CPS1985)
t.test(log(wage) ~ gender, data = CPS1985)


## ---- Conditional inference trees --------------------------------------------

## Refit CTree from motivation section
library("partykit")
ct_ttnc <- ctree(Survived ~ Gender + Age + Class, data = ttnc)
plot(ct_ttnc)
print(ct_ttnc)

## Predictions
ndm <- data.frame(Gender = "Male", Age = "Adult", Class = c("1st", "2nd", "3rd"))
predict(ct_ttnc, newdata = ndm, type = "node")
predict(ct_ttnc, newdata = ndm, type = "response")
predict(ct_ttnc, newdata = ndm, type = "prob")

## Women and children first?
ndf <- data.frame(Gender = "Female", Age = "Adult", Class = c("1st", "2nd", "3rd"))
ndc <- data.frame(Gender = "Male", Age = "Child", Class = c("1st", "2nd", "3rd"))
cbind(
  Male   = predict(ct_ttnc, newdata = ndm, type = "prob")[, 2],
  Female = predict(ct_ttnc, newdata = ndf, type = "prob")[, 2],
  Child  = predict(ct_ttnc, newdata = ndc, type = "prob")[, 2]
)

## Refined tree (allowing small subgroups for children)
ct_ttnc2 <- ctree(Survived ~ Gender + Age + Class, data = ttnc,
  alpha = 0.01, minbucket = 5, minsplit = 15, maxdepth = 4)
plot(ct_ttnc2)

## New predictions
predict(ct_ttnc2, newdata = ndc, type = "prob")

## Fitted class labels and nodes from the tree
ttnc$Fit <- predict(ct_ttnc2, type = "response")
ttnc$Group <- factor(predict(ct_ttnc2, type = "node"))

## Confusion matrix
xtabs(~ Fit + Survived, data = ttnc)

## Group-specific survival proportions
tab <- xtabs(~ Group + Survived, data = ttnc)
prop.table(tab, 1)

## ROC analysis
## - Set up predictions
library("ROCR")
pred <- prediction(predict(ct_ttnc, type = "prob")[, 2], ttnc$Survived)
## - Accurcacies
plot(performance(pred, "acc"))
## - ROC curve
plot(performance(pred, "tpr", "fpr"))
abline(0, 1, lty = 2)


## ---- Recursive partitioning -------------------------------------------------


## Fit RPart/CART tree to Titanic survival data
library("rpart")
rp_ttnc <- rpart(Survived ~ Gender + Age + Class, data = ttnc)

## Display from rpart package
plot(rp_ttnc)
text(rp_ttnc)
print(rp_ttnc)

## Leverage display from partykit package
py_ttnc <- as.party(rp_ttnc)
plot(py_ttnc)
print(py_ttnc)

## Results from pruning and cross-validation
rp_ttnc$cptable

## Prune tree (suboptimally here)
prune(rp_ttnc, cp = 0.1)


## ---- Model-based recursive partitioning -------------------------------------

## Add preferential treatment variable to data
ttnc <- transform(ttnc,
  Treatment = factor(Gender == "Female" | Age == "Child", 
    levels = c(FALSE, TRUE), labels = c("Male&Adult", "Female|Child")
  )
)

## Fit and display model-based tree estimating heterogenous treatment effects
mob_ttnc <- glmtree(Survived ~ Treatment | Class + Gender + Age,
  data = ttnc, family = binomial, alpha = 0.01)
plot(mob_ttnc)
print(mob_ttnc)


## ---- Evolutionary learning of globally optimal trees ------------------------

## Package and random seed for reproducibility
library("evtree")
set.seed(1)

## Fit and display tree
ev_ttnc <- evtree(Survived ~ Gender + Age + Class, data = ttnc)
plot(ev_ttnc)
ev_ttnc


## ---- ggplot2 visualizations -------------------------------------------------

## Load package and set theme
library("ggparty")
theme_set(theme_minimal())
## Trees from above
autoplot(ct_ttnc2)
autoplot(ct_cps)
autoplot(py_ttnc)
autoplot(ev_ttnc)
