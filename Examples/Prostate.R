# Recreate left-hand subplot of Figure 6.6 from thesis,
# as well as numerical performance summaries reported in the text.

library(Hmisc)
prostate_full = read.table("./Examples/prostate.data")
prostate = subset(prostate_full, train)[, -10]
prostate_test = subset(prostate_full, !train)[, -10]
head(prostate)


## Try to mimic Fig. 3.7 of Hastie et al.'s ESL
V = 10
n = nrow(prostate)
p = ncol(prostate) - 1
MSEtrain_FS = matrix(0, V, p + 1)
MSEtest_FS = matrix(0, V, p + 1)
system.time({
  folds = sample(rep(1:V, length.out = n))
  for(vv in 1:V){
    ii_train = which(folds != vv)
    train = prostate[ii_train, ]
    test = prostate[-ii_train, ]
    lmNow = lm(lpsa ~ 1, data = train)
    lmFull = lm(lpsa ~ ., data = train)
    MSEtrain_FS[vv, 1] = mean((predict(lmNow, newdata = train) - train$lpsa)^2)
    MSEtest_FS[vv, 1]  = mean((predict(lmNow, newdata = test) - test$lpsa)^2)
    for(jj in 1:p) {
      addTerm = addterm(lmNow, scope = lmFull, sorted = TRUE, k = 0)
      newVar = ifelse(rownames(addTerm)[1] == "<none>",
                      rownames(addTerm)[2], rownames(addTerm)[1])
      lmNow = update(lmNow, formula(paste(". ~ . +", newVar)))
      MSEtrain_FS[vv, jj + 1] = mean((predict(lmNow, newdata = train) - train$lpsa)^2)
      MSEtest_FS[vv, jj + 1]  = mean((predict(lmNow, newdata = test) - test$lpsa)^2)
    }
  }
})
colSEs = function(mat) {
  apply(mat, 2, function(x) sd(x)/sqrt(length(x)))
}

pdf("./Examples/Prostate.pdf", width = 4, height = 4, pointsize = 10)
layout(1)
errbar(0:p, colMeans(MSEtest_FS), colMeans(MSEtest_FS) + colSEs(MSEtest_FS),
       colMeans(MSEtest_FS) - colSEs(MSEtest_FS), type = 'p',
       ylab = "test error", xlab = "subset size",
       pch = 0, lty = 0)
lines(0:p, colMeans(MSEtest_FS), lty = 2)
errbar(0:p, colMeans(MSEtest_FS), colMeans(MSEtest_FS) + colSEs(MSEtest_FS),
       colMeans(MSEtest_FS) - colSEs(MSEtest_FS), add = TRUE, pch = NA)
abline(h = (colMeans(MSEtest_FS)+colSEs(MSEtest_FS))[which.min(colMeans(MSEtest_FS))],
       lty = 3)
abline(v = 2, lty = 3)
legend("topright", legend = c("FS", "1SE rule"),
       pch = c(0, NA), lty = c(2, 3))
dev.off()


## Replicate parts of Table 3.3;
## and also check out FS with k=7, the global minimizer of CV error
## (in addition to k=2 which wins according to SeqCV or the 1SE rule)
lm0 = lm(lpsa ~ 1, data = prostate)

## Full OLS model
lmFull = lm(lpsa ~ ., data = prostate)
mean((predict(lmFull, newdata = prostate_test) - prostate_test$lpsa)^2)
## 0.521, matches ESL

## SeqCV or CV+1SE
lm2 = stepAIC(lm0, scope = list(lower=lm0, upper=lmFull),
              direction="forward", k = 0, steps = 2)
mean((predict(lm2, newdata = prostate_test) - prostate_test$lpsa)^2)
## 0.492, matches ESL

## Standard CV (global min)
lm7 = stepAIC(lm0, scope = list(lower=lm0, upper=lmFull),
              direction="forward", k = 0, steps = 7)
mean((predict(lm7, newdata = prostate_test) - prostate_test$lpsa)^2)
## 0.517, not reported in ESL



