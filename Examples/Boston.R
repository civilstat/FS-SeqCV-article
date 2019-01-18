# Recreate right-hand subplots of Figure 6.6 from thesis

## Tong Zhang shared his FoBa paper's R code here:
## http://tongzhang-ml.org/software.html
## Source it and follow examples
## to replicate his graphs and add my own FS line.

source("./Examples/greedy.R")

library(MASS)
data(Boston)


## Mimic "forward-greedy" i.e. OMP in Fig. 5-6 of Zhang's FoBa paper
Boston$intercept = 1
B = 500 # Zhang used 50 reps, but I found it's too few to be stable
n = nrow(Boston)
p = ncol(Boston) - 1
nc = 50
nv = n - nc
MSEtrain_FS = matrix(0, B, 10)
MSEtest_FS = matrix(0, B, 10)
MSEtrain_FG = matrix(0, B, 10)
MSEtest_FG = matrix(0, B, 10)
MSEtrain_FoBa = matrix(0, B, 10)
MSEtest_FoBa = matrix(0, B, 10)
system.time({
  for(bb in 1:B){
    ii_train = sample(n, nc)
    train = Boston[ii_train, ]
    test = Boston[-ii_train, ]
    lmNow = lm(medv ~ 0, data = train)
    lmFull = lm(medv ~ ., data = train)
    model.FG = foba(as.matrix(train[, -14]),
                    train$medv, type = "for", steps = 20, intercept = FALSE) 
    model.FoBa = foba(as.matrix(train[, -14]),
                    train$medv, type = "foba.a", steps = 20, intercept = FALSE) 
    for(jj in 1:10) {
      addTerm = addterm(lmNow, scope = lmFull, sorted = TRUE, k = 0)
      newVar = ifelse(rownames(addTerm)[1] == "<none>",
                      rownames(addTerm)[2], rownames(addTerm)[1])
      lmNow = update(lmNow, formula(paste(". ~ . +", newVar)))
      MSEtrain_FS[bb, jj] = mean((predict(lmNow, newdata = train) - train$medv)^2)
      MSEtest_FS[bb, jj]  = mean((predict(lmNow, newdata = test) - test$medv)^2)

      MSEtrain_FG[bb, jj] = mean((predict(model.FG, as.matrix(train[, -14]),
                                          k = jj, type = "fit")$fit - train$medv)^2)
      MSEtest_FG[bb, jj] = mean((predict(model.FG, as.matrix(test[, -14]),
                                         k = jj, type = "fit")$fit - test$medv)^2)
      MSEtrain_FoBa[bb, jj] = mean((predict(model.FoBa, as.matrix(train[, -14]),
                                          k = jj, type = "fit")$fit - train$medv)^2)
      MSEtest_FoBa[bb, jj] = mean((predict(model.FoBa, as.matrix(test[, -14]),
                                         k = jj, type = "fit")$fit - test$medv)^2)
    }
  }
})

pdf("./Examples/Boston.pdf", width = 6, height = 4, pointsize = 10)
layout(matrix(1:2, 1))
plot(1:10, colMeans(MSEtrain_FS), type = 'b',
     ylab = "training error", xlab = "sparsity", log = "y",
     pch = 0, lty = 2)
lines(1:10, colMeans(MSEtrain_FG), type = 'b',
      pch = 3, lty = 3)
lines(1:10, colMeans(MSEtrain_FoBa), type = 'b',
      pch = 20, lty = 1)
legend("topright", legend = c("FoBa", "forward-greedy", "FS"),
       pch = c(20, 3, 0), lty = c(1, 3, 2))
plot(1:10, colMeans(MSEtest_FS), type = 'b',
     ylab = "test error", xlab = "sparsity", log = "y",
     pch = 0, lty = 2)
lines(1:10, colMeans(MSEtest_FG), type = 'b',
      pch = 3, lty = 3)
lines(1:10, colMeans(MSEtest_FoBa), type = 'b',
      pch = 20, lty = 1)
legend("topright", legend = c("FoBa", "forward-greedy", "FS"),
       pch = c(20, 3, 0), lty = c(1, 3, 2))
dev.off()

