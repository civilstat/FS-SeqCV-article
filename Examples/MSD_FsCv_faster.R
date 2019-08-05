####
## Recreate parts of Table 6.1 from thesis
## as well as numerical results mentioned inline in Section 6.2.2
## using a 10:90 train:test split, and also a 90:10 split,
## both within the original learning-subset of the data



## July 29-30, 2019

## Revising this code based on reviewer feedback that it's VERY slow.
## Instead of doing it sequentially with R's basic MASS:addterm()
## and all the machinery of creating lm objects,
## now we use leaps package to generate the whole path very quickly
## and extract only the betahats needed for predictions
## (without all the other lm object stuff), which is much faster.
## Also, now we run both the 90:10 and the 10:90 splits in this same file.



#### Setup ####

library(leaps) ## for regsubsets()

## Load in full learning and holdout sets
load("./Examples/MSD_learn.Rdata") ## learning (we will sub-split this further)
load("./Examples/MSD_holdout.Rdata")  ## holdout (used for reporting holdout MSEs)

## 90 covariates (1st column is the response Y)
p = ncol(MSD_learn) - 1

X_MSD_learn = as.matrix(MSD_learn[, -1])
Y_MSD_learn = as.numeric(MSD_learn[, 1])
X_MSD_holdout = as.matrix(MSD_holdout[, -1])
Y_MSD_holdout = as.numeric(MSD_holdout[, 1])

## Make my own 1:9 i.e. 1/10 split on the training subset:
## Train on the first 46371,
## test on the next 417344 to choose a trained model,
## finally evaluate it on the recommended holdout of 51630.
MSD_train_1090 = head(MSD_learn, 46371)
MSD_test_1090 = tail(MSD_learn, 417344)

## Matrix/vector versions for use with regsubsets()
X_MSD_train_1090 = as.matrix(MSD_train_1090[, -1])
Y_MSD_train_1090 = as.numeric(MSD_train_1090[, 1])
X_MSD_test_1090 = as.matrix(MSD_test_1090[, -1])
Y_MSD_test_1090 = as.numeric(MSD_test_1090[, 1])



## Then, reverse the roles
## so that we also have a 9:1 i.e. 9/10 split on the training subset
MSD_train_9010 = MSD_test_1090
MSD_test_9010 = MSD_train_1090

## Matrix/vector versions for use with regsubsets() etc
X_MSD_train_9010 = as.matrix(MSD_train_9010[, -1])
Y_MSD_train_9010 = as.numeric(MSD_train_9010[, 1])
X_MSD_test_9010 = as.matrix(MSD_test_9010[, -1])
Y_MSD_test_9010 = as.numeric(MSD_test_9010[, 1])



# Clean up unneeded dataframes
rm("MSD_learn", "MSD_holdout")



#### Helper functions for evaluating trained models ####

## Considerably faster if we don't have to cbind(1, test) for each new K...
## So, pre-load versions of each X matrix with an intercept column

X_MSD_holdoutint = cbind(1, X_MSD_holdout)
colnames(X_MSD_holdoutint)[1] = "(Intercept)"

X_MSD_test_1090int = cbind(1, X_MSD_test_1090)
colnames(X_MSD_test_1090int)[1] = "(Intercept)"

X_MSD_test_9010int = cbind(1, X_MSD_test_9010)
colnames(X_MSD_test_9010int)[1] = "(Intercept)"

testRMSEholdout = function(Khat, myFS) {
  Bhat = coef(myFS, Khat)
  XBhat = X_MSD_holdoutint[, names(Bhat)] %*% Bhat
  sqrt(mean((XBhat - Y_MSD_holdout)^2))
}

testRMSE1090 = function(Khat, myFS) {
  Bhat = coef(myFS, Khat)
  XBhat = X_MSD_test_1090int[, names(Bhat)] %*% Bhat
  sqrt(mean((XBhat - Y_MSD_test_1090)^2))
}

testRMSE9010 = function(Khat, myFS) {
  Bhat = coef(myFS, Khat)
  XBhat = X_MSD_test_9010int[, names(Bhat)] %*% Bhat
  sqrt(mean((XBhat - Y_MSD_test_9010)^2))
}



## For FullCV, also write a function to speed up test-error calc:
## revise all of the coef.regsubsets() output into a matrix of betahats,
## so we can just load the large 1090 test-set ONCE
## and calculate all the XBhat columns with a single matrix-multiply.
## (By default, regsubsets gives only the nonzero entries of betahat.)
betahat0 = vector("numeric", p+1)
names(betahat0) = c("(Intercept)", colnames(X_MSD_holdout))
fillbeta = function(Khat, myFS) {
  betahat0[names(coef(myFS, Khat))] = coef(myFS, Khat)
  return(betahat0)
}



## For SeqCV, going completely sequentially (train & test one Khat at a time)
## involves a lot of overhead moving between R and FORTRAN and back,
## which artificially slows us down (compared to doing it all in FORTRAN).
## We will instead go in batches: train a path of Kmod=30 Khats at a time,
## then test each of them sequentially one at a time,
## and move on to the next 30 Khats if no local min is found.
## This should be conservative for the timings we would get
## if we had built the sequential-testing directly into the FORTRAN code.




#### Actual runs for final timing, Khat, and RMSE ####

## Nr times to replicate each timing,
## to get more stable timing estimates and have Margins Of Error
nreps = 10

## Batch size for fitting SeqCV paths:
## run up to first K=30 at once, then check each of those RMSEs on test split;
## if haven't hit a local min yet, continue on to fit path for K=31:60, etc
Kmod = 30



## SeqCV, 10:90
rep_Seq1090 = replicate(nreps, {
  time_Seq1090 = system.time({
    intercept = mean(Y_MSD_train_1090)
    rmse_old = sqrt(mean((intercept - Y_MSD_test_1090)^2))  # 10.982
    fs = regsubsets(x = X_MSD_train_1090, y = Y_MSD_train_1090,
                    nvmax = Kmod, method = "forward")
    K = 1
    rmse_new = testRMSE1090(K, fs)
    while(rmse_new < rmse_old & K < p) {
      if(K %% Kmod == 0) {
        fs = regsubsets(x = X_MSD_train_1090, y = Y_MSD_train_1090,
                        nvmax = K + Kmod, method = "forward",
                        force.in = which(tail(summary(fs)$which, 1)[-1]))
      }
      K = K + 1
      rmse_old = rmse_new
      rmse_new = testRMSE1090(ifelse(K %% Kmod == 0, Kmod, K %% Kmod), fs)
    }
    Khat_Seq1090 = K - 1
  })
  return(c(time = time_Seq1090[3], Khat_Seq1090 = Khat_Seq1090))
})
mean(rep_Seq1090[1, ])  # 1.38 sec
2*sd(rep_Seq1090[1, ])/sqrt(nreps)  # 0.04 sec
(Khat_Seq1090 = rep_Seq1090[2,1])  # 23
# Fit model up to this Khat on *full* learning dataset and test on holdout data
fslearn_Seq1090 = regsubsets(x = X_MSD_learn, y = Y_MSD_learn,
                             nvmax = Khat_Seq1090, method = "forward")
testRMSEholdout(Khat_Seq1090, fslearn_Seq1090)  # 9.598



## SeqCV, 90:10
rep_Seq9010 = replicate(nreps, {
  time_Seq9010 = system.time({
    intercept = mean(Y_MSD_train_9010)
    rmse_old = sqrt(mean((intercept - Y_MSD_test_9010)^2))  # 10.561
    fs = regsubsets(x = X_MSD_train_9010, y = Y_MSD_train_9010,
                    nvmax = Kmod, method = "forward")
    K = 1
    rmse_new = testRMSE9010(K, fs)
    while(rmse_new < rmse_old & K < p) {
      if(K %% Kmod == 0) {
        fs = regsubsets(x = X_MSD_train_9010, y = Y_MSD_train_9010,
                        nvmax = K + Kmod, method = "forward",
                        force.in = which(tail(summary(fs)$which, 1)[-1]))
      }
      K = K + 1
      rmse_old = rmse_new
      rmse_new = testRMSE9010(ifelse(K %% Kmod == 0, Kmod, K %% Kmod), fs)
    }
    Khat_Seq9010 = K - 1
  })
  return(c(time = time_Seq9010[3], Khat_Seq9010 = Khat_Seq9010))
})
mean(rep_Seq9010[1, ])  # 3.39 sec
2*sd(rep_Seq9010[1, ])/sqrt(nreps)  # 0.06 sec
(Khat_Seq9010 = rep_Seq9010[2,1])   # 29
# Fit model up to this Khat on *full* learning dataset and test on holdout data
fslearn_Seq9010 = regsubsets(x = X_MSD_learn, y = Y_MSD_learn,
                             nvmax = Khat_Seq9010, method = "forward")
testRMSEholdout(Khat_Seq9010, fslearn_Seq9010)  # 9.562



## FullCV, 10:90
rep_Full1090 = replicate(nreps, {
  time_Full1090 = system.time({
    fs1090 = regsubsets(x = X_MSD_train_1090, y = Y_MSD_train_1090,
                        nvmax = p, method = "forward")
    Bhat1090 = sapply(1:p, function(x) fillbeta(x, fs1090))
    XBhat1090 = X_MSD_test_1090int %*% Bhat1090
    RMSEs1090 = sqrt(colMeans((XBhat1090 - Y_MSD_test_1090)^2))
    Khat_Full1090 = which.min(RMSEs1090)
  })
  return(c(time = time_Full1090[3], Khat_Full1090 = Khat_Full1090))
})
mean(rep_Full1090[1, ])  # 4.80 sec
2*sd(rep_Full1090[1, ])/sqrt(nreps)  # 0.06 sec
(Khat_Full1090 = rep_Full1090[2,1]) # 60
# Fit model up to this Khat on *full* learning dataset and test on holdout data
fslearn_Full1090 = regsubsets(x = X_MSD_learn, y = Y_MSD_learn,
                              nvmax = Khat_Full1090, method = "forward")
testRMSEholdout(Khat_Full1090, fslearn_Full1090)  # 9.515



## FullCV, 90:10
rep_Full9010 = replicate(nreps, {
  time_Full9010 = system.time({
    fs9010 = regsubsets(x = X_MSD_train_9010, y = Y_MSD_train_9010,
                        nvmax = p, method = "forward")
    Bhat9010 = sapply(1:p, function(x) fillbeta(x, fs9010))
    XBhat9010 = X_MSD_test_9010int %*% Bhat9010
    RMSEs9010 = sqrt(colMeans((XBhat9010 - Y_MSD_test_9010)^2))
    Khat_Full9010 = which.min(RMSEs9010)
  })
  return(c(time = time_Full9010[3], Khat_Full9010 = Khat_Full9010))
})
mean(rep_Full9010[1, ])  # 5.00 sec
2*sd(rep_Full9010[1, ])/sqrt(nreps)  # 0.06 sec
(Khat_Full9010 = rep_Full9010[2,1]) # 76
# Fit model up to this Khat on *full* learning dataset and test on holdout data
fslearn_Full9010 = regsubsets(x = X_MSD_learn, y = Y_MSD_learn,
                              nvmax = Khat_Full9010, method = "forward")
testRMSEholdout(Khat_Full9010, fslearn_Full9010)  # 9.511



#### RMSEs for Null and Full models ####

## This is only for RMSEs.
## (The timings above were done only for selection time,
##  not incl final refit to full learning set at Khat and holdout RMSE calc.)

## Null model with Khat = 0
intercept = mean(Y_MSD_learn)
sqrt(mean((intercept - Y_MSD_holdout)^2))  # 10.852

## Full model with Khat = p
fslearn_FullAll = regsubsets(x = X_MSD_learn, y = Y_MSD_learn,
                             nvmax = p, method = "forward")
testRMSEholdout(p, fslearn_FullAll)  # 9.510

