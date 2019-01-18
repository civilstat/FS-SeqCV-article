# Recreate Figure 4.1 in the thesis.

# Simulating small example datasets
# to show that, with high SNR,
# CV with large p helps account for
# post-selection effect of greedy FS



#### Define functions ####

L2 = function(x) sum(x^2)

genData = function(n = 30, p = 2, B0 = 10) {
  X = data.frame(matrix(rnorm(p*n), nrow = n))
  epsilon = rnorm(n)
  Y = B0 + epsilon
  return(list(X = X, Y = Y))
}

# FS with CV (really just a single sample split),
# training and then reporting test RSS for comparing just these 3 models:
# 0 (intercept-only model),
# 1Fixed (intercept and X_1, no matter what the other Xs are), and 
# 1Greedy (intercept and X_i, where X_i is greedily chosen for max |corr| with Y)
#
# Compute training fits and test RSSs for a given Train/Test split:
FS1CV_OnSplits = function(Y, X, iiTrain, iiTest) {
  dataset = data.frame(Y, X)
  
  # Greedily choose X_i based on training data
  # whichX = which.max(abs(cor(cbind(Y[iiTrain], X[iiTrain, ])))[1, -1])
  whichX = which.max(abs(cor(Y[iiTrain], X[iiTrain, ])))
  
  # Fit all 3 models on training data
  trainFit0 = lm(Y ~ 1, dataset[iiTrain, ])
  trainFit1Fixed = lm(Y ~ ., dataset[iiTrain, c(1, 2)])
  trainFit1Greedy = lm(Y ~ ., dataset[iiTrain, c(1, whichX+1)])
  
  # Compute test error
  testRss0 = L2(Y[iiTest] - predict(trainFit0, newdata = X[iiTest, ]))
  testRss1Fixed = L2(Y[iiTest] - predict(trainFit1Fixed, newdata = X[iiTest, ]))
  testRss1Greedy = L2(Y[iiTest] - predict(trainFit1Greedy, newdata = X[iiTest, ]))
  
  # Return test errors to be compared
  c(testRss1Fixed = testRss1Fixed,
    testRss0 = testRss0,
    testRss1Greedy = testRss1Greedy)
}



#### Setup ####

# Sometimes, we'll find
# testRss1Fixed < testRss0,
# so if we only had X1 in the dataset,
# we'd declare it to be significant
# and choose model size Khat >= 1.
# (And frequency of this stays constant with p.)
# 
# But almost always, we'll find
# testRss0 < testRss1Greedy,
# so if we had several X variables in the dataset,
# we'd declare the first-chosen to be spurious
# and choose model size Khat = 0.
# (And frequency of this rises with p.)
#
# So higher p makes CV work better, in this particular setting.

# For example:
set.seed(2001)
n = 30; p = 2; B0 = 10; iiTrain = 1:(n/2); iiTest = (n/2 + 1):n
dataset = genData(n, p, B0)
(OneRepP2 = FS1CV_OnSplits(dataset$Y, dataset$X, iiTrain, iiTest))

#### Show an example in detail: ####
oldPar = par(no.readonly = TRUE)

pdf("./SmallSims/HighP_17-10-02.pdf",
    width = 6.5, height = 5.5, pointsize = 10)
attach(dataset)
par(mfrow = c(2, 2))
for(cc in 1:p){
  plot(X[iiTrain, cc], Y[iiTrain], las = 1,
       xlab = ifelse(cc==1,
                     expression(X[list(1,train)]),
                     expression(X[list(FS,train)])),
       ylab = expression(Y[train]),
       xlim = c(-3, 3),
       ylim = range(Y))
  if(cc == 1){
    title(main = "Best fit with low p is likely flat")
  } else {
    title(main = "Best fit with high p may be steep")
  }
  abline(lm(Y[iiTrain] ~ X[iiTrain, cc]))
  abline(lm(Y[iiTrain] ~ 1), col = "grey35", lty = 2)
  text(-1.3, 9, paste0("Corr(Y, X", 
                       #ifelse(cc == 1, 1, "_FS"),
                       ") = ",
                      signif(cor(Y[iiTrain], X[iiTrain, cc]), 2)))
  plot(X[iiTest, cc], Y[iiTest], las = 1,
       xlab = ifelse(cc==1,
                     expression(X[list(1,test)]),
                     expression(X[list(FS,test)])),
       ylab = expression(Y[test]),
       xlim = c(-3, 3),
       ylim = range(Y))
  if(cc == 1){
    title(main = "Test data may accept a flat fit")
  } else {
    title(main = "Test data will reject a steep fit")
  }
  abline(lm(Y[iiTrain] ~ X[iiTrain, cc]))
  abline(lm(Y[iiTrain] ~ 1), col = "grey35", lty = 2)
  legend('topright',
         c(paste0("RSS = ", signif(OneRepP2[-2][cc], 3)),
           paste0("RSS = ", signif(OneRepP2[2], 3))),
         lty = 1:2, col = c("black", "grey35"), bty = "n")
}
par(oldPar)
detach(dataset)
dev.off()

