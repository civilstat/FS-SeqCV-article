# Create the simulated data used for Figure S2 in supplement


#### SIMULATE ####

## * P(overfit) and P(underfit) vs n_c/n, at various SNRs,
##   at least for a correct fixed path;
##   note that the theory (from Zhang) assumes multiple-folds CV-a,
##   so we're using 20 reps of multifold / MC-CV

## k=5, p=6: ran for 4223 sec = 70 min
## k=10, p=11: ran for 8637 sec = 144 min

library(MASS)
for(p in c(6, 11)) {
  k = p-1
  sigma = 1
  ## "SNR" = beta^2/(sigma^2/n)
  
  n = 500 # want at least 300, so have nc=15 for 5:95 ratio on k=10, p=11
  M = 600 # nr datasets to generate, so MOEs are +/- 0.02
  B = 20 # nr reps of CV per dataset
  ## last year, sims showed it's stable by 20 folds / CV rep,
  ## but maybe that's not as true for very high ratios, nearly LOO?
  
  ratios = 1:19/20 ## Try at train:test ratios from 0.05 to 0.95

  ## Orig. tried SNRs around 20, 40, 80, but 40 and 80 looked identical,
  ## and all were pretty low except for the smallest ratios,
  ## so let's try some lower SNRs too
  SNRs = c(1, 5, 10, 20, 40)
  betamins = sqrt(SNRs / n)

  nRs = length(ratios)
  nSNRs = length(SNRs)
  Punders = matrix(NA, nRs, nSNRs,
                   dimnames = list(ratios, SNRs))
  Povers = Punders
  
  runtime = system.time({
    for(ii_SNRs in 1:nSNRs) {
      for(ii_rats in 1:nRs) {
        beta = c(rep(betamins[ii_SNRs], k), rep(0, p-k))
        nc = ceiling(n*ratios[ii_rats])
        
        under = 0; over = 0
        for(mm in 1:M){
          X = mvrnorm(n, mu = rep(0, p), Sigma = diag(p))
          eps = rnorm(n, 0, sigma)
          Y = drop(X %*% beta + eps)
          
          MSEs = rep(0, p)
          datanow = data.frame(Y = Y, X = X)
          for(bb in 1:B){
            ii_train = sample(n, nc)
            train = datanow[ii_train, ]
            test = datanow[-ii_train, ]
            for(pp in 1:p){
              Xtrain = cbind(1, as.matrix(train[, 2:(pp+1)]))
              betahat = solve(t(Xtrain)%*%Xtrain) %*% t(Xtrain)%*%as.matrix(train[, 1, drop=FALSE])
              prednow = cbind(1, as.matrix(test[, 2:(pp+1)])) %*% betahat
              MSEnow = mean((prednow - test$Y)^2)
              MSEs[pp] = MSEs[pp] + MSEnow
            }
          }
          MSEs = MSEs/B
          under = under + any(MSEs[k] > MSEs[1:(k-1)])
          over = over + (MSEs[k] > MSEs[k+1])
        }
        Punders[ii_rats, ii_SNRs] = under/M
        Povers[ii_rats, ii_SNRs] = over/M
      }
    }
  })
  dump(c("Punders", "Povers"), paste0("./SmallSims/FigS2_ProbsAtK", k, "_20170904.R"))
  print(runtime)
  print(k)
  print(Punders)
  print(Povers)
}

