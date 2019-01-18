## Recreate parts of Table 6.1 from thesis
## as well as numerical results mentioned inline in Section 6.2.2
## using a 10:90 train:test split within the original training data


#### Setup ####

library(MASS)  ## for addterm()

## Load in full training and holdout sets
load("./Examples/MSD_train.Rdata")
load("./Examples/MSD_test.Rdata")

## 90 covariates (1st column is the response Y)
p = ncol(MSD_test) - 1

## Make my own 1:9 i.e. 1/10 split on the training subset:
## Train on the first 46371,
## test on the next 417344 to choose a trained model,
## finally evaluate it on the recommended holdout of 51630.
MSD_train_1090 = head(MSD_train, 46371)
MSD_test_1090 = tail(MSD_train, 417344)

## Function for evaluating trained models on holdout set.
testRMSE = function(mylm, testdata = MSD_test) {
  sqrt(mean((predict(mylm, newdata = testdata) - testdata$Y)^2))
}



#### Baselines: null and full models ####

## Train null and full models on subsplit of training data,
## so we can use them to specify scope for stepAIC()
newTime = system.time({
  lm0_1090 = lm(Y ~ 1, data = MSD_train_1090)
})
lm1_1090 = lm(Y ~ ., data = MSD_train_1090)



#### FS+CV, 10:90 train:test split ####

## Take one step on training data,
## evaluate on testing data;
## continue up to full p-variable model.
## At end, choose global minimizer of testMSE.
## Also keep track of timings,
## so we can also do SeqCV (find 1st local minimizer)
## and see how much faster it would have been.
system.time({
  ## Initialize
  MSEs = rep(0, p+1)
  Vars = rep(0, p+1)
  Times = rep(0, p+1)
  jj = 0
  lmOld = lm0_1090
  Times[jj + 1] = newTime["elapsed"]
  MSEs[jj + 1] = testRMSE(lmOld, MSD_test_1090)^2
  Vars[jj + 1] = "1"
  print(c(jj, "1", MSEs[jj + 1]))
  
  jj = jj + 1
  newTime = system.time({
    addTerm = addterm(lmOld, scope = lm1_1090, sorted = TRUE, k = 0)
    newVar = rownames(addTerm)[1]
    Vars[jj + 1] = newVar
    lmNew = update(lmOld, formula(paste(". ~ . +", newVar)))
    MSEs[jj + 1] = testRMSE(lmNew, MSD_test_1090)^2
  })
  Times[jj + 1] = newTime["elapsed"]
  print(c(jj, newVar, MSEs[jj + 1], MSEs[jj + 1] - MSEs[jj]))
  
  ## Repeat as needed
  while(jj < p) {
    jj = jj + 1
    newTime = system.time({
      lmOld = lmNew
      addTerm = addterm(lmOld, scope = lm1_1090, sorted = TRUE, k = 0)
      newVar = ifelse(rownames(addTerm)[1] != "<none>",
	                  rownames(addTerm)[1],
	                  rownames(addTerm)[2])
      Vars[jj + 1] = newVar
      (lmNew = update(lmOld, formula(paste(". ~ . +", newVar))))
      MSEs[jj + 1] = testRMSE(lmNew, MSD_test_1090)^2
    })
    Times[jj + 1] = newTime["elapsed"]
    print(c(jj, newVar, MSEs[jj + 1], MSEs[jj + 1] - MSEs[jj]))
  }
})
## Runs in 2001 sec = 33 min
sum(Times)

## Store all MSEs for later plotting
dump("MSEs", ""); dump("MSEs", "./Examples/MSD_FsCV_1090_results.R")
MSEs <-
  c(120.612352933778, 114.341693747119, 109.637849181229, 105.981385676691, 
    103.356527120918, 102.062654120613, 99.3791626898993, 98.8183165421964, 
    98.2382468177767, 97.6092225882649, 96.9324530578183, 96.407428704165, 
    95.9098747786536, 95.5797434836676, 95.3201109890231, 95.125570380643, 
    94.840358894258, 94.7841723591399, 94.6243134823434, 94.4191102048902, 
    94.2391716258156, 94.0352143335164, 94.0016020312707, 93.9471875664691, 
    93.9704903915423, 93.759575898655, 93.6789141527421, 93.5293673279041, 
    93.3437450224954, 93.2599062975889, 93.1578440795041, 93.1558463707902, 
    93.0261704785298, 93.0260785662193, 92.9943735214749, 92.926045417545, 
    92.9006440089861, 92.8787929262866, 92.8594362866407, 92.8023000088047, 
    92.7744462927719, 92.7858738912922, 92.7501262600225, 92.7060195309718, 
    92.6404379411498, 92.6042906024894, 92.5404090890593, 92.5455548956212, 
    92.5389392579048, 92.5118095610024, 92.4846762102917, 92.4630543032223, 
    92.4608092035434, 92.4884251948643, 92.4695096087756, 92.4732311696567, 
    92.471195045236, 92.4208114243418, 92.4291722516833, 92.4275753576215, 
    92.4133482374419, 92.4206997496024, 92.4622848900912, 92.4769264971545, 
    92.4925573094741, 92.4962930825615, 92.5257057764136, 92.5233616734751, 
    92.5197238137441, 92.5038076606689, 92.5016385886808, 92.5003265223934, 
    92.5014041916341, 92.4904080205219, 92.4949565996541, 92.497027177226, 
    92.4992959934878, 92.498914561959, 92.4898757432342, 92.4909811155766, 
    92.491398122265, 92.5042517405085, 92.5026365798013, 92.4952616802373, 
    92.4924670359203, 92.4888021094335, 92.4903785557321, 92.4910034320319, 
    92.4912239355907, 92.491578540134, 92.4916796931566)

## Store all Times for later plotting
dump("Times", ""); dump("Times", "./Examples/MSD_FsCV_1090_results.R", append = TRUE)
Times <-
  c(0.0500000000000007, 1.339, 1.578, 1.512, 1.443, 1.606, 2.271, 
    2.895, 4.45800000000001, 5.072, 6.176, 7.708, 8.81, 8.568, 10.174, 
    11.232, 11.096, 11.64, 12.088, 12.947, 9.738, 13.692, 15.59, 
    15.907, 16.588, 17.492, 16.951, 16.091, 15.941, 19.706, 16.93, 
    22.591, 23.255, 24.852, 25.059, 20.439, 27.185, 25.379, 26.378, 
    30.293, 29.712, 26.7570000000001, 31.798, 32.5790000000001, 33.707, 
    27.679, 34.799, 35.745, 29.45, 36.011, 36.207, 36.91, 37.605, 
    33.1569999999999, 36.357, 34.1899999999998, 38.021, 38.057, 38.251, 
    37.557, 39.0250000000001, 32.7760000000001, 38.646, 38.1610000000001, 
    38.8719999999998, 31.623, 31.8780000000002, 36.6510000000001, 
    38.0159999999998, 36.7719999999999, 31.1859999999999, 29.298, 
    35.1470000000002, 30.846, 28.777, 27.8150000000001, 24.4770000000001, 
    25.4360000000001, 22.5639999999999, 23.925, 16.8190000000002, 
    15.0820000000001, 17.6309999999999, 18.934, 18.3320000000001, 
    11.7559999999999, 11.8489999999999, 11.598, 11.2530000000002, 
    9.83999999999992, 8.31500000000005)



## Find the global minimizer
(kFsCv = which.min(MSEs) - 1)
formulaFsCv = as.formula(paste("Y ~ ",
                               paste(Vars[1:(kFsCv+1)], collapse = " + ")))
## Selects 60 variables out of 90
## (so ~34 sec/variable)
dump("formulaFsCv", ""); dump("formulaFsCv", "./Examples/MSD_FsCV_1090_results.R", append = TRUE)
formulaFsCv <-
  quote(Y ~ 1 + Mean1 + Var2 + Mean3 + Mean2 + Mean6 + Var11 + 
          Cov17 + Cov33 + Cov3 + Var8 + Var1 + Cov12 + Cov54 + Var4 + 
          Cov61 + Cov24 + Cov35 + Var7 + Cov40 + Mean11 + Cov1 + Var6 + 
          Cov64 + Mean8 + Mean9 + Var3 + Cov16 + Cov34 + Cov14 + Mean5 + 
          Cov63 + Cov45 + Cov29 + Cov65 + Cov15 + Cov50 + Cov26 + Cov52 + 
          Cov47 + Cov48 + Cov18 + Cov9 + Cov11 + Cov51 + Cov23 + Cov41 + 
          Cov44 + Cov49 + Cov22 + Cov39 + Cov46 + Cov10 + Cov62 + Var5 + 
          Cov38 + Cov60 + Var12 + Cov32 + Cov28 + Mean10)


## Refit to full training set
newTime = system.time({ lmFsCv = lm(formulaFsCv, data = MSD_train) })
newTime + sum(Times)
## Whole process took 2013 sec = 34 min
## Evaluate on holdout set
testRMSE(lmFsCv)
## RMSE = 9.52



## Find the first local minimizer
(kFsSeqCv = which.max(sign(diff(MSEs))) - 1)
formulaFsSeqCv = as.formula(paste("Y ~ ",
                                  paste(Vars[1:(kFsSeqCv+1)], collapse = " + ")))
## Selects 23 variables out of 90
sum(Times[1:(kFsSeqCv+1)])
## in 178 sec = 3 min
## (so ~8 sec/variable)
dump("formulaFsSeqCv", ""); dump("formulaFsSeqCv", "./Examples/MSD_FsCV_1090_results.R", append = TRUE)
formulaFsSeqCv <-
  quote(Y ~ 1 + Mean1 + Var2 + Mean3 + Mean2 + Mean6 + Var11 + 
          Cov17 + Cov33 + Cov3 + Var8 + Var1 + Cov12 + Cov54 + Var4 + 
          Cov61 + Cov24 + Cov35 + Var7 + Cov40 + Mean11 + Cov1 + Var6 + 
          Cov64)

## Refit to full training set
newTime = system.time({ lmFsSeqCv = lm(formulaFsSeqCv, data = MSD_train) })
newTime + sum(Times[1:(kFsSeqCv+1)])
## Whole process took 180 sec = 3 min
## Evaluate on holdout set
testRMSE(lmFsSeqCv)
## RMSE = 9.61

