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

## Make my own 9:1 i.e. 9/10 split on the training subset,
## simply swapping roles of subsets from MSD_FsCv_1090.R:
## Train on the LAST 417344,
## test on the FIRST 46371 to choose a trained model,
## finally evaluate it on the recommended holdout of 51630.
MSD_train_9010 = tail(MSD_train, 417344)
MSD_test_9010 = head(MSD_train, 46371)

## Function for evaluating trained models on holdout set.
testRMSE = function(mylm, testdata = MSD_test) {
  sqrt(mean((predict(mylm, newdata = testdata) - testdata$Y)^2))
}




#### Baselines: null and full models ####

## Train null and full models on subsplit of training data
## so we can use them to specify scope for stepAIC()
newTime = system.time({
  lm0_9010 = lm(Y ~ 1, data = MSD_train_9010)
})
lm1_9010 = lm(Y ~ ., data = MSD_train_9010)



#### FS+CV, 90:10 train:test split ####

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
  lmOld = lm0_9010
  Times[jj + 1] = newTime["elapsed"]
  MSEs[jj + 1] = testRMSE(lmOld, MSD_test_9010)^2
  Vars[jj + 1] = "1"
  print(c(jj, "1", MSEs[jj + 1]))
  
  jj = jj + 1
  newTime = system.time({
    addTerm = addterm(lmOld, scope = lm1_9010, sorted = TRUE, k = 0)
    newVar = rownames(addTerm)[1]
    Vars[jj + 1] = newVar
    lmNew = update(lmOld, formula(paste(". ~ . +", newVar)))
    MSEs[jj + 1] = testRMSE(lmNew, MSD_test_9010)^2
  })
  Times[jj + 1] = newTime["elapsed"]
  print(c(jj, newVar, MSEs[jj + 1], MSEs[jj + 1] - MSEs[jj]))
  
  ## Repeat as needed
  while(jj < p) {
    jj = jj + 1
    newTime = system.time({
      lmOld = lmNew
      addTerm = addterm(lmOld, scope = lm1_9010, sorted = TRUE, k = 0)
      newVar = ifelse(rownames(addTerm)[1] != "<none>",
	                  rownames(addTerm)[1],
	                  rownames(addTerm)[2])
      Vars[jj + 1] = newVar
      (lmNew = update(lmOld, formula(paste(". ~ . +", newVar))))
      MSEs[jj + 1] = testRMSE(lmNew, MSD_test_9010)^2
    })
    Times[jj + 1] = newTime["elapsed"]
    print(c(jj, newVar, MSEs[jj + 1], MSEs[jj + 1] - MSEs[jj]))
  }
})
## Runs in 16320 sec = 272 min
sum(Times)

## Store all MSEs for later plotting
dump("MSEs", ""); dump("MSEs", "./Examples/MSD_FsCV_9010_results.R")
MSEs <-
  c(111.535377245901, 106.569142781289, 102.350277248978, 99.0694900606472, 
    97.4696400593366, 96.5869633427118, 94.3816329337353, 93.714312942435, 
    93.1295123583896, 92.5462459278345, 91.9307240381637, 91.3163386926909, 
    91.2475019972275, 91.0123181541998, 90.9683107652144, 90.8241980611048, 
    90.79794987503, 90.4270781445618, 90.2205254433965, 90.068392719242, 
    89.7952477006051, 89.6175098782025, 89.3820575960944, 89.2969965295444, 
    89.09088026604, 89.0094859678329, 88.9138766844929, 88.8095365855988, 
    88.7490974805363, 88.6872167260627, 88.7105146061802, 88.6227410481742, 
    88.5800585426573, 88.4734982218131, 88.4994891537648, 88.5257282927206, 
    88.4809537391096, 88.413893097602, 88.3492333551454, 88.3532362131067, 
    88.3134174697004, 88.2842339709034, 88.2433495383411, 88.1894312418918, 
    88.1554669185524, 88.0962059971146, 88.1090695729836, 88.1104005971853, 
    88.152989459283, 88.1385541993411, 88.1406441333862, 88.1007618032768, 
    88.0885963214717, 88.0718712316456, 88.0304537609237, 88.0030938072487, 
    88.0167915656645, 88.0085903088807, 87.9728419751828, 87.9802947673344, 
    87.9959963553795, 87.9856988046771, 87.9752638049245, 87.9793545272116, 
    87.9727911711825, 87.9628891240382, 87.9541332034122, 87.9517717959082, 
    87.9528368539397, 87.9419944958304, 87.9290974602133, 87.9213666867734, 
    87.9289131080322, 87.924873987217, 87.9209714655132, 87.916257628467, 
    87.9133644078147, 87.9198331676613, 87.9155924138457, 87.9180137481436, 
    87.92082942105, 87.9216925033874, 87.9199240530881, 87.9184519942761, 
    87.9208471716391, 87.9185059878767, 87.9198780119541, 87.9217502254488, 
    87.9208110181652, 87.9206648787287, 87.9209054933615)

## Store all Times for later plotting
dump("Times", ""); dump("Times", "./Examples/MSD_FsCV_9010_results.R", append = TRUE)
Times <-
  c(0.491999999999997, 13.218, 15.973, 19.121, 18.649, 26.429, 
    29.944, 33.048, 35.985, 43.055, 53.503, 53.733, 60.025, 73.323, 
    74.998, 81.687, 86.467, 101.16, 99.349, 107.552, 120.094, 122.853, 
    121.504, 132.803, 120.183, 104.821, 97.155, 142.926, 138.492, 
    152.769, 174.428, 162.139, 232.453, 182.708, 197.876, 229.491, 
    319.013, 311.89, 323.322999999999, 206.864, 220.845, 232.5, 218.074000000001, 
    239.925999999999, 233.032, 243.915, 243.28, 265.396000000001, 
    266.496, 274.569, 259.497, 272.13, 274.645, 284.703, 274.644999999999, 
    274.832, 273.690999999999, 290.414000000001, 268.079, 262.35, 
    282.290999999999, 283.82, 271.537, 289.346, 284.814, 272.068000000001, 
    276.515000000001, 267.512000000001, 264.960000000001, 277.364000000001, 
    266.556, 265.085000000001, 245.348, 226.412, 247.381000000001, 
    236.690999999999, 225.534, 210.678, 312.843000000001, 286.241, 
    227.815000000001, 147.690000000001, 143.766, 127.503000000001, 
    111.365, 109.389000000001, 103.525, 86.8990000000013, 77.4300000000003, 
    58.2540000000008, 39.0080000000016)



## Find the global minimizer
(kFsCv = which.min(MSEs) - 1)
formulaFsCv = as.formula(paste("Y ~ ",
                               paste(Vars[1:(kFsCv+1)], collapse = " + ")))
## Selects 76 variables out of 90
## (so ~215 sec/variable)
dump("formulaFsCv", ""); dump("formulaFsCv", "./Examples/MSD_FsCV_9010_results.R", append = TRUE)
formulaFsCv <-
  quote(Y ~ 1 + Mean1 + Mean3 + Var2 + Mean2 + Mean6 + Var11 + 
          Cov14 + Var8 + Var1 + Cov33 + Cov3 + Mean9 + Cov12 + Cov34 + 
          Cov16 + Cov24 + Cov17 + Cov54 + Mean11 + Mean8 + Cov61 + 
          Var6 + Cov37 + Cov35 + Var4 + Var3 + Var7 + Cov1 + Mean5 + 
          Cov45 + Cov40 + Cov15 + Cov64 + Cov41 + Var12 + Cov9 + Cov52 + 
          Cov47 + Cov59 + Cov65 + Cov51 + Cov26 + Cov50 + Cov23 + Cov63 + 
          Var9 + Cov21 + Cov5 + Cov22 + Var5 + Cov28 + Mean10 + Cov39 + 
          Cov46 + Cov49 + Cov42 + Cov10 + Cov29 + Cov4 + Cov13 + Cov19 + 
          Cov8 + Cov20 + Cov58 + Cov60 + Cov38 + Cov7 + Mean7 + Cov18 + 
          Cov48 + Cov11 + Mean4 + Cov25 + Cov32 + Cov53 + Var10)

## Refit to full training set
newTime = system.time({ lmFsCv = lm(formulaFsCv, data = MSD_train) })
newTime + sum(Times)
## Whole process took 16336 sec = 272 min
## Evaluate on holdout set
testRMSE(lmFsCv)
## RMSE = 9.51



## Find the first local minimizer
(kFsSeqCv = which.max(sign(diff(MSEs))) - 1)
formulaFsSeqCv = as.formula(paste("Y ~ ",
                                  paste(Vars[1:(kFsSeqCv+1)], collapse = " + ")))
## Selects 29 variables out of 90
sum(Times[1:(kFsSeqCv+1)])
## in 2285 sec = 38 min
## (so ~79 sec/variable)
dump("formulaFsSeqCv", ""); dump("formulaFsSeqCv", "./Examples/MSD_FsCV_9010_results.R", append = TRUE)
formulaFsSeqCv <-
  quote(Y ~ 1 + Mean1 + Mean3 + Var2 + Mean2 + Mean6 + Var11 + 
          Cov14 + Var8 + Var1 + Cov33 + Cov3 + Mean9 + Cov12 + Cov34 + 
          Cov16 + Cov24 + Cov17 + Cov54 + Mean11 + Mean8 + Cov61 + 
          Var6 + Cov37 + Cov35 + Var4 + Var3 + Var7 + Cov1 + Mean5)

## Refit to full training set
newTime = system.time({ lmFsSeqCv = lm(formulaFsSeqCv, data = MSD_train) })
newTime + sum(Times[1:(kFsSeqCv+1)])
## Whole process took 2285 sec = 38 min
## Evaluate on holdout set
testRMSE(lmFsSeqCv)
## RMSE = 9.57

