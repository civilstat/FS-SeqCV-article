SimData <-
structure(list(nrSims = c(1000L, 1000L, 1000L, 1000L, 1000L, 
1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 
1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 
1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L), p = c(11L, 
11L, 11L, 11L, 11L, 50L, 50L, 50L, 50L, 50L, 250L, 250L, 250L, 
250L, 250L, 11L, 11L, 11L, 11L, 11L, 50L, 50L, 50L, 50L, 50L, 
250L, 250L, 250L, 250L, 250L), n = c(6250L, 6250L, 6250L, 6250L, 
6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 
6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 
6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L), K = c(5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L), betamin = c(0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2), betamax = c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L), beta.fun = c("beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform"
), Sigma.fun = c("Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr"), epsilon.fun = c("epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal", "epsilon.Normal", "epsilon.Normal", "epsilon.Normal", 
"epsilon.Normal"), doNormalize = c(FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), rescaleByN = c(FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE), mu = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), procedure = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "FS_K5", class = "factor"), 
    outcome = structure(c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
    5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
    5L, 1L, 2L, 3L, 4L, 5L), .Label = c("foundTrueModel", "nrFalsePos", 
    "nrFalseNegs", "estimationSSE", "elapsed"), class = "factor"), 
    mean = c(1, 0, 0, 0.00100171982655864, 0.025231, 1, 0, 0, 
    0.000988540746055899, 0.108086, 1, 0, 0, 0.000991821248151434, 
    0.509877000000001, 1, 0, 0, 0.00148254886927966, 0.0250909999999999, 
    1, 0, 0, 0.00146388176443612, 0.107797, 1, 0, 0, 0.00149323724077577, 
    0.513163999999998), sd = c(0, 0, 0, 0.000576707874987201, 
    0.018054573010919, 0, 0, 0, 0.000571695314111554, 0.015058132399184, 
    0, 0, 0, 0.00057043408890184, 0.0266040531995892, 0, 0, 0, 
    0.000993876568116858, 0.0178285887713407, 0, 0, 0, 0.000917769849111677, 
    0.014785682514088, 0, 0, 0, 0.000921430582353427, 0.0284947548691884
    )), .Names = c("nrSims", "p", "n", "K", "betamin", "betamax", 
"beta.fun", "Sigma.fun", "epsilon.fun", "doNormalize", "rescaleByN", 
"mu", "procedure", "outcome", "mean", "sd"), row.names = c(NA, 
-30L), class = "data.frame")
mySeed <-
1709105
indir <-
"~/CV/"
outdir <-
"~/CV/"
simSettingsList <-
NULL
whichProcedures <-
"FS_K5"
useParallel <-
TRUE
nCores <-
6
