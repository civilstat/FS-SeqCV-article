SimData <-
structure(list(nrSims = c(1000L, 1000L, 1000L, 1000L, 1000L, 
1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 
1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 
1000L, 1000L, 1000L, 1000L, 1000L, 1000L, 1000L), n = c(6250L, 
6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 
6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 
6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 6250L, 
6250L, 6250L), p = c(11L, 11L, 11L, 11L, 11L, 50L, 50L, 50L, 
50L, 50L, 250L, 250L, 250L, 250L, 250L, 11L, 11L, 11L, 11L, 11L, 
50L, 50L, 50L, 50L, 50L, 250L, 250L, 250L, 250L, 250L), K = c(10L, 
10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 
10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 10L, 
10L, 10L, 10L), betamin = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 
0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), betamax = c(2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), beta.fun = c("beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform", "beta.Uniform", "beta.Uniform", "beta.Uniform", 
"beta.Uniform"), Sigma.fun = c("Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", 
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr"
), epsilon.fun = c("epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2"), doNormalize = c(FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
), rescaleByN = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE), mu = c(0.05, 0.05, 0.05, 
0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 
0.25, 0.25, 0.25, 0.25, 0.25), procedure = structure(c(1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "FS_FixedK", class = "factor"), 
    outcome = structure(c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
    5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
    5L, 1L, 2L, 3L, 4L, 5L), .Label = c("foundTrueModel", "nrFalsePos", 
    "nrFalseNegs", "estimationSSE", "elapsed"), class = "factor"), 
    mean = c(0.989, 0.011, 0.011, 0.0268328270174824, 0.030665, 
    0.918, 0.094, 0.094, 0.0354482242246255, 0.15993, 0.895, 
    0.131, 0.131, 0.0557351665451622, 0.607306, 0.988, 0.012, 
    0.012, 0.0321376756810961, 0.0310879999999999, 0.913, 0.104, 
    0.104, 0.125105324005605, 0.160993, 0.833, 0.207, 0.207, 
    0.0700738255537239, 0.614469999999998), sd = c(0.104354635210372, 
    0.104354635210372, 0.104354635210372, 0.0581967783877949, 
    0.034168908225012, 0.274502006097135, 0.336566898665447, 
    0.336566898665447, 0.0862205525981664, 0.0395716187818114, 
    0.306706812883361, 0.44725724034278, 0.44725724034278, 0.336145960894485, 
    0.056751291378019, 0.108939744206914, 0.108939744206914, 
    0.108939744206914, 0.0653530921772707, 0.0346859246578403, 
    0.281976081451088, 0.428214160633867, 0.428214160633867, 
    2.60764395105608, 0.0399342922649032, 0.373162498451077, 
    0.548134519485365, 0.548134519485365, 0.470822129842142, 
    0.057250298660329)), .Names = c("nrSims", "n", "p", "K", 
"betamin", "betamax", "beta.fun", "Sigma.fun", "epsilon.fun", 
"doNormalize", "rescaleByN", "mu", "procedure", "outcome", "mean", 
"sd"), row.names = c(NA, -30L), class = "data.frame")
mySeed <-
170924
indir <-
"~/CV/"
outdir <-
"~/CV/"
simSettingsList <-
NULL
whichProcedures <-
"FS_FixedK"
useParallel <-
TRUE
nCores <-
6
