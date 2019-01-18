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
50L, 50L, 50L, 50L, 50L, 250L, 250L, 250L, 250L, 250L), K = c(5L, 
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
"Sigma.ConstCorr", "Sigma.ConstCorr", "Sigma.ConstCorr"), epsilon.fun = c("epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2", 
"epsilon.t2", "epsilon.t2", "epsilon.t2", "epsilon.t2"), doNormalize = c(FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE), rescaleByN = c(FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), mu = c(0.1, 
0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
0.1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
0.5, 0.5, 0.5), procedure = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "FS_FixedK", class = "factor"), 
    outcome = structure(c(1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
    5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 5L, 1L, 2L, 3L, 4L, 
    5L, 1L, 2L, 3L, 4L, 5L), .Label = c("foundTrueModel", "nrFalsePos", 
    "nrFalseNegs", "estimationSSE", "elapsed"), class = "factor"), 
    mean = c(0.973, 0.028, 0.028, 0.0179455341459552, 0.029535, 
    0.93, 0.074, 0.074, 0.0220223348026431, 0.148816, 0.874, 
    0.135, 0.135, 0.0368565006011572, 0.598063999999999, 0.929, 
    0.076, 0.076, 0.0294639736452097, 0.0310320000000001, 0.803, 
    0.208, 0.208, 0.0619442115549897, 0.15235, 0.709, 0.315, 
    0.315, 0.0654382414660157, 0.600452999999999), sd = c(0.162164414398774, 
    0.171012412547292, 0.171012412547292, 0.0696564459559341, 
    0.0327423401093422, 0.255274685711618, 0.276768135088924, 
    0.276768135088924, 0.0610216670865698, 0.0372785438802779, 
    0.332015412645609, 0.372711569331989, 0.372711569331989, 
    0.214451428171302, 0.055195293498852, 0.256953351846254, 
    0.28338014098434, 0.28338014098434, 0.0909646119719137, 0.0347833992137178, 
    0.397931337480914, 0.441509831069401, 0.441509831069401, 
    0.627692939114566, 0.0383218745113815, 0.454450795494436, 
    0.527307361083698, 0.527307361083698, 0.256898668263704, 
    0.0552628788317873)), .Names = c("nrSims", "n", "p", "K", 
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
