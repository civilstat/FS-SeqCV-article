# Create the simulated data used for Figure S5 in supplement:
# specifically, the simulations at the largest n=6250.



# 5 sims at nCores=6 took 430+23=453 sec,
# so 100 sims should take 9060 sec = 151 min,
# so let's do 1000 sims in 25.5 hrs.
# Starting on hydra2 around 3:35pm Sept 24,
# expect to finish by 5pm Sept 25.

# Repeating 17-09-14 "small n" with "bad" t2 noise,
# but now adding the largest n=6250.
# 
# Running
# (SeqCVa vs FullCVa) x (1S vs VF),
# at 3 train:test ratios (4:1, 1:1, and 1:4).
# Also add KnownK,
# all at several combos of n, mu, k, p.


mySeed = 170924 # actually running on Sept 24 not 14
mySuffix = "17-09-14_LargeN"

# Keep doNormalize = FALSE.

# Use data settings:
# Sigma = ConstCorr
# K = 3,5,10
# p = 11,50,250
# n = 6250
# mu = 1/(2K) * [1, 5]
# beta_min = 1/5
# epsilon ~ t2
#
# and use either 1S or 5F CV-a variants, with Seq vs Full CV,
# at train/test ratios 4:1, 1:1, and 1:4.
# Also KnownK with K=3,5,10
# (run separately, so we can tailor the alg to K)



### SPECIFY FILE PATH, LOAD FUNCTIONS
indir = "./Functions/"
outdir = "./MainSims/"

# Flag whether or not to use parallel processing,
# and if so, how many cores to use
useParallel = TRUE; nCores = 6        # on department servers
# useParallel = TRUE; nCores = 2        # on desktop
# useParallel = FALSE; nCores = NULL    # for debugging w/o running in parallel



#### FIRST, all the estimated-K methods ####

### SPECIFY SIMULATION SETTINGS
# Specify the simulation settings for which we want all cross-combinations,
# or any constants we want reported such as nrSims
# (but define whichProcedures separately, so each new dataset runs all procedures)
simSettingsList = list(nrSims = 1000L,
                       n = c(6250L),
                       p = c(11L, 50L, 250L),
                       K = c(3L, 5L, 10L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.t2",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("SeqCVa_4to1S", "FullCVa_4to1S", "SeqCVa_5F", "FullCVa_5F",
                    "SeqCVa_1to1S", "FullCVa_1to1S", "SeqCVa_2F", "FullCVa_2F",
                    "SeqCVa_1to4S", "FullCVa_1to4S", "SeqCVa_i5F", "FullCVa_i5F")

# Specify procedure variants with non-default options.
# 4:1 ratios:
SeqCVa_4to1S = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1, trainProportion = 0.8, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_4to1S = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1, trainProportion = 0.8, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
SeqCVa_5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 5, trainProportion = 0.8, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 5, trainProportion = 0.8, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
# 1:1 ratios:
SeqCVa_1to1S = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1, trainProportion = 0.5, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_1to1S = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1, trainProportion = 0.5, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
SeqCVa_2F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 2, trainProportion = 0.5, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_2F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 2, trainProportion = 0.5, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
# 1:4 ratios:
SeqCVa_1to4S = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1, trainProportion = 0.2, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_1to4S = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1, trainProportion = 0.2, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
SeqCVa_i5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1/5, trainProportion = 0.2, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_i5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1/5, trainProportion = 0.2, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}



### LOAD FUNCTIONS
source(paste0(indir, "CV_Functions.R"))

# Run expand.grid0 manually, so we can tweak some of the settings
simSettingsDf = do.call(expand.grid0, simSettingsList)

# Compute mu from muFactor * 1/(2K)
simSettingsDf$mu = simSettingsDf$muFactor / (2 * simSettingsDf$K)
# Delete unused argument to avoid plyr errors
simSettingsDf$muFactor = NULL


### RUN AND SAVE SIMULATION RESULTS TO FILE
# (By default, file is .../SimData_YY-MM-DD.R)
# Pass in the expanded DF of settings, NOT the original list
SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed, mySuffix = mySuffix,
                        simSettingsDf = simSettingsDf)




#### NEXT, all the known-K methods ####

# Known K
FS_K3 = function(Y, X) {
  FS_KnownK(Y, X, K = 3)
}
FS_K5 = function(Y, X) {
  FS_KnownK(Y, X, K = 5)
}
FS_K10 = function(Y, X) {
  FS_KnownK(Y, X, K = 10)
}

# K=3
simSettingsList = list(nrSims = 1000L,
                       n = c(6250L),
                       p = c(11L, 50L, 250L),
                       K = c(3L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.t2",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("FS_FixedK")
FS_FixedK = FS_K3

simSettingsDf = do.call(expand.grid0, simSettingsList)
simSettingsDf$mu = simSettingsDf$muFactor / (2 * simSettingsDf$K)
simSettingsDf$muFactor = NULL

SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed,
                        mySuffix = paste0(mySuffix, "_K3"),
                        simSettingsDf = simSettingsDf)

# K=5
simSettingsList = list(nrSims = 1000L,
                       n = c(6250L),
                       p = c(11L, 50L, 250L),
                       K = c(5L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.t2",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("FS_FixedK")
FS_FixedK = FS_K5

simSettingsDf = do.call(expand.grid0, simSettingsList)
simSettingsDf$mu = simSettingsDf$muFactor / (2 * simSettingsDf$K)
simSettingsDf$muFactor = NULL

SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed,
                        mySuffix = paste0(mySuffix, "_K5"),
                        simSettingsDf = simSettingsDf)

# K=10
simSettingsList = list(nrSims = 1000L,
                       n = c(6250L),
                       p = c(11L, 50L, 250L),
                       K = c(10L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.t2",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("FS_FixedK")
FS_FixedK = FS_K10

simSettingsDf = do.call(expand.grid0, simSettingsList)
simSettingsDf$mu = simSettingsDf$muFactor / (2 * simSettingsDf$K)
simSettingsDf$muFactor = NULL

SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed,
                        mySuffix = paste0(mySuffix, "_K10"),
                        simSettingsDf = simSettingsDf)
