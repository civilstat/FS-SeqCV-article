# Create the simulated data used for Figure 6.4 in thesis:
# specifically, the simulations at the smaller n = 50, 250, 1250.



# 5 sims at nCores=6 took 199 sec,
# so 100 sims should take 4000 sec = 67 min,
# so let's do 1000 sims in 11.25 hrs.
# Starting on hydra2 around 4:35pm Sept 14,
# expect to finish by 4am Sept 15.

# Repeating 17-09-09 small-n sims, and still using ConstCorr Sigma,
# but now with "bad" t2 noise.
# 
# Running
# (SeqCVa vs FullCVa) x (1S vs VF),
# at 3 train:test ratios (4:1, 1:1, and 1:4).
# Also add KnownK,
# all at several combos of n, mu, k, p.


mySeed = 170914
mySuffix = "17-09-14_SmallN"


# Keep doNormalize = FALSE.

# Use data settings:
# Sigma = ConstCorr
# K = 3,5,10
# p = 11,50,250
# n = 50,250,1250
# mu = 1/(2K) * [1, 5]
# beta_min = 1/5
# epsilon ~ t2
#
# and use either 1S or 5F CV-a variants, with Seq vs Full CV,
# at train/test ratios 4:1, 1:1, and 1:4.
# Also KnownK with K=3,5,10
# (these are fast; we can drop the wrong-K runs later)



### SPECIFY FILE PATH, LOAD FUNCTIONS
indir = "./Functions/"
outdir = "./MainSims/"

# Flag whether or not to use parallel processing,
# and if so, how many cores to use
useParallel = TRUE; nCores = 6        # on department servers
# useParallel = TRUE; nCores = 2        # on desktop
# useParallel = FALSE; nCores = NULL    # for debugging w/o running in parallel




### SPECIFY SIMULATION SETTINGS
# Specify the simulation settings for which we want all cross-combinations,
# or any constants we want reported such as nrSims
# (but define whichProcedures separately, so each new dataset runs all procedures)
simSettingsList = list(nrSims = 1000L,
                       n = c(50L, 250L, 1250L),
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
                    "SeqCVa_1to4S", "FullCVa_1to4S", "SeqCVa_i5F", "FullCVa_i5F",
                    "FS_K3", "FS_K5", "FS_K10")

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

