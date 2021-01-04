# Create the simulated data used for Figure 3 in paper



# 5 sims at nCores=6 took 11 sec,
# so 100 sims should take 220 sec = 4 min,
# so let's do 1000 sims in 37 min.
# Starting on hydra4 around 2:15pm Sept 25,
# expect to finish by 3pm Sept 25.

# Using "bad" Sigma, specifically the special case (K=10, p=11)
# where we saw a clear difference between good and bad.
# 
# Running
# (SeqCVa vs FullCVa) x (1S vs VF),
# at 3 train:test ratios (4:1, 1:1, and 1:4).
# Also add KnownK,
# all at several combos of n, mu,
# but setting k=10, p=11.


mySeed = 170925
mySuffix = "17-09-25_WorstCase"


# Keep doNormalize = FALSE.

# Use data settings:
# Sigma = Sigma.NegKPosOthers with SigmaApprox = FALSE,
# K = 10
# p = 11
# n = 50,250,1250,6250
# mu = 1/(2K) * [1, 1.9]
# beta_min = 1/5
# epsilon ~ N(0,1)
#
# and use either 1S or 5F CV-a variants, with Seq vs Full CV,
# at train/test ratios 4:1, 1:1, and 1:4.
# Also KnownK with K=10



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
                       n = c(50L, 250L, 1250L, 6250L),
                       p = c(11L),
                       K = c(10L),
                       muFactor = c(1.0, 1.9),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.NegKPosOthers",
                       epsilon.fun = "epsilon.Normal",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("SeqCVa_4to1S", "FullCVa_4to1S", "SeqCVa_5F", "FullCVa_5F",
                    "SeqCVa_1to1S", "FullCVa_1to1S", "SeqCVa_2F", "FullCVa_2F",
                    "SeqCVa_1to4S", "FullCVa_1to4S", "SeqCVa_i5F", "FullCVa_i5F",
                    "FS_K10")

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
FS_K10 = function(Y, X) {
  FS_KnownK(Y, X, K = 10)
}



### LOAD FUNCTIONS
source(paste0(indir, "CV_Functions.R"))

# Run expand.grid0 manually, so we can tweak some of the settings
simSettingsDf = do.call(expand.grid0, simSettingsList)

# Compute mu from muFactor * 1/(2K-1)
simSettingsDf$mu = simSettingsDf$muFactor / (2 * simSettingsDf$K - 1)
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
