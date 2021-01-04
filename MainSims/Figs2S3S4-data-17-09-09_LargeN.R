# Create the simulated data used for Figures 2, S3, and S4 from paper and supplement:
# specifically, the estimated-K simulations at the largest n = 6250,
# for K=5, p=50 and 250, as well as some other settings not used in final paper.



# 5 sims at nCores=6 took 431 sec,
# so 100 sims should take 8620 sec = 144 min,
# so let's do 1000 sims in 24 hrs.
# Starting on hydra2 around 9:50am Sept 9,
# expect to finish by 10am Sept 10.

# Following on 17-09-08 test,
# do a full run now:
# (SeqCVa vs FullCVa) x (1S vs VF),
# at 3 train:test ratios (4:1, 1:1, and 1:4),
# plus KnownK,
# all at several combos of n, mu, k, p.
# 
# But splitting it up into several runs:
# one for all the settings with n=6250, and
# another for all the remaining settings, and
# a third for KnownK.
#
# THIS SCRIPT, LargeN, is for the "n=6250, CV" settings.

mySeed = 170909
mySuffix = "17-09-09_LargeN"


# Keep doNormalize = FALSE,
# Sigma = ConstCorr.

# Use data settings:
# K = 3,5,10
# p = 11,50,250
# n = 6250
# mu = 1/(2K) * [1, 5]
# beta_min = 1/5
# epsilon ~ Normal
#
# and use either 1S or 5F CV-a variants, with Seq vs Full CV,
# at train/test ratios 4:1, 1:1, and 1:4.



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
                       p = c(11L, 50L, 250L),
                       n = c(6250L),
                       K = c(3L, 5L, 10L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = c("epsilon.Normal"),
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


