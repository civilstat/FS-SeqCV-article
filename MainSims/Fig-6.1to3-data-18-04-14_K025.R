# Create the simulated data used for Figures 6.1, 6.2, 6.3 in thesis:
# specifically, the known-K AND estimated-K simulations for K=25,
# at every n = 50, 250, 1250, 6250,
# at every p = 10, 50, 250, 1250.



# 4 sims at nCores=6 took 6750 sec,
# so 400 sims should take 675000 sec = 11250 min = 188 hrs = 7.8 days.
# (1000 sims' precision may be higher than we really need.)
# Starting on hydra3 around 1pm Sat April 14,
# expect to finish by 8am Sun April 22.

# Running JUST K=25 combos of the below sims

# Running
# SeqCVa VF vs FullCVa VF,
# at 2 train:test ratios (4:1 and 1:4).
# Also add KnownK,
# all at several combos of n, mu, k, p.
# Aiming for especially large k and p:
# k = 5, 25, 125,
# p = 10, 50, 250, 1250...
# REMEMBER to drop all combinations where k>p or k>n

mySeed = 180414.025
mySuffix = "18-04-14_K025"


# Keep doNormalize = FALSE.

# Use data settings:
# Sigma = ConstCorr
# K = 25
# p = 10, 50, 250, 1250
# n = 50, 250, 1250, 6250
# mu = 1/(2K), 5/(2K)
# beta_min = 1/5
# epsilon ~ N(0,1)
#
# and use 5F CV-a variants, with Seq CV and FullCV,
# at train/test ratios 4:1 and 1:4.
# Also KnownK with K=25



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
simSettingsList = list(nrSims = 400L,
                       n = c(50L, 250L, 1250L, 6250L),
                       p = c(10L, 50L, 250L, 1250L),
                       K = c(25L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.Normal",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("SeqCVa_5F", "SeqCVa_i5F",
                    "FullCVa_5F", "FullCVa_i5F",
                    "FS_K25")

# Specify procedure variants with non-default options.
# 4:1 ratios:
SeqCVa_5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 5, trainProportion = 0.8, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 5, trainProportion = 0.8, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
# 1:4 ratios:
SeqCVa_i5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1/5, trainProportion = 0.2, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = TRUE)
}
FullCVa_i5F = function(Y, X) {
  FS_CV_K(Y, X, nFolds = 1/5, trainProportion = 0.2, nSplits = 1,
          voteFolds = FALSE, voteSplits = TRUE, SeqCV = FALSE)
}
# Known K
FS_K25 = function(Y, X) {
  FS_KnownK(Y, X, K = 25)
}



### LOAD FUNCTIONS
source(paste0(indir, "CV_Functions.R"))

# Run expand.grid0 manually, so we can tweak some of the settings
simSettingsDf = do.call(expand.grid0, simSettingsList)

# Compute mu from muFactor * 1/(2K)
simSettingsDf$mu = simSettingsDf$muFactor / (2 * simSettingsDf$K)
# Delete unused argument to avoid plyr errors
simSettingsDf$muFactor = NULL

# Drop all rows where K>n (0% chance we can estimate true model)
# or where K>p (nonsensical to define)
simSettingsDf = subset(simSettingsDf,
                       K < n & K < p)

# Try to rearrange rows,
# so that you don't have all four huge sims in the same core:
# (max-p + max-n) with (min-p + min-n),
# (max-p + min-n) with vice versa,
# and all four moderate-p together.
rowshuffle = c(9:10, 3:8, 1:2, 11:12)
simSettingsDf = simSettingsDf[c(rowshuffle, 12+rowshuffle), ]

### RUN AND SAVE SIMULATION RESULTS TO FILE
# (By default, file is .../SimData_YY-MM-DD.R)
# Pass in the expanded DF of settings, NOT the original list
SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed, mySuffix = mySuffix,
                        simSettingsDf = simSettingsDf)
