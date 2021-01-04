# Create the simulated data used for Figures 2, S3, and S4 from paper and supplement:
# specifically, the known-K AND estimated-K simulations for K=5,
# at every n = 50, 250, 1250, 6250,
# at the p = 10 and 1250 not already simulated for K=5.



# 4 sims at nCores=6 took 5720 sec,
# so 400? sims should take 572000 sec = 9530 min = 159 hrs = 6.6 days.
# (1000 sims' precision may be higher than we really need.)
# Repeating 18-04-14, but doing all 400 sims at once.
# Hope to have 160hrs from 5:30pm Tues --> 9:30am Tues:
# Starting on hydra5 around 5:30pm Tues April 17,
# expect to finish by 9:30am Tues April 24.

# Running JUST K=5, p=(10,1250) combos of the below sims,
# since we have the other p settings for K=5 already.

# Running
# SeqCVa VF vs FullCVa VF,
# at 2 train:test ratios (4:1 and 1:4).
# Also add KnownK,
# all at several combos of n, mu, k, p.
# Aiming for especially large k and p:
# k = 5, 25, 125,
# p = 10, 50, 250, 1250...
# REMEMBER to drop all combinations where k>p or k>n

mySeed = 180417.005
mySuffix = "18-04-17_K005"


# Keep doNormalize = FALSE.

# Use data settings:
# Sigma = ConstCorr
# K = 5
# p = 10, 1250
# n = 50, 250, 1250, 6250
# mu = 1/(2K), 5/(2K)
# beta_min = 1/5
# epsilon ~ N(0,1)
#
# and use 5F CV-a variants, with Seq CV and FullCV,
# at train/test ratios 4:1 and 1:4.
# Also KnownK with K=5



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
                       p = c(10L, 1250L),
                       n = c(50L, 250L, 1250L, 6250L),
                       K = c(5L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.Normal",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("SeqCVa_5F", "SeqCVa_i5F",
                    "FullCVa_5F", "FullCVa_i5F",
                    "FS_K5")

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
FS_K5 = function(Y, X) {
  FS_KnownK(Y, X, K = 5)
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
# so that you don't have all huge sims in the same core;
# here there are 16 sims, so let's actually repeat 2 small ones --> 18
# so we can control exactly which 3 sims go in each of the 6 cores.
rowshuffle = c(1,1,8, 3,4,6, 2,5,7)
simSettingsDf = simSettingsDf[c(rowshuffle, 8+rowshuffle), ]


### RUN AND SAVE SIMULATION RESULTS TO FILE
# (By default, file is .../SimData_YY-MM-DD.R)
# Pass in the expanded DF of settings, NOT the original list
SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed, mySuffix = mySuffix,
                        simSettingsDf = simSettingsDf)

