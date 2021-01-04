# Create the simulated data used for Figures Figures 2, S3, and S4 from paper and supplement:
# specifically, the known-K AND estimated-K simulations for K=125,
# at every n = 50, 250, 1250, 6250,
# at every p = 10, 50, 250, 1250.



# 4 sims at nCores=6 took 6540 sec,
# so 400 sims should take 654000 sec = 10900 min = 182 hrs = 7.6 days.
# (1000 sims' precision may be higher than we really need.)
# Starting on hydra2 around 1pm Sat April 14,
# expect to finish by 3am Sun April 22.

# Running JUST K=125 combos of the below sims

# Running
# SeqCVa VF vs FullCVa VF,
# at 2 train:test ratios (4:1 and 1:4).
# Also add KnownK,
# all at several combos of n, mu, k, p.
# Aiming for especially large k and p:
# k = 5, 25, 125,
# p = 10, 50, 250, 1250...
# REMEMBER to drop all combinations where k>p or k>n

mySeed = 180414.125
mySuffix = "18-04-14_K125"


# Keep doNormalize = FALSE.

# Use data settings:
# Sigma = ConstCorr
# K = 125
# p = 10, 50, 250, 1250
# n = 50, 250, 1250, 6250
# mu = 1/(2K), 5/(2K)
# beta_min = 1/5
# epsilon ~ N(0,1)
#
# and use 5F CV-a variants, with Seq CV and FullCV,
# at train/test ratios 4:1 and 1:4.
# Also KnownK with K=125



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
                       p = c(10L, 50L, 250L, 1250L),
                       n = c(50L, 250L, 1250L, 6250L),
                       K = c(125L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.Normal",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("SeqCVa_5F", "SeqCVa_i5F",
                    "FullCVa_5F", "FullCVa_i5F",
                    "FS_K125")

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
FS_K125 = function(Y, X) {
  FS_KnownK(Y, X, K = 125)
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
# so that you don't have both huge-n sims in the same core:
# (max-p + max-n) with (min-p + min-n),
# (max-p + min-n) with vice versa,
# and both moderate-n together.
simSettingsDf = simSettingsDf[c(1,6,3,4,5,2,
                                7,12,9,10,11,8), ]


### RUN AND SAVE SIMULATION RESULTS TO FILE
# (By default, file is .../SimData_YY-MM-DD.R)
# Pass in the expanded DF of settings, NOT the original list
SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed, mySuffix = mySuffix,
                        simSettingsDf = simSettingsDf)
