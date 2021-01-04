# Create the simulated data used for Figures 2, S3, and S4 from paper and supplement:
# specifically, the KnownK simulations at the largest n = 6250,
# for K=5, p=50 and 250, as well as some other settings not used in final paper.



# 5 sims at nCores=6 took 8 sec,
# so 1000 sims should take 1600 sec = 27 min
# Starting around 2:30pm Sept 9,
# expect to finish 3pm Sept 9.

# We need to supplement the KnownK sims from 16-07-19 and 16-07-20,
# just for extra-large n, at every combo of K, p, and muFactor.
# THIS SCRIPT just does K=5 at muFactor=1,5.
# The other two scripts today do K=3 and then K=10,
# at each of muFactor=1,5.

mySeed = 1709105
mySuffix = "17-09-10_K5"

# Keep doNormalize = FALSE from now on.
# Also still keeping fixed for now:
# Sigma = ConstCorr, epsilon ~ Normal,
# voting across splits instead of avg'ing testRSS.

# Compare data settings:
# K = 5
# p = 11, 50, 250
# n = 6250
# mu = 1/(2K) * [1, 5]
# beta_min = 1/5
#
# and use just the FS alg variants
# FS_KnownK
# with K = 5



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
                       K = c(5L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.Normal",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("FS_K5")

# Specify procedure variants with non-default options
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


### RUN AND SAVE SIMULATION RESULTS TO FILE
# (By default, file is .../SimData_YY-MM-DD.R)
# Pass in the expanded DF of settings, NOT the original list
SimData = runAndSaveSim(indir, outdir,
                        simSettingsList = NULL, whichProcedures,
                        useParallel, nCores,
                        mySeed = mySeed, mySuffix = mySuffix,
                        simSettingsDf = simSettingsDf)

