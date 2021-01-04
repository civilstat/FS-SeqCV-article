# Create the simulated data used for Figures 2, S3, and S4 from paper and supplement:
# specifically, the KnownK simulations at the smaller n = 50, 250, 1250,
# for K=5, p=50 and 250, as well as some other settings not used in final paper.

# Unfortunately, I didn't set a seed or datafile suffix when I ran this;
# but this script was how I generated SimData_16-07-19.R



# 4 sims at nCores=6 took 24.5/3 sec,
# so 1200 sims should take 2450 sec = 40 min
# Starting around 2:40pm July 19,
# expect to finish 3:30pm Jul 19.

# Rerun (almost) same settings as 16-06-22 sims,
# but with KNOWN true value of K.
# Basically, freq of finding true model
# will just tell us how often FS find the right model *path*
# separately from how often CV chooses the right K.
# This will let us calibrate how often
# each algorithm's mistakes are due to bad K, vs to bad model path.

# Because of how the multiple-sim code is structured,
# easiest if we just run all 3 procedures (fixed K values)
# on every data setting, then drop the "wrong K" procedures afterwards.
# So this wastes computation (3x as long as needed),
# but ought to run fast enough that this won't be major issue,
# since not doing any CV.

# ALSO, drop the high betamin setting,
# since we found it's equiv to using high n on low betamin.
# And drop muFactor = 1/5 setting, since equiv to muFactor = 1.


# Keep doNormalize = FALSE from now on.
# Also still keeping fixed for now:
# Sigma = ConstCorr, epsilon ~ Normal,
# voting across splits instead of avg'ing testRSS.

# Compare data settings:
# K = 3, 5, 10
# p = 11, 50, 250
# n = 50, 250, 1250
# mu = 1/(2K) * [1, 5]
# beta_min = 1/5
#
# and use just the FS alg variants
# FS_KnownK
# with each K = 3, 5, 10



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
simSettingsList = list(nrSims = 1200L,
                       p = c(11L, 50L, 250L),
                       n = c(50L, 250L, 1250L),
                       K = c(3L, 5L, 10L),
                       muFactor = c(1L, 5L),
                       betamin = c(1/5),
                       betamax = 2L,
                       beta.fun = "beta.Uniform",
                       Sigma.fun = "Sigma.ConstCorr",
                       epsilon.fun = "epsilon.Normal",
                       doNormalize = FALSE, rescaleByN = FALSE)
whichProcedures = c("FS_K3", "FS_K5", "FS_K10")

# Specify procedure variants with non-default options
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
                        # mySeed = mySeed, mySuffix = mySuffix,
                        simSettingsDf = simSettingsDf)

