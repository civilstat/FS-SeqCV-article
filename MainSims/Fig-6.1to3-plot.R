# Recreate Figures 6.1, 6.2, 6.3 from thesis

# Data-generating scripts are a bit messy,
# because they were run at different times,
# and I cobbled together the relevant results from several different runs
# for the final graphs in the thesis.
# In summary...

# 16-07-19 and 17-09-10_K5:
#   for smaller n (50, 250, 1250) or largest n (6250) respectively,
#   these are results for oracle FS with *known* K,
#   for K=5 and p=50 or 250,
#   from an earlier phase when we had other K and p settings not used in final thesis.
# 17-09-09_SmallN and 17-09-09_LargeN:
#   for smaller n (50, 250, 1250) or largest n (6250) respectively,
#   these are results for FS *estimating* K,
#   for K=5 and p=50 or 250,
#   from an earlier phase when we had other K and p settings not used in final thesis.
# 18-04-17_K005, 18-04-14_K025, and 18-04-14_K125:
#   for K = 5, 25, or 125 respectively,
#   these are results for FS with both *known* and *estimated* K,
#   at all n (50, 250, 1250, 6250) and at the newer p = 10, 1250 settings,
#   as well as the older p = 50, 250 settings for the newer K = 25, 125.

# Taken together, these 7 files generated the results for...
#   known K vs estimated K,
#   K = 5, 25, 125,
#   n = 50, 250, 1250, 6250,
#   p = 10, 50, 250.

# (For future papers based on this work,
# might be worthwhile to redo 18-04-17_K005 run at *all* settings of p,
# so we only need the equivalent of those three 18-04-1*_K*** files
# rather than the 7 files we have here.)

#################


# Revision of plots for the dissertation
# to show more of the story, incl. big-K and huge-P results
# as well as showing FP and FN rates.


#### GOOD SIM SETTINGS ####

### TABLES AND PLOTS

indir = "./Functions/"
outdir = "./MainSims/"
source(paste0(indir, "CV_Functions.R"))

# Load old sims:
# We only need these for K=5,
# since we are redoing sims for larger K separately.

# Specify the date/version whose sims to use:
mySuffix = "17-09-09_SmallN"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
SimData_SmallN = subset(SimData, K == 5)

# Redo with LargeN dataset
indir = "./Functions/"
outdir = "./MainSims/"
mySuffix = "17-09-09_LargeN"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
SimData_LargeN = subset(SimData, K == 5)


# Now add the various KnownK datasets:
# old sims at smaller n, plus new sims at largest n
# Get the KnownK cases for huge n=6250
# K=5
indir = "./Functions/"
outdir = "./MainSims/"
mySuffix = "17-09-10_K5"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
SimData_K5 = SimData

# KnownK at smaller n
indir = "./Functions/"
outdir = "./MainSims/"
mySuffix = "16-07-19"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
SimData_K = subset(SimData, K == 5)

# Drop all procedures whose K doesn't match true K:
SimData_K = subset(SimData_K, paste0("FS_K", K) == procedure)
# Combine KnownK datasets
SimData = rbind(SimData_K, SimData_K5)
# Redo the procedure name
SimData$procedure = factor("FS_KnownK")

# Combine with *all* datasets
SimData = rbind(SimData, SimData_SmallN, SimData_LargeN)

# Drop p=11
SimData = subset(SimData, p > 11)
# Drop the procedures we don't care about;
# this is all the "old" simulation data we need
SimData_old = subset(SimData,
                     procedure %in% c("FS_KnownK",
                                      "SeqCVa_5F", "FullCVa_5F",
                                      "SeqCVa_i5F", "FullCVa_i5F"))

# 18-04-17:
# Now also grab the later simulations, with new p=10,1250.
# (For K=5, the 04-17 results are same thing as 04-14 results
#  but with more sim repetitions.
#  Also add in the 04-14 results for K=25 and K=125.)
indir = "./Functions/"
outdir = "./MainSims/"
mySuffix = "18-04-17_K005"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
SimData_new = SimData
levels(SimData_new$procedure)[5] = "FS_KnownK"
summary(SimData_new$procedure)
# WAIT -- we need to remove part of the n=50,p=10 sims,
# because we intentionally doubled them up
# in order to ensure we can manually balance splitting sims over cores
# but now I don't have any unique ID I can use to drop the "extra" sims
# and they show up as double lines on the plot :(
# ...
# There are 450 rows: for each of the two muFactors 1 and 5,
# there are 225 rows, or 9 groups of 25,
# where each group of 25 is one of the 9 n*p settings,
# but the 1st two are repeats of each other.
# So I think we can just remove rows 1:25 and 226:250,
# and we'll be left with a single copy of those sim-settings.
SimData_new = SimData_new[-c(1:25, 226:250), ]

indir = "./Functions/"
outdir = "./MainSims/"
mySuffix = "18-04-14_K025"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
levels(SimData$procedure)[5] = "FS_KnownK"
summary(SimData$procedure)
SimData_new = rbind(SimData_new, SimData)

indir = "./Functions/"
outdir = "./MainSims/"
mySuffix = "18-04-14_K125"
source(paste0(outdir, "SimData_", mySuffix, ".R"))
levels(SimData$procedure)[5] = "FS_KnownK"
summary(SimData$procedure)
SimData_new = rbind(SimData_new, SimData)

# Combine with older sims
SimData = rbind(SimData_new, SimData_old)
SimData$procedure = factor(SimData$procedure)

# Fix mySuffix back to today, for saving plots
mySuffix = "18-04-17"
indir = "./Functions/"
outdir = "./MainSims/"

### PLOTTING
# Prep to plot
theme_set(theme_bw())
label_eq = function(labels) { label_both(labels, sep = " = ") }

# Calculate MOEs (2*SE)
SimData$MOE = with(SimData, calcMOEs(sd, nrSims))

simToPlot = SimData

# Drop columns we won't need for plotting
# because they (mostly?) only have one value anyway
colsToDrop = which(names(simToPlot) %in%
                     c("nrSims", "betamin", "betamax", "beta.fun",
                       "Sigma.fun", "epsilon.fun",
                       "doNormalize", "rescaleByN"))
simToPlot = simToPlot[, -colsToDrop]

# For plotting, let's recreate muFactor
simToPlot$muFactor = simToPlot$mu * (2 * simToPlot$K)


# What's left?
summary(simToPlot)


# Reorder procedures, s.t. Seq>Full and F>S at each ratio
simToPlot$procedure = factor(simToPlot$procedure,
  levels = levels(simToPlot$procedure)[c(5,1,3,2,4)])

levels(simToPlot$procedure) = 
  c("Known K",
    "Seq 5-fold CV", "Full 5-fold CV",
    "Seq inv.-5-fold CV", "Full inv.-5-fold CV")

# Plotting colors...
# Let's use Paired 6 to get 2 colors, each in dark & light,
# for the 2 ratios and 2 fold-vs-split options.
# Then use 2 linetypes, dashed vs dotted,
# for Seq vs Full.
brewerColors = c("black", brewer.pal(12, "Paired")[c(2,1,6,5)])
myLty = c(1, 2,2, 3,3)

# For sake of plotting, redefine mu
# to say 1/(2K) or 5/(2K) directly in the labels
simToPlot$mu = factor(simToPlot$muFactor, levels = c("1", "5"),
                      labels = c("1/(2K)", "5/(2K)"))

## 18-04-17:
## Make 3 plots, one for each outcome:
## foundTrueModel, nrFalsePos, nrFalseNegs
  
  
# foundTrueModel:
thisTitle = expression("X ~ N(0, "*Sigma[1]*"("*mu*")),  "*epsilon*" ~ N(0, 1),  "*beta*" ~ Unif(0.2, 2)")
# Change colors and plotting characters for black-and-white figure
myPch = c(20, 83, 70, 83, 70)
myCol = grey(c(0, .3, .7, .3, .7))
p = ggplot(subset(simToPlot, outcome == "foundTrueModel"),
           aes(x = n, y = mean,
               ymin = mean - MOE, ymax = mean + MOE,
               color = procedure, linetype = procedure, shape = procedure))
p = p + geom_point(size = 3, fill = "black") + geom_line() + 
  geom_errorbar(width = .05) +
  scale_x_log10(breaks = c(50,250,1250,6250)) +
  xlab("n") +
  ylab("Prob. of finding true model") +
  ggtitle(thisTitle) +
  scale_color_manual(values = myCol, name = "Stopping rule") +
  scale_linetype_manual(values = myLty, name = "Stopping rule") +
  scale_shape_manual(values = myPch, name = "Stopping rule") +
  facet_grid(K ~ mu + p, labeller = "label_eq") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

# Add spacing between just two of the facets,
# to better separate low- vs high-mu groups
# https://stackoverflow.com/questions/49123019/add-space-between-specific-facets-in-ggplot2-facet-grid
gt = ggplot_gtable(ggplot_build(p))
gt$widths[12] = 3*gt$widths[12]
grid.draw(gt)

ggsave(paste0(outdir, "SimPlot_", mySuffix, "_ProbSuccess.pdf"),
       grid.draw(gt),
       width = 9, height = 5)



# nrFalsePos:
p = ggplot(subset(simToPlot, outcome == "nrFalsePos"),
           aes(x = n, y = mean,
               ymin = mean - MOE, ymax = mean + MOE,
               color = procedure, linetype = procedure))
p = p + geom_point() + geom_line() + 
  geom_errorbar(width = .05) +
  scale_x_log10(breaks = c(50,250,1250,6250)) +
  xlab("n") +
  ylab("Avg. nr. of false positives") +
  ggtitle(thisTitle) +
  scale_color_manual(values = brewerColors, name = "Stopping rule") +
  scale_linetype_manual(values = myLty, name = "Stopping rule") +
  facet_grid(K ~ mu + p, labeller = "label_eq", scales = "free_y") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

gt = ggplot_gtable(ggplot_build(p))
gt$widths[12] = 2*gt$widths[12]
grid.draw(gt)

ggsave(paste0(outdir, "SimPlot_", mySuffix, "_FalsePos.pdf"),
       grid.draw(gt),
       width = 9, height = 5)



# nrFalseNegs:
p = ggplot(subset(simToPlot, outcome == "nrFalseNegs"),
           aes(x = n, y = mean,
               ymin = mean - MOE, ymax = mean + MOE,
               color = procedure, linetype = procedure))
p = p + geom_point() + geom_line() + 
  geom_errorbar(width = .05) +
  scale_x_log10(breaks = c(50,250,1250,6250)) +
  xlab("n") +
  ylab("Avg. nr. of false negatives") +
  ggtitle(thisTitle) +
  scale_color_manual(values = brewerColors, name = "Stopping rule") +
  scale_linetype_manual(values = myLty, name = "Stopping rule") +
  facet_grid(K ~ mu + p, labeller = "label_eq", scales = "free_y") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

gt = ggplot_gtable(ggplot_build(p))
gt$widths[12] = 2*gt$widths[12]
grid.draw(gt)

ggsave(paste0(outdir, "SimPlot_", mySuffix, "_FalseNegs.pdf"),
       grid.draw(gt),
       width = 9, height = 5)

