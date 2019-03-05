# Recreate Figure 6.5 from thesis


#### BAD SIGMA SIM SETTINGS ####

# Using "bad" Sigma, specifically the special case (K=10, p=11)
# where we saw a clear difference between good and bad.
# 
# Running
# (SeqCVa vs FullCVa) x (1S vs VF),
# at 3 train:test ratios (4:1, 1:1, and 1:4).
# Also add KnownK,
# all at several combos of n, mu,
# but setting k=10, p=11.

indir = "./Functions/"
outdir = "./MainSims/"
source(paste0(indir, "CV_Functions.R"))

# Specify the date/version whose sims to use:
mySuffix = "17-09-25_WorstCase"

source(paste0(outdir, "SimData_", mySuffix, ".R"))

# Rename with today's date/version for saving plots:
mySuffix = "18-01-09_WorstCase"

  
# May need to respecify indir and outdir after this,
# since it'll load what they were where the sims ran (e.g. hydra)
# and overwrite where you are plotting them (e.g. chromebook)
indir = "./Functions/"
outdir = "./MainSims/"

### PLOTTING
# Prep to plot
theme_set(theme_bw())
label_eq = function(labels) { label_both(labels, sep = " = ") }

# Calculate MOEs (2*SE)
SimData$MOE = with(SimData, calcMOEs(sd, nrSims))

# Plot just the main outcome, "foundTrueModel"
simToPlot = subset(SimData, outcome == "foundTrueModel")

# Drop columns we won't need for plotting
# because they only have one value anyway
colsToDrop = which(names(simToPlot) %in%
                     c("nrSims", "betamin", "betamax", "beta.fun",
                       "Sigma.fun", "epsilon.fun",
                       "doNormalize", "rescaleByN", "outcome"))
simToPlot = simToPlot[, -colsToDrop]

# For plotting, let's recreate muFactor
simToPlot$muFactor = simToPlot$mu * (2 * simToPlot$K - 1)

# What's left?
summary(simToPlot)


# Rename the KnownK procedure
summary(simToPlot$procedure)
levels(simToPlot$procedure) = c(levels(simToPlot$procedure)[1:12], "FS_KnownK")
summary(simToPlot$procedure)
summary(simToPlot)

# Reorder procedures, s.t. Seq>Full and F>S at each ratio
simToPlot$procedure = factor(simToPlot$procedure,
                             levels = levels(simToPlot$procedure)[c(13,
                                                                    3,4,1,2,
                                                                    7,8,5,6,
                                                                    11,12,9,10)])
summary(simToPlot$procedure)

# Plotting colors...
# We have 1+12 procedures, but it's
# 3 (ratios) * 2 (folds vs single-split) * 2 (Seq vs Full)
# So let's use Paired 6 to get 3 colors, each in dark & light,
# for the 3 ratios and 2 fold-vs-split options.
# Then use 2 linetypes, dash vs solid,
# for Seq vs Full.
# And black lines for KnownK
brewerColors = c("black", brewer.pal(12, "Paired")[c(2,2,1,1,
                                                     10,10,9,9,
                                                     6,6,5,5)])
myLty = c(1, rep(c(3,1), 6))

# For sake of plotting, redefine mu
# to say 1/(2K-1) or 1.9/(2K-1) directly in the labels
simToPlot$mu = factor(simToPlot$muFactor, levels = c(1, 1.9),
                      labels = c("1/(2K-1)", "1.9/(2K-1)"))

thisTitle = expression("X ~ N(0, "*Sigma[2]*"("*mu*")),  "*epsilon*" ~ N(0, 1),  "*beta*" ~ Unif(0.2, 2)")


## making a new version:
## all subplots for just a few lines (Seq, Folds)
##   to show big-picture effects.

brewerColors = c("black", brewer.pal(12, "Paired")[c(2,6)])
myLty = 1:3

simAllSubplots = subset(simToPlot,
                        procedure %in% c("FS_KnownK", "SeqCVa_5F",
                                         "SeqCVa_i5F"))
simAllSubplots$procedure = factor(simAllSubplots$procedure)
levels(simAllSubplots$procedure) =
  c("Known K", "SeqCV 5-fold", "SeqCV inv.-5-fold")

# Change colors and plotting characters for black-and-white figure
myPch = c(20, 83, 83)
myCol = grey(c(0, .3, .3))
p = ggplot(simAllSubplots, aes(x = n, y = mean,
                               ymin = mean - MOE, ymax = mean + MOE,
                               color = procedure,
                               linetype = procedure, shape = procedure))
p = p + geom_point(size = 3, fill = "black") + geom_line() + 
  geom_errorbar(width = .05) +
  scale_x_log10(breaks = c(50,250,1250,6250)) +
  xlab("n (log scale)") +
  ylab("Prob. of finding true model") +
  ggtitle(thisTitle) +
  facet_grid(K ~ mu + p, labeller = "label_eq", space = "free_x", scales = "free_x") + 
  scale_color_manual(values = myCol,
                     name = "Stopping rule") +
  scale_linetype_manual(values = myLty,
                        name = "Stopping rule") +
  scale_shape_manual(values = myPch,
                     name = "Stopping rule") +
  theme(legend.key.width = unit(2, "line"),
        panel.spacing.x = unit(.7, "line"))


ggsave(paste0(outdir, "SimPlot_", mySuffix, "_SeqFolds.pdf"), p,
       width = 6.5, height = 5)

