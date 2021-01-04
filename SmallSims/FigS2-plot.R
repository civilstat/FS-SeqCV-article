# Recreate Figure S2 from supplement


#### PLOT THE RESULTS ####
# From the run on Sept 4, 2017,
# using n=500, M=600, SNRs = c(1, 5, 10, 20, 40).

n = 500
ratios = 1:19/20
SNRs = c(1, 5, 10, 20, 40)
betamins = sqrt(SNRs / n)


pdf("./SmallSims/FigS2_PhatUnderOverFit_BothKs.pdf", width = 7.5, height = 4)
layout(matrix(1:2, 1))
par(mar = c(4, 3, 2, 1) + .1)

## k=5, p=6: ran for 4223 sec = 70 min
source("./SmallSims/FigS2_ProbsAtK5_20170904.R")
plot(ratios, Punders[, 1], type = 'l', xlim = 0:1, ylim = 0:1,
     ylab = "", xlab = expression(n[c]/n), las = 0)
title(ylab = "Probability", line = 2)
lines(ratios, Povers[, 1], lty = 2)
for(ii in 2:length(SNRs)) {
  lines(ratios, Punders[, ii], lwd = ii)
  lines(ratios, Povers[, ii], lty = 2)
}
legend("topleft", lty=1:2, bg = "white",
       legend = c(expression(paste(hat(P), "(underfit)")),
                  expression(paste(hat(P), "(overfit)"))))
title(main = expression(paste(k == 5)))
text(0.1, .05, "40", cex = 0.6, pos = 2)
text(0.2, .09, "20", cex = 0.6)
text(0.4, .16, "10", cex = 0.6)
text(0.6, .38, "5", cex = 0.6)
text(0.83, .88, expression(paste(nb[k-1]/sigma^2, " = 1")), cex = 0.6)




## k=10, p=11: ran for 8637 sec = 144 min
source("./SmallSims/FigS2_ProbsAtK10_20170904.R")
plot(ratios, Punders[, 1], type = 'l', xlim = 0:1, ylim = 0:1,
     ylab = "", xlab = expression(n[c]/n), las = 0)
lines(ratios, Povers[, 1], lty = 2)
for(ii in 2:length(SNRs)) {
  lines(ratios, Punders[, ii], lwd = ii)
  lines(ratios, Povers[, ii], lty = 2)
}
title(main = expression(paste(k == 10)))
text(0.11, .05, "40", cex = 0.6, pos = 2)
text(0.21, .09, "20", cex = 0.6)
text(0.4, .18, "10", cex = 0.6)
text(0.6, .42, "5", cex = 0.6)
text(0.93, .89, "1", cex = 0.6)

dev.off()
