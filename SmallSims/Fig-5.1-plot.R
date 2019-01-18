# Recreate Figure 5.1 from thesis

# Compute training ratio n_c/n for a given alpha,
# using only k=K_0-1,
# and SNR = b1/(sigma^2/n)
underfitTrainRatio = function(alpha, SNR, K0) {
  Z = qnorm(1-alpha/2)
  Chi = qchisq(alpha/2, df=1)
  num = (sqrt(SNR) - Z)^2 + Chi - Z^2
  trainRatio = 1/((num/K0) - 1)
  # If negative or over 1, nonsensical answer,
  # so n must be too small, so use LOO i.e. n_c = n-1
  trainRatio = ifelse(trainRatio <= 0 | trainRatio >= 1,
                      1, trainRatio)
  return( trainRatio )
}

alphas = seq(.0001, .1999, by = .0001)
SNRs = c(40, 50, 60, 100, 200, 400)




pdf("./SmallSims/ProbUnderOverFit_BothKs.pdf", width = 7.5, height = 4)
layout(matrix(1:2, 1))
par(mar = c(4, 3, 2, 1) + .1)
k = 5
plot(underfitTrainRatio(alphas, 400, k),
     alphas, xlim = c(-0.04, .99), ylim = c(0, 0.2),
     type = "l", lwd = 6,
     ylab = "", xlab = expression(n[c]/n))
title(ylab = "Probability", line = 2)
curve(1 - pchisq(1+1/x, df=1), from = 0.01, to = 1,
      add = TRUE, lty = 2)
for(ii in 1:5){
  mycurve = cbind(underfitTrainRatio(alphas, SNRs[ii], k), alphas)
  mycurve = mycurve[mycurve[, 1]<.99, ]
  lines(mycurve[, 1], mycurve[, 2], lwd = ii)
}
xmax = max(underfitTrainRatio(alphas, 400, k))
lines(c(xmax, 1), c(0, 0), lwd = 6)
legend("topright", lty=1:2, bg = "white",
       legend = c("P(underfit)", "P(overfit)"))
title(main = expression(paste(k == 5)))
text(-0.03, .003, "400", cex = 0.6)
text(0.085, .006, "200", cex = 0.6)
text(0.15, .014, "100", cex = 0.6)
text(0.3, .020, "60", cex = 0.6)
text(0.4, .028, "50", cex = 0.6)
text(0.65, .045, expression(paste(nb[k-1]/sigma^2, " = 40")), cex = 0.6)

k = 10
alphas = seq(.0001, .1999, by = .0001)
plot(underfitTrainRatio(alphas, 400, k),
     alphas, xlim = c(-0.04,.99), ylim = c(0, 0.2),
     type = "l", lwd = 6,
     ylab = "", xlab = expression(n[c]/n))
curve(1 - pchisq(1+1/x, df=1), from = 0.01, to = 1,
      add = TRUE, lty = 2)
for(ii in 1:5){
  mycurve = cbind(underfitTrainRatio(alphas, SNRs[ii], k), alphas)
  mycurve = mycurve[mycurve[, 1]<.99, ]
  lines(mycurve[, 1], mycurve[, 2], lwd = ii)
}
xmax = max(underfitTrainRatio(alphas, 400, k))
lines(c(xmax, 1), c(0, 0), lwd = 6)
title(main = expression(paste(k == 10)))
text(0.045, .007, "400", cex = 0.6, pos = 2)
text(0.13, .01, "200", cex = 0.6)
text(0.29, .014, "100", cex = 0.6)
text(0.59, .04, "60", cex = 0.6)
text(0.79, .06, "50", cex = 0.6)
text(0.95, .13, "40", cex = 0.6)
dev.off()
