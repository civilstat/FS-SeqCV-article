# Recreate Figure 4.2 in the thesis.

## Using n_c=sqrt(n) not n^(2/3),
## because this will let us argue that both
## p^2 = n_c and p^2 = n_v/n_c
## are contours for Wrapper, but not at all for FS+CV,
## so that *neither* restriction
## (p^2/n_c -> 0 nor p^2*n_c/n_v -> 0)
## is necessary for FS+CV to be consistent.

library(RColorBrewer)
library(plotrix)

###
# 07-16:
# Used n_c = sqrt(n)
load("./SmallSims/CVWrapper_20170716_SqrtN.Rdata")
p = 1:8
n = ceiling((20*(1:20))^(2))
nc_fun = function(x) floor(x^(1/2))
nrN = length(n)

# Try getting smoother contours:
# use 2d splines or such to fit a response surface,
# then do contour() with that smoothed-surface as z

# Clear out the attr garbage
ProbSuccess = as.data.frame(as.matrix(ProbSuccess))
# Add nc explicitly, for cleaner use of loess()
ProbSuccess$nc = nc_fun(ProbSuccess$n)

# Fit 2D loess smoother to each outcome
lo.wrap = loess(ProbWrapper ~ nc*p, data = ProbSuccess,
                span = .4)
contour(nc_fun(n), p, matrix(predict(lo.wrap, ProbSuccess), nrN),
        nlevels = 5)
lo.FS = loess(ProbFS ~ nc*p, data = ProbSuccess,
              span = .4)
contour(nc_fun(n), p, matrix(predict(lo.FS, ProbSuccess), nrN),
        levels = c(.9, .93, .96, .99))


# These look reasonably smooth enough to be pretty,
# but not crazy oversmoothed so you can't see anything.
# Just remember to report use of loess smoothing
# in the plot captions.

myColors = brewer.pal(9, "OrRd")[c(1:6, 8)]
# the 7th and 9th colors look funky when printed B&W









pdf("./SmallSims/CVBoth_20170716_SqrtN.pdf", width = 7, height = 8)
layout(matrix(1:2, nrow = 2))

## Top plot for Wrapper FS
par(mar = c(4, 3, 5, 7.5))
image(nc_fun(n), p, matrix(ProbSuccess$ProbWrapper, nrN),
      col = myColors,
      las = 1,
      xlab = "",
      ylab=  "",
      main = "Wrapper FS avoids overfitting\nwith fewer spurious variables & larger samples")
title(ylab = "Number of spurious variables, p", line = 2)
contour(nc_fun(n), p, matrix(predict(lo.wrap, ProbSuccess), nrN),
        nlevels = 5, add = TRUE,
        lwd = 2, labcex = .8)
color.legend(425, 1, 465, 8,
             legend = formatC(range(ProbSuccess$ProbWrapper), 2, format = "f"),
             rect.col = myColors,
             gradient = "y", align = "rb")
mtext("Prob(WrapperFS correctly picks null model)", side=4, line=5.7)

## Bottom plot for FS+CV
par(mar = c(5, 3, 4, 7.5))
image(nc_fun(n), p, matrix(ProbSuccess$ProbFS, nrN),
      col = myColors,
      las = 1,
      xlab = "",
      ylab=  "",
      main = "Forward Selection + Cross-Validation avoids overfitting\nwith more spurious variables & larger samples")
title(ylab = "Number of spurious variables, p", line = 2)
title(xlab = expression(paste("Training-sample size, ", n[c] == sqrt(n) %~~% frac(n[v], n[c]) )),
      line = 4)
contour(nc_fun(n), p, matrix(predict(lo.FS, ProbSuccess), nrN),
        levels = c(.9, .93, .96, .99), add = TRUE,
        lwd = 2, labcex = .8)
color.legend(425, 1, 465, 8,
             legend = formatC(range(ProbSuccess$ProbFS), 2, format = "f"),
             rect.col = myColors,
             gradient = "y", align = "rb")
mtext("Prob(FS+CV correctly picks null model)", side=4, line=5.7)

dev.off()
