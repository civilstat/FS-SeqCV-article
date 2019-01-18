# Create the simulated data used for Figure 4.2 in thesis

# 7/12/17
# First draft of simulation for
# Prob(no overfit), for the "simple case"
# of orthogonal design and empty true model.
# Plot heatmaps (showing P(success)) on
# axes with p vs n_c.
# Show the contours follow p=sqrt(n_c)
# for Shao / wrapper CV.
# See what they look like for FS+CV.


# nr sims per setting
B = 500

p = 1*(1:8)
n = (20*(1:20))^2

nrP = length(p)
nrN = length(n)

nc_fun = function(x) floor(sqrt(x))


ProbSuccess = expand.grid(n = n, p = p,
                          ProbFS = NA, ProbWrapper = NA)

system.time({
  for(ii in 1:nrow(ProbSuccess)) {
    print(ii)
    nn = ProbSuccess$n[ii]
    pp = ProbSuccess$p[ii]
    ProbSuccess[ii, c("ProbFS", "ProbWrapper")] = rowMeans(replicate(B, {
      n_c = nc_fun(nn)
      n_v = nn - n_c
      X = replicate(pp, rnorm(nn))
      Y = rnorm(nn)
      X_c = X[1:n_c, , drop = FALSE]; X_v = X[(n_c+1):nn, , drop = FALSE]
      Y_c = Y[1:n_c]  ; Y_v = Y[(n_c+1):nn]
      
      SSE_0 = sum((Y_v - mean(Y_c))^2)
      
      h = which.max(abs(cor(Y_c, X_c)))
      df_c = data.frame(Y = Y_c, X = X_c[, h])
      df_v = data.frame(Y = Y_v, X = X_v[, h])
      lm_c = lm(Y ~ X, data = df_c)
      SSE_h = sum((Y_v - predict(lm_c, newdata = df_v))^2)
      
      SSE_wrapper = min(sapply(1:pp, function(hh) {
        df_c = data.frame(Y = Y_c, X = X_c[, hh])
        df_v = data.frame(Y = Y_v, X = X_v[, hh])
        lm_c = lm(Y ~ X, data = df_c)
        return(sum((Y_v - predict(lm_c, newdata = df_v))^2))
      }))
      return(c(SSE_0 <= SSE_h, SSE_0 <= SSE_wrapper))
    }))
  }
})
save(ProbSuccess, file = "./SmallSims/CVWrapper_20170716_SqrtN.Rdata")


