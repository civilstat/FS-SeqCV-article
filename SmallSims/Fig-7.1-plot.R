# Recreate Figure 7.1 from thesis


#### NOTES ####

# Simulations for CV+1SE, vs KWW and HF inspired confidence sets:
# 1SE: any model whose point est is
#      less than the winner's 1SE away from the winner's point est.
# ZSE: similar but use estimated correlation with winner's neighbors
#      to form a single estimate of SE(difference between models),
#      then retain any model for which
#      Bonf-corrected CI for "diff from winner" overlaps with 0.
# KWW: any model whose Bonf-corrected CI
#      (for the mean) overlaps with the winner's.
# HF:  any model whose Bonf-corrected CI
#      *for the difference with winner* overlaps 0.


# Most basic thing we can try:
# using single-split CV and naive SE, and Bonferroni corrections,
# and a fixed model path so we can ignore post-sel'n inference...
# Vary SNR by changing n (keep same beta_min and sigma^2),
# also change k, p, mu,
# also try different split ratios.
# For each of our 3 conf.sets, record:
# - What is the actual coverage, say for nominal 95% sets?
#   Or if we figure out empirically approx. coverage for 1SE rule,
#   then *set* the KWW and HF rules to that *nominal* coverage,
#   what *actual* coverage do they get?
# - What is the average size of each confidence set?
#   Do we tend to get sizes ranked as 1SE < HF < KWW, as expected?
# - P(correct model) if we choose smallest model in the set?




#### FUNCTIONS ####

# We need functions to generate a new dataset,
# to train & test models along the (fixed, correct) model path,
# to construct each confidence set type from the matrix of test errors
#   (n_test rows by p columns, if model path is fixed),
# and finally to evaluate each confidence set's coverage, size, and
#   whether its smallest model was in fact the true model.
# We are going with the FullCV winner (global optimizer of test MSE),
#   NOT the SeqCV winner (first local optimizer).

library(matrixStats) # for colSds

make_data = function(n = 100, p = 10, k = 5, sigma = 1, beta = NULL) {
  stopifnot(k >= 1)  # haven't yet prepared for true models with k=0
  if(is.null(beta)) {
    beta = c(k:1, rep(0, p-k))
  } else {
    stopifnot(length(beta) == p)
  }
  X = as.data.frame(matrix(rnorm(n*p), nrow = n))
  names(X) = paste0("X", 1:p)
  E = rnorm(n, sd = sigma)
  Y = as.matrix(X) %*% beta + E
  return(data.frame(Y = Y, X))
}

train_and_test = function(dataset, train_prop = 0.8) {
  n = nrow(dataset)
  p = ncol(dataset) - 1
  n_train = ceiling(n * train_prop)
  n_test = n - n_train
  train = head(dataset, n_train)
  test = tail(dataset, n_test)
  
  # Fit a path of p+1 models, with up to p+1 parameters incl intercept
  betahats = matrix(0, nrow = p+1, ncol = p+1)
  colnames(betahats) = c("(Intercept)", names(dataset)[-1])
  betahats[1, ] = c(coef(lm(Y ~ 1, data = train["Y"])), rep(0, p))
  for(jj in 1:p){
    betahats[jj+1, ] = c(coef(lm(Y ~ .,
                                 data = train[, 1:(1+jj)])),
                         rep(0, p - jj))
  }
  
  Yhats = cbind(1, as.matrix(test[, -1])) %*% t(betahats)
  Resids = Yhats - test$Y
  colnames(Resids) = 0:p
  # Return a matrix with n_test rows by p+1 columns;
  # (j+1)'th column has all test residuals for model of size j
  return(Resids)  
}

SE_naive = function(mat) {
  stopifnot(length(dim(mat)) == 2)
  n = nrow(mat)
  return(colSds(mat) / sqrt(n))
}

ConfSet_1SE = function(Resids) {
  n_test = nrow(Resids)
  p = ncol(Resids) - 1
  MSEs = colMeans(Resids^2)
  col_winner = which.min(MSEs)
  SE_winner = SE_naive(Resids[, col_winner, drop = FALSE]^2)
  
  ColsConfSet = which(MSEs <= MSEs[col_winner] + SE_winner)
  # Should we also make this ConfSet contiguous:
  # trim any model with a gap past winner?
  # Hmm, glmnet doesn't, so let's not bother either.
  # Just subtract 1 from column values to get model sizes
  ModelSizeConfSet = ColsConfSet - 1
  return(unname(ModelSizeConfSet))
}

ConfSet_ZSEdiff = function(Resids, conf_level = 0.95) {
  n_test = nrow(Resids)
  p = ncol(Resids) - 1
  MSEs = colMeans(Resids^2)
  col_winner = which.min(MSEs)
  SE_winner = SE_naive(Resids[, col_winner, drop = FALSE]^2)
  
  # Avoid errors if min or max models are winner;
  # else, choose smaller of neighbors' cors, to get larger SEdiff
  if(col_winner == 1) {
    cor_neighbor = cor(Resids[, col_winner]^2, Resids[, col_winner+1]^2)
  } else if(col_winner == p+1) {
    cor_neighbor = cor(Resids[, col_winner]^2, Resids[, col_winner-1]^2)
  } else {
    cor_neighbor = min(cor(Resids[, col_winner]^2, Resids[, col_winner-1]^2),
                       cor(Resids[, col_winner]^2, Resids[, col_winner+1]^2))
  }
  
  # SE for difference, if all SEs and corrs were identical
  SE_diff = SE_winner * sqrt(2 * (1 - cor_neighbor))
  # Correct for p simultaneous CIs for diff. from winner
  BonfZ = -qnorm((1-conf_level)/(p*2))
  
  ColsConfSet = which(MSEs - MSEs[col_winner] <= BonfZ * SE_diff)
  ModelSizeConfSet = ColsConfSet - 1
  return(unname(ModelSizeConfSet))
}

ConfSet_KWW = function(Resids, conf_level = 0.95) {
  n_test = nrow(Resids)
  p = ncol(Resids) - 1
  # Correction for p+1 simultaneous CIs
  BonfZ = -qnorm((1-conf_level)/((p+1)*2))
  MSEs = colMeans(Resids^2)
  MSE_SEs = SE_naive(Resids^2)
  MSE_MOEs = MSE_SEs * BonfZ
  col_winner = which.min(MSEs)

  ColsConfSet = which(MSEs - MSE_MOEs <= (MSEs + MSE_MOEs)[col_winner])
  ModelSizeConfSet = ColsConfSet - 1
  return(unname(ModelSizeConfSet))
}

ConfSet_HF = function(Resids, conf_level = 0.95, Bonf = TRUE, two_sided = TRUE) {
  n_test = nrow(Resids)
  p = ncol(Resids) - 1
  sides = ifelse(two_sided, 2, 1)
  b = ifelse(Bonf, p, 1)
  # Correction for p simultaneous CIs for diff. from winner
  BonfZ = -qnorm((1-conf_level)/(b*sides))
  MSEs = colMeans(Resids^2)
  col_winner = which.min(MSEs)
  MSE_diffs = MSEs - MSEs[col_winner]
  Res_diffs_sq = Resids^2 - Resids[, col_winner]^2
  MSE_diff_SEs = SE_naive(Res_diffs_sq)
  MSE_diff_MOEs = MSE_diff_SEs * BonfZ

  ColsConfSet = which(MSE_diffs - MSE_diff_MOEs <= 0)
  ModelSizeConfSet = ColsConfSet - 1
  return(unname(ModelSizeConfSet))
}

eval_ConfSet = function(set, k) {
  if(length(set) > 0) {
    return(c(coverage = k %in% set,
             selection = k == min(set),
             size = length(set)))
  } else {
    return(c(coverage = FALSE, selection = FALSE, size = 0))
  }
}


#### SIMS ####

B = 1000

## Set up a structure for the output dataframe
out_df = data.frame(method = factor(c("1SE", "KWW", "HF", "HFnaive", "ZSEdiff")),
               coverage = 0.1, selection = 0.1, size = 0.1,
               k = 1L, n = 1L, p = 1L,
               train_prop = 0.1, conf_level = 0.1)

system.time({
  out_df_all = head(out_df, 0)
  for(cf_lv in c(.90, .95)) {
    for(kk in c(3, 10)) {
      for(nn in c(50, 250, 1250)) {
        for(pp in c(11, 50)) {
          for(tr_pr in c(.2, .8)) {
            out = replicate(B, {
              d = make_data(n = nn, p = pp, k = kk, sigma = 1,
                            beta = c(seq(2, .2, length.out = kk), rep(0, pp-kk)))
              r = train_and_test(d, train_prop = tr_pr)
              cs1 = ConfSet_1SE(r)
              e1 = eval_ConfSet(cs1, k = kk)
              cs2 = ConfSet_KWW(r, conf_level = cf_lv)
              e2 = eval_ConfSet(cs2, k = kk)
              cs3 = ConfSet_HF(r, conf_level = cf_lv)
              e3 = eval_ConfSet(cs3, k = kk)
              cs4 = ConfSet_HF(r, conf_level = cf_lv, Bonf = FALSE, two_sided = FALSE)
              e4 = eval_ConfSet(cs4, k = kk)
              cs5 = ConfSet_ZSEdiff(r, conf_level = cf_lv)
              e5 = eval_ConfSet(cs5, k = kk)
              return(c(e1, e2, e3, e4, e5))
            })
            out_matrix = matrix(colMeans(t(out)), 5, byrow = TRUE)
            colnames(out_matrix) = rownames(out)[1:3]
            out_df = cbind(data.frame(method = c("1SE", "KWW", "HF", "HFnaive", "ZSEdiff")),
                           out_matrix,
                           k = kk, n = nn, p = pp,
                           train_prop = tr_pr, conf_level = cf_lv)
            out_df_all = rbind(out_df_all, out_df)
          }
        }
      }
    }
  }
})
out_df_all$nr_sims = B

dump("out_df_all", file = "./SmallSims/SimResults_20180321_2confs.R")


#### PLOTS ####

library(ggplot2)
library(gridExtra)
source("./SmallSims/SimResults_20180321_2confs.R")
head(out_df_all)

out_df_all$covg_MOE = with(out_df_all,
                           2*sqrt(coverage*(1-coverage)/nr_sims))
out_df_all$seln_MOE = with(out_df_all,
                           2*sqrt(selection*(1-selection)/nr_sims))

out_df_all = subset(out_df_all, train_prop == 0.8)
out_df_all$method = factor(out_df_all$method,
                           levels = levels(out_df_all$method)[c(1,4,2,3,5)])

out_df_all$ConfLevel = out_df_all$conf_level

p1 = ggplot(out_df_all, aes(x = n, y = coverage,
                           ymin = coverage - covg_MOE,
                           ymax = coverage + covg_MOE,
                           color = method,
                           linetype = method))
p1 = p1 + geom_hline(aes(yintercept = ConfLevel)) +
  geom_point(shape = 1) + geom_line() + 
  geom_errorbar(width = .05) +
  scale_x_log10(breaks = unique(out_df_all$n)) +
  xlab("n") +
  ylab("Coverage") +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0.7, 1)) +
  facet_grid(p ~ ConfLevel + k, labeller = "label_both") +
  theme(legend.position = "top") +
  labs(color = "Method: ", linetype = "Method: ")




p2 = ggplot(out_df_all, aes(x = n, y = selection,
                           ymin = selection - seln_MOE,
                           ymax = selection + seln_MOE,
                           color = method,
                           linetype = method))
p2 = p2 + geom_point(shape = 1) + geom_line() + 
  geom_errorbar(width = .05) +
  scale_x_log10(breaks = unique(out_df_all$n)) +
  xlab("n") +
  ylab("Prob(correct selection)") +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(p ~ ConfLevel + k, labeller = "label_both") +
  theme(legend.position = "top") +
  labs(color = "Method: ", linetype = "Method: ")




p3 = ggplot(out_df_all, aes(x = n, y = size,
                           color = method,
                           linetype = method))
p3 = p3 + geom_point(shape = 1) + geom_line() + 
  scale_x_log10(breaks = unique(out_df_all$n)) +
  xlab("n") +
  ylab("Mean size of conf. set") +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(p ~ ConfLevel + k, labeller = "label_both", scales = "free_y") +
  theme(legend.position = "top", panel.spacing.x = unit(1, "lines")) +
  labs(color = "Method: ", linetype = "Method: ")




p1a = p1 + theme(legend.position = "none", panel.spacing.x = unit(1, "lines")) + xlab("")
p2a = p2 + theme(legend.position = "none", panel.spacing.x = unit(1, "lines")) + xlab("")
p3a = p3 + theme(legend.position = "bottom")
ggsave("./SmallSims/SimResults_All_20180321.pdf",
       grid.arrange(p1a, p2a, p3a, nrow = 3,
                    heights = c(5, 5, 6)),
       width = 6.5, height = 9)
