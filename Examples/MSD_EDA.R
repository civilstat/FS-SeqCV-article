## Recreate Figure 6.7 from thesis
## as well as numerical results mentioned in Section 6.2.2,
## inline and in Table 6.1


#### Split up raw dataset ####

## These first few commands are RUN IN THE COMMAND LINE, not in R,
## in order to quickly split up full dataset (after unzipping)
## into the recommended learning:holdout split.

## Download from here (zipped as 201MB) and unzip (to 427MB):
## http://archive.ics.uci.edu/ml/datasets/YearPredictionMSD

## Then split the data so we will learn on the first 463,715 examples,
## and test (holdout) on the last 51,630 examples, as recommended.

## RUN IN THE COMMAND LINE:
# head -n 463715 YearPredictionMSD.txt > MSD_learn.txt
# tail -n 51630 YearPredictionMSD.txt > MSD_holdout.txt


#### Resave as .Rdata files ####

## Now go back to running code in R, not in command line.

## Load in each dataset, then save as .Rdata for faster loading next time.
## Also rename the data columns: Y, M1:12, V1:12, C1:66

# MSD_colnames = c("Y",
#                  paste0("Mean", 1:12),
#                  paste0("Var", 1:12),
#                  paste0("Cov", 1:66))

# system.time({
#   MSD_holdout = read.csv("./Examples/MSD_holdout.txt", header = FALSE)
#   names(MSD_holdout) = MSD_colnames
#   save(MSD_holdout, file = "./Examples/MSD_holdout.Rdata")
# }) ## Took 12 sec to read.csv() and save()
# rm("MSD_holdout")
# system.time(load("./Examples/MSD_holdout.Rdata")) ## Only takes 0.2 sec now

# system.time({
#   MSD_learn = read.csv("./Examples/MSD_learn.txt", header = FALSE)
#   names(MSD_learn) = MSD_colnames
#   save(MSD_learn, file = "./Examples/MSD_learn.Rdata")
# }) ## Took 130 sec to read.csv() and save()
# rm("MSD_learn")
# system.time(load("./Examples/MSD_learn.Rdata")) ## Only takes 2 sec now


#### Load in data and look at the correlations ####

load("./Examples/MSD_learn.Rdata")
load("./Examples/MSD_holdout.Rdata")

## learning data EDA:
## Check correlation matrix, among Xs as well as between Xs and Y.
## There *are* a few high corrs, but rare,
## so incoherence is NOT met but doesn't seem badly violated either...

## Condition number (of Xs alone):
eigs = eigen(cor(MSD_learn[, -1]), symmetric = FALSE, only.values = TRUE)
range(eigs$values)
sqrt(max(eigs$values) / min(eigs$values))
## 13.26, seems reasonable.
## For corr matrix, condition nr of 30 is "high" so 14 isn't so bad:
## http://faculty.cas.usf.edu/mbrannick/regression/Collinearity.html


## There is a tiny "long tail" of high corrs among Xs,
## but not much to worry about (and not much we CAN do, with OLS)
myBlue = blues9[7]
pdf("./Examples/MSD_Corrs_Original.pdf", width = 6, height = 1.5, pointsize = 10)
par(mar = c(5, 4, 0.5, 2.1))
layout(matrix(1:2, 1))
tmp = cor(MSD_learn[, -1])
diag(tmp) = NA
hist(tmp, seq(-1, 1, len = 81), col = myBlue, border = NA, xlim = c(-1, 1),
     xlab = "Correlations among predictors", main = "",
     freq = FALSE, ylim = c(0, 9), las = 1)
## And no X has a high corr with Y
hist(cor(MSD_learn[, 1], MSD_learn[, 2:91]), seq(-1, 1, len = 81),
     col = myBlue, border = NA, xlim = c(-1, 1),
     xlab = "Correlations between predictor and Year", main = "",
     freq = FALSE, ylim = c(0, 9), las = 1)
dev.off()



#### Check the null and full OLS models ####

## Function for evaluating trained models on holdout set.
testRMSE = function(mylm, testdata = MSD_holdout) {
  sqrt(mean((predict(mylm, newdata = testdata) - testdata$Y)^2))
}

## Fit both the null and full models
system.time({ lm0 = lm(Y ~ 1, data = MSD_learn) })  ## < 1 sec
system.time({ lm1 = lm(Y ~ ., data = MSD_learn) })  ## about 5 sec

## Evaluate them on holdout set
testRMSE(lm0)
## RMSE = 10.85
testRMSE(lm1)
## RMSE = 9.51
summary(lm1)$r.squared ## 0.24

## Difference in RMSEs between null and full models:
(testRMSE(lm0) - testRMSE(lm1))
## 1.34 is around 1 year and 4 months on scale of data,
## so even the full model is not substantially better than null.


#### Using Sec 5.2 rule of thumb for choosing a training ratio ####

## Say we only want to use a sparse model if there is substantial sparsity:
## we want no more than 1/3rd of the 90 predictors, so hope for k <= 30.

## First, what seems to be a reasonable betamin?
## (on the scale of standardized X data, i.e. mean 0 and SD 1 in each column of X)
MSD_learn_stdized = data.frame(Y = MSD_learn$Y, scale(MSD_learn[, -1]))
lm1_stdized = lm(Y ~ ., data = MSD_learn_stdized)
rev(sort(abs(coef(lm1_stdized)[-1])))
sum(abs(coef(lm1_stdized)[-1]) > 0.3)
## The 33 largest absolute betahats are each above 0.35,
## and two more are still above 0.3,
## so assuming |betamin| = 0.3 seems reasonable with k=30.

## Based on Section 5.2 heuristics,
## we can safely use a training ratio of n_c/n = 1/10
## if the following is at least sqrt(1 + 10/1) which is approx 3.32:
n = nrow(MSD_learn)
k = 30
betamin = 0.3
sigmahat = summary(lm1)$sigma
sqrt(n * betamin^2 / (k * sigmahat^2))  ## 3.9
## So indeed, using 1/10th of the training split for fitting the FS path
## and the other 9/10ths for deciding when to stop
## seems reasonable here.
