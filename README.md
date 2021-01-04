# FS-SeqCV-article
Simulations and data analyses from Jerzy Wieczorek and Jing Lei's article "Model-Selection Properties of Forward Selection and Sequential Cross-Validation for High-Dimensional Regression," to appear in *The Canadian Journal of Statistics* (accepted on Dec 22, 2020).

## Main simulations
The `MainSims` folder contains code (and saved simulation results) to replicate the main simulations from Section 5.1 (Figures 2 and 3) and Supplemental Section S3.3 (Figures S3, S4, and S5).

These simulations rely on functions defined (and libraries loaded) in `Functions/CV_Functions.R`.

## Examples using real data
The `Examples` folder contains code and data to replicate the real-data example from Section 5.2, including Figure 4 and Table 1.

The 201MB data file of the Million Song Dataset is not included, but `Fig4_MSD_EDA.R` contains instructions to download and pre-process it.

## Other smaller simulations
The `SmallSims` folder contains code (and saved simulation results) to replicate other simulations from throughout the article: Figure 1 in Section 3.3, and Figures S1 and S2 in Supplemental Section S3.2.

## `R` packages used for the article
The primary packages used are `leaps` for implementing Forward Selection, `doParallel` for parallel computing of the main simulations, and `ggplot2` for generating most of the figures.

A complete list of the `R` packages called in these scripts is: `doParallel`, `ggplot2`, `grid`, `gridExtra`, `Hmisc`, `leaps`, `MASS`, `Matrix`, `matrixStats`, `plotrix`, `plyr`, `R.utils`, and `RColorBrewer`.
