# Robocov

An R package for Robust Estimation of covariance and precision matrices when the data has missing entries.

Authors: [Kushal K Dey](http://kkdey.github.io/) and [Rahul Mazumder](http://www.mit.edu/~rahulmaz/)
Contact: kdey@hsph.harvard.edu

## How to cite

If you find Robocov useful, please cite:

Dey, K.K. and Mazumder, R. (2020). A convex optimization framework for gene-level tissue network estimation 
with missing data and its application in understanding disease architecture. bioRxiv.

Robocov software is a companion software to [CorShrink](https://github.com/kkdey/CorShrink) approach that uses adaptive shrinkage. If you use CorShrink, please cite

Dey, K.K. and Stephens, M. (2019). CorShrink : Empirical Bayes shrinkage estimation of correlations, with applications.
bioRxiv. Cold Spring Harbor Laboratory. 10.1101/368316. https://www.biorxiv.org/content/early/2018/07/24/368316.full.pdf

## Installation

`Robocov` requires the CVXR package that is available on CRAN. One can install Robocov from Github as follows.

```
install.packages("CVXR")
install_github('kkdey/Robocov')
```

Then load Robocov into R

```
library(Robocov)
```

## Licenses

The Robocov package is distributed under [GPL - General Public License (>= 2)]

## Demo

Lets start by loading an example samples by features data matrix (X)

```
data("sample_by_feature_data")
```

This data matrix has 544 samples, 53 variables and you will see a large number (~70%) of missing entries (NA) in this matrix. One can either use the standard sample correlation matrix.

```
standard = cor(sample_by_feature_data, use = "pairwise.complete.obs")
```
<embed src="vignettes/standard.pdf" width="800px" height="2100px" />

Alternatively, you can use the Robocov correlation estimator.

```
robocov = Robocov_cor(sample_by_feature_data, loss = "lasso")
```

<embed src="vignettes/robocov.pdf" width="800px" height="2100px" />

Analogously, you can estimate a sparse partial correlation matrix to better understand the causal structure.

<embed src="vignettes/probocov.pdf" width="800px" height="2100px" />

Observe that Robocov correlation and partial correlation estimators are visually more interpretable and less cluttered than the standard estimator. We also show in our manuscript that our method has better false positivity rate than standard approach and CorShrink.




