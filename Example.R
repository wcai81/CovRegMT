## The following packages are required: Rcpp and RcppArmadillo

source("Functions.R")



## Illustrative example with simulated returns data

Rets <- read.table("Returns.txt")
Rets <- as.matrix(Rets)



## Procedures for k-FWER Control

# Required Inputs:
#   - data: A data matrix of dimensions T x N, where T is the number of observations and N is the number of variables.
#   - center: When set to TRUE, the data will be centered by subtracting means. When set to FALSE, the data is assumed to be already mean-zero before calling the function.
#   - k: An integer value such that k >= 1, defining the k-FWER to be controlled.
#   - step: Either set to "SS" for a single-step adjustment or "SD" for a step-down adjustment.
#   - B: Determines the number of resampling draws.
#   - alpha: A significance level in the interval (0, 1).

# Note: 
# - The value of B should be chosen so that alpha * B is an integer.

# Example:
#   - To control the 1-FWER with a single-step adjustment using B=100, and then setting to zero the correlations that are not significant at the 5% level:

kFWERresults <- reg_kFWER(data = Rets, center = TRUE, k = 1, step = "SS", B = 100, alpha = 0.05)

# Output:
# In the example above, 'kFWERresults' contains:
# - kFWERresults$covmat: the sample covariance matrix
# - kFWERresults$cormat: the sample correlation matrix
# - kFWERresults$pvalues: the matrix of kFWER-adjusted p-values 
# - kFWERresults$regcovmat: the regularized covariance matrix



## Procedures for FDP Control

# Required Inputs:
#   - data: A data matrix of dimensions T x N, where T is the number of observations and N is the number of variables.
#   - center: When set to TRUE, the data will be centered by subtracting means. When set to FALSE, the data is assumed to be already mean-zero before calling the function.
#   - gamma: A value in the interval [0, 1) representing the FDP exceedance threshold.
#   - step: Either set to "SS" for a single-step adjustment or "SD" for a step-down adjustment.
#   - B: Determines the number of resampling draws.
#   - alpha: A significance level in the interval (0, 1).

# Note: 
# - The value of B should be chosen so that alpha * B is an integer.

# Example:
#   - To control the probability that the FDP exceeds 10% using single-step adjusted p-values, computed with B=100, and then setting to zero the correlations that are not significant at the 5% level:

FDPresults <- reg_FDP(data = Rets, center = TRUE, gamma = 0.1, step = "SS", B = 100, alpha = 0.05)

# Output:
# In the example above, 'FDPresults' contains:
# - FDPresults$covmat: the sample covariance matrix
# - FDPresults$cormat: the sample correlation matrix
# - FDPresults$pvalues: the matrix of FDP-adjusted p-values 
# - FDPresults$regcovmat: the regularized covariance matrix

