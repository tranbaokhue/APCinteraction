
# APCinteraction <img src="https://img.shields.io/badge/R-≥4.0.0-blue" align="right" />

A nonparametric R package for detecting interaction in two-way ANOVA designs with balanced replication using **all possible crossed comparisons (APC)**. The test statistics in this package extend the methods of Hartlaub, Dean, and Wolfe (1999) to cases with multiple replications, using rank-based alignment techniques.

## Installation

You can install the development version of `APCinteraction` from [GitHub](https://github.com/) with:


``` r
# install.packages("pak")
pak::pak("tranbaokhue/APCinteraction")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(APCinteraction)
# Generate sample data with interaction and normal error
A <- factor(rep(1:3, each = 9))   # 3 levels
B <- factor(rep(rep(1:3, each = 3), times = 3))  # 3 levels, repeated
interaction <- as.numeric(as.character(A)) * as.numeric(as.character(B))
value <- interaction + rnorm(27)
data <- data.frame(value, A, B)

# Run the APCSSM test
APCSSA(value ~ A + B, data = data)
```


## Citation

Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024).  
*Nonparametric tests for interaction in two-way ANOVA with balanced replications*.  
arXiv preprint [arXiv:2410.04700](https://arxiv.org/abs/2410.04700).

