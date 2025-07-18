
# APCinteraction <img src="https://img.shields.io/badge/R-≥4.0.0-blue" align="right" />

A nonparametric R package for detecting interaction in two-way ANOVA designs with balanced replication using **all possible crossed comparisons (APC)**. The test statistics in this package extend the methods of Hartlaub, Dean, and Wolfe (1999) to cases with multiple replications, using rank-based alignment techniques.

## Installation

You can install the development version of `APCinteraction` from [GitHub](https://github.com/) with:


``` r
# Install devtools if you don't have it
install.packages("devtools")

# Install APCinteraction from GitHub
devtools::install_github("tranbaokhue/APCinteraction")
```

## Overview

The **APCinteraction** R package provides nonparametric tests for detecting interaction in two-way ANOVA designs with balanced replications. The core test statistics — **APCSSA** and **APCSSM** — are based on *All Possible Crossed Comparisons* (APC) and extend the methods of Hartlaub, Dean, and Wolfe (1999) to settings with replication, as recommended by Salazar-Alvarez et al. (2014).

These tests are designed to be robust and flexible, avoiding reliance on restrictive parametric assumptions. The methods utilize rank-based procedures with alignment to eliminate nuisance effects, offering enhanced power and interpretability in detecting interactions.

Key features include:

- Nonparametric testing of interaction in two-way designs
- Support for balanced designs with replication
- Methods based on mean or median alignment (APCSSA and APCSSM)
- Functions for simulation of null distributions when needed
- Integrated plotting and summary outputs

This package is particularly useful for researchers seeking robust alternatives to traditional ANOVA F-tests when normality assumptions may not hold.


## Examples

The following are basic examples which show you how to solve a common problem of determining whether there is interaction within a dataset.


### Example 1: Data with normal error (recommend: aov or APCSSA)

``` r
library(APCinteraction)
# Set the seed for reproducibility
set.seed(206)

# Parameters - number of levels for factors A and B, and replications
nA <- 3; nB <- 3; nrep <- 3

# Generate levels and create the full design
A_vals <- 1:nA
B_vals <- 1:nB
design <- expand.grid(A = A_vals, B = B_vals)
design <- design[rep(seq_len(nrow(design)), each = nrep), ]

# Compute response as interaction (product of numeric factor levels)
mu <- A_vals[design$A] * B_vals[design$B]
value <- mu + rnorm(length(mu))  # Add normal noise

# Assemble the data
data <- data.frame(
  value = value,
  A = factor(design$A),
  B = factor(design$B)
)

# Run the APCSSA test
APCSSA(value ~ A + B, data = data)
```


### Example 2: Data with Cauchy error (recommend: APCSSM)

``` r
# Set seed for reproducibility
set.seed(36)

# Parameters - number of levels for factos A and B, replications, and interaction effect
nA <- 3; nB <- 4; nrep <- 2; c <- 1.25

# Generate levels and create the full design
A_vals <- seq(-2, 2, length.out = nA)
B_vals <- seq(-1.5, 1.5, length.out = nB)
design <- expand.grid(A = 1:nA, B = 1:nB)
design <- design[rep(seq_len(nrow(design)), each = nrep), ]

# Create the specific interaction matrix: alternate ±c in top rows
specInt <- matrix(0, nA, nB)
specInt[1:2, 1:nB] <- matrix(c(c, -c, -c, c), nrow = 2, byrow = TRUE)

# Compute response
mu <- A_vals[design$A] + B_vals[design$B] + specInt[cbind(design$A, design$B)]
value <- mu + rt(length(mu), df = 1)

# Assemble the data
data <- data.frame(
  value = value,
  A = factor(design$A, labels = round(A_vals, 2)),
  B = factor(design$B, labels = round(B_vals, 2))
)

# With Cauchy errors, we opt for APCSSM to check for interaction
APCSSM(value ~ A + B, data = data)
```


## Citation

Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024).  *Nonparametric tests for interaction in two-way ANOVA with balanced replications*.  arXiv preprint [arXiv:2410.04700](https://arxiv.org/abs/2410.04700).

