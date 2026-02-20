
# APCinteraction <img src="https://img.shields.io/badge/R-≥4.0.0-blue" align="right" />

A nonparametric R package for detecting interaction in two-way ANOVA designs with balanced replication using **all possible crossed comparisons (APC)**. The test statistics in this package extend the methods of Hartlaub, Dean, and Wolfe (1999) to cases with multiple replications, using rank-based alignment techniques.

## Installation

You can install the development version of `APCinteraction` from [GitHub](https://github.com/) with:

``` r
# Install devtools if you don't have it
install.packages("devtools")

# Install APCinteraction from GitHub
devtools::install_github("tranbaokhue/APCinteraction")

# Load the package
library(APCinteraction)
```

## Companion data package

The companion package **APCinteractionData** provides 90 pre-computed null distributions for each test statistic (100,000 simulations each), covering two-way designs with 2 to 6 factor levels and 1 to 10 replications. Install it for fast _p_-value estimation:

``` r
devtools::install_github("tranbaokhue/APCinteractionData")
```

This companion package is optional. Two null distribution files are bundled in **APCinteraction** for the examples, and users can generate their own with `sim_nullAPCSSA()` / `sim_nullAPCSSM()`.

## Overview

The **APCinteraction** R package provides nonparametric tests for detecting interaction in two-way ANOVA designs with balanced replications. The core test statistics — **APCSSA** and **APCSSM** — are based on *All Possible Comparisons* (APC) and extend the methods of Hartlaub, Dean, and Wolfe (1999) to settings with replication, as recommended by Salazar-Alvarez et al. (2014).

These tests are designed to be robust and flexible, avoiding reliance on restrictive parametric assumptions. The methods utilize rank-based procedures with alignment to eliminate nuisance effects, offering enhanced power in detecting interactions.

Key functions the package provide are:

### Functions for calculating the test statistics

These statistics use All Possible Comparisons (APC) to detect interaction effects. The final test statistic is the maximum of the two standardized statistics. A _p_-value is estimated by comparing this statistic to a pre-simulated null distribution specific to the design dimensions. In addition to providing the statistics and estimated _p_-value in the summary table, the function will also output an interaction plot to help visualize potential interaction.

``` r
# Test using mean alignment
APCSSA(formula, data, numSim = 1e+05)

# Test using median alignment
APCSSM(formula, data, numSim = 1e+05)
```

If the design setting you are hoping to check for interaction is not among the pre-simulated ones in **APCinteractionData**, we also provide functions to simulate the null distribution on your own device.

### Functions related to simulating null distributions

While the companion package **APCinteractionData** provides an extensive collection of 90 pre-simulated null distributions (with 100,000 replications), there are definitely still scenarios where that is not sufficient. For example, a researcher might have an experiment with more factor levels or want to guarantee the accuracy of the estimated _p_-value by comparing it with null distributions generated from 150,000 and 250,000 replications.

Thus, we also provide functions for generating your own null distributions given a design and replication number. These functions are designed with the same structure as the function to calculate the test statistics for ease of use.

``` r
# Simulation for APCSSA null distributions
sim_nullAPCSSA(formula, data, numSim = 1e+05)

# Simulation for APCSSM null distributions
sim_nullAPCSSM(formula, data, numSim = 1e+05)
```

Accompanying the two simulation functions above, we provide users a function to quickly save the simulated values to a specified directory so that they are accessible by **APCSSA/APCSSM** even in other R sessions. This is optional, but is highly recommended!

``` r
save_null(type, i, j, k, path)
```

## Examples

The following are basic examples which show you how to solve a common problem of determining whether there is interaction within a dataset.


### Example 1: Data with normal error (recommend: aov or APCSSA)

``` r
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

# Run the APCSSA test and view summary
result <- APCSSA(value ~ A + B, data = data)
summary.APC(result)
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
result <- APCSSM(value ~ A + B, data = data)
summary.APC(result)
```


## Citation

Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024).  *Nonparametric tests for interaction in two-way ANOVA with balanced replications*.  arXiv preprint [arXiv:2410.04700](https://arxiv.org/abs/2410.04700).

