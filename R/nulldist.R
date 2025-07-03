# The following functions help user simulate their own null distributions
# to estimate the p-value depending on the 2-way ANOVA design they have.

#' Simulate the null distribution for APCSSA
#'
#' @description
#' A complete function that simulate the null distribution for APCSSA and automatically save it to the working directory
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores in parallel or not
#'
#' @details
#' While we have provided null distributions for various settings in this package, there are still two-way settings that we don't readily have the null distribution to evaluate the significance of the test statistics (see References for the complete list of settings available). For example, one might wish to estimate the *p*-value with a higher `numSim`, ie. 250,000.
#'
#' The parameter `parallel` is defaulted to TRUE and will use 1 less than the total number of cores on user's device. This setting helps speed up the simulation by a factor of the number of cores utilized, so please refrain from changing the default unless the device do not work well with parallel-computing.
#'
#' @references
#' Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024). Nonparametric tests for interaction in two-way ANOVA with balanced replications. *arXiv preprint* arXiv:2410.04700.
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRA and APCRCA) that will get standardized into APCSSA and a numeric vector with length equals to the `numSim` of all APCSSA statistics on the null data sets
#'
#' @seealso
#' \code{\link{APCSSA}}, \code{\link{sim_nullAPCSSM}}
#'
#' @importFrom magrittr %>%
#' @export
sim_nullAPCSSA <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  .null1APCSSA(i, j, k, numSim, parallel)
  .null2APCSSA(i, j, k, numSim, parallel)
}

#' Simulate the null distribution for APCSSM
#'
#' @description
#' A complete function that simulate the null distribution for APCSSM and automatically save it to the working directory
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores in parallel or not
#'
#' @details
#' While we have provided null distributions for various settings in this package, there are still two-way settings that we don't readily have the null distribution to evaluate the significance of the test statistics (see References for the complete list of settings available). For example, one might wish to estimate the *p*-value with a higher `numSim`, ie. 250,000.
#'
#' The parameter `parallel` is defaulted to TRUE and will use 1 less than the total number of cores on user's device. This setting helps speed up the simulation by a factor of the number of cores utilized, so please refrain from changing the default unless the device do not work well with parallel-computing.
#'
#' @references
#' Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024). Nonparametric tests for interaction in two-way ANOVA with balanced replications. *arXiv preprint* arXiv:2410.04700.
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRM and APCRCM) that will get standardized into APCSSM and a numeric vector with length equals to the `numSim` of all APCSSM statistics on the null data sets
#'
#' @seealso
#' \code{\link{APCSSM}}, \code{\link{sim_nullAPCSSA}}
#'
#' @importFrom magrittr %>%
#' @export
sim_nullAPCSSM <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  .null1APCSSM(i, j, k, numSim, parallel)
  .null2APCSSM(i, j, k, numSim, parallel)
}
