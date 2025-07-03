# The following functions help user simulate their own null distributions
# to estimate the p-value depending on the 2-way ANOVA design they have.

#' This function helps simulate the null for APCSSA test statistics.
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores in parallel or not
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRA and APCRCA) that will get standardized into APCSSA and a numeric vector with length equals to the numSim of all APCSSA statistics on the null data sets
#' @importFrom magrittr %>%
#' @export
sim_nullAPCSSA <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  .null1APCSSA(i, j, k, numSim, parallel)
  .null2APCSSA(i, j, k, numSim, parallel)
}

#' This function helps simulate the null for APCSSM test statistics.
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores in parallel or not
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRM and APCRCM) that will get standardized into APCSSM and a numeric vector with length equals to the numSim of all APCSSA statistics on the null data sets
#' @importFrom magrittr %>%
#' @export
sim_nullAPCSSM <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  .null1APCSSM(i, j, k, numSim, parallel)
  .null2APCSSM(i, j, k, numSim, parallel)
}
