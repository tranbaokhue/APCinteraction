## First column contains the values of the observations;
## second contains the break-down of factor A (i of them);
## third contains the break-down of factor B (j of them);

# This file contains the core codes to run different tests for interaction
# in a balanced two-way ANOVA

######## APCSSA/APCSSM ########

# @examples
# library(Rfit)
# data(BoxCox) # From Rfit Library
# APCSSAts(BoxCox)

#' APCSSA
#'
#' @description This test uses all possible crossed comparisons (APC) based on aligned ranks with alignment done using the average (mean).
#'
#' @param dataFrame dataFrame should be in long format with the observed values in the first column, Factor A in the second, and Factor B in the third.
#' @param numTrial Since the critical values are generated using simulations, this number lets the function know which critical values we are comparing the test statistics to.
#' @return The output is the APCSSA test statistics, which is the maximum of standardized statistics APCCRA and APCRCA.
#'
#' @export
## Test statistic for APCSSA ----
APCSSA <- function(dataFrame, numTrial) {
  ## Get the scaled, unstandardized test statistics
  APCCRAD <- .APCCRAD(dataFrame)
  APCRCAD <- .APCRCAD(dataFrame)

  ## Standardize and pick the max of the two
  APCCRADstar <- (APCCRAD - APCSSnullDist[1])/APCSSnullDist[2]
  APCRCADstar <- (APCRCAD - APCSSnullDist[3])/APCSSnullDist[4]

  APCSSA <- max(APCCRADstar, APCRCADstar)
  return(APCSSA)
}


#' APCSSM
#'
#' @description This test uses all possible crossed comparisons (APC) based on aligned ranks with alignment done using the median.
#'
#' @param dataFrame dataFrame should be in long format with the observed values in the first column, Factor A in the second, and Factor B in the third.
#' @param numTrial Since the critical values are generated using simulations, this number lets the function know which critical values we are comparing the test statistics to.
#' @return The output is the APCSSM test statistics, which is the maximum of standardized statistics APCCRM and APCRCM.
#'
#' @export
## Test statistic for APCSSM ----
APCSSM <- function(dataFrame, numTrial) {
  ## Get the scaled, unstandardized test statistics
  APCCRMD <- .APCCRMD(dataFrame)
  APCRCMD <- .APCRCMD(dataFrame)

  ## Standardize and pick the max of the two
  APCCRMDstar <- (APCCRMD - APCSSnullDist[1])/APCSSnullDist[2]
  APCRCMDstar <- (APCRCMD - APCSSnullDist[3])/APCSSnullDist[4]

  APCSSM <- max(APCCRMDstar, APCRCMDstar)
  return(APCSSM)
}
