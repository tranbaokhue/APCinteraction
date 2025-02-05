## First column contains the values of the observations;
## second contains the break-down of factor A (i of them);
## third contains the break-down of factor B (j of them);

# This file contains the core codes to run different tests for interaction
# in a balanced two-way ANOVA

######## APCSSA/APCSSM ########

## Test statistic for APCSSA ----

#' APCSSA
#'
#' @description This tests check All Possible Comparisons using rank-based calculations with alignment done using the average (mean).
#'
#' @param dataFrame Dataframe in long format with the observed values in the first column, Factor A in the second, and Factor B in column 3.
#' @param numTrial Since the critical values are generated using simulations, this number lets the functions know which critical values we are comparing the test statistics to.
#' @return The output is the APCSSA test statistics, which is the maximum of standardized statistics APCCRA and APCRCA.
#' @examples
#' library(Rfit)
#' data(BoxCox) # From Rfit Library
#' APCSSAts(BoxCox)
#'
#' @export
APCSSAts <- function(dataFrame, numTrial) {
  ## Get the scaled, unstandardized test statistics
  APCCRAD <- APCCRADts(dataFrame)
  APCRCAD <- APCRCADts(dataFrame)

  ## Standardize and pick the max of the two
  APCCRADstar <- (APCCRAD - APCSSnullDist[1])/APCSSnullDist[2]
  APCRCADstar <- (APCRCAD - APCSSnullDist[3])/APCSSnullDist[4]

  APCSSA <- max(APCCRADstar, APCRCADstar)
  return(APCSSA)
}

## Test statistic for APCSSM ----
#' @export
APCSSMts <- function(dataFrame, numTrial) {
  ## Get the scaled, unstandardized test statistics
  APCCRMD <- APCCRMDts(dataFrame)
  APCRCMD <- APCRCMDts(dataFrame)

  ## Standardize and pick the max of the two
  APCCRMDstar <- (APCCRMD - APCSSnullDist[5])/APCSSnullDist[6]
  APCRCMDstar <- (APCRCMD - APCSSnullDist[7])/APCSSnullDist[8]

  APCSSM <- max(APCCRMDstar, APCRCMDstar)
  return(APCSSM)
}
