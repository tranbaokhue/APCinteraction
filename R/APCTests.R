## First column contains the values of the observations;
## second contains the break-down of factor A (i of them);
## third contains the break-down of factor B (j of them);

# This file contains the core codes to run different tests for interaction
# in a balanced two-way ANOVA

######## APCSSA/APCSSM ########

## Test statistic for APCSSA ----

#' Title: APCSSA
#'
#' @description A more detailed explanation of what the function does.
#'
#' @param x Description of input parameter x
#' @param y Description of input parameter y
#' @return Description of the output
#' @examples
#' A(10, 20) # Example usage
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
