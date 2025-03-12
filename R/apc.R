## First column contains the values of the observations;
## second contains the break-down of factor A (i of them);
## third contains the break-down of factor B (j of them);

# This file contains the codes for running different tests of interaction
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
#' @param numSim Since the critical values are generated using simulations, this number lets the function know which critical values we are comparing the test statistics to.
#' @return The output is the APCSSA test statistics, which is the maximum of standardized statistics APCCRA and APCRCA.
#'
#' @export
## Test statistic for APCSSA ----
APCSSA <- function(formula, data, numSim) {
  # Extract the data frame from the formula
  df <- model.frame(formula, data)

  # Ensure the correct column structure (values, Factor A, Factor B)
  if (ncol(df) != 3) {
    stop("Formula must include exactly one response and two factors.")
  }

  # Prepare data in the required format
  dataFrame <- data.frame(
    value = df[, 1],
    A = df[, 2],
    B = df[, 3]
  )
  # Extract i, j, k from factor levels
  i <- nlevels(as.factor(dataFrame[, 2]))
  j <- nlevels(as.factor(dataFrame[, 3]))
  k <- nrow(dataFrame) / (i * j)

  # Check for available null distributions in the global environment
  nullXXA_name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  nullSSA_name <- paste0("nullAPCSSA_", i, "x", j, "x", k)

  if (!exists(nullXXA_name, envir = .GlobalEnv) ||
      !exists(nullSSA_name, envir = .GlobalEnv)) {

    # If no stored null distributions, check if numSim == 100000
    if (numSim == 100000) {
      # Attempt to generate null distributions
      tryCatch({
        nullAPCXXA(i, j, k, numSim)
        nullAPCSSA(i, j, k, numSim)
      }, error = function(e) {
        stop("Error generating null distributions. Try running `sim_nullAPCSSA()` manually.")
      })
    } else {
      stop(paste("No valid null distribution found. Run `sim_nullAPCSSA()` to generate it."))
    }
  }

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
#' @param numSim Since the critical values are generated using simulations, this number lets the function know which critical values we are comparing the test statistics to.
#' @return The output is the APCSSM test statistics, which is the maximum of standardized statistics APCCRM and APCRCM.
#'
#' @export
## Test statistic for APCSSM ----
APCSSM <- function(dataFrame, numSim) {
  ## Get the scaled, unstandardized test statistics
  APCCRMD <- .APCCRMD(dataFrame)
  APCRCMD <- .APCRCMD(dataFrame)

  ## Standardize and pick the max of the two
  APCCRMDstar <- (APCCRMD - APCSSnullDist[1])/APCSSnullDist[2]
  APCRCMDstar <- (APCRCMD - APCSSnullDist[3])/APCSSnullDist[4]

  APCSSM <- max(APCCRMDstar, APCRCMDstar)
  return(APCSSM)
}
