## First column contains the values of the observations;
## second contains the break-down of factor A (i of them);
## third contains the break-down of factor B (j of them);

# This file contains the codes for running different tests of interaction
# in a balanced two-way ANOVA

######## APCSSA/APCSSM ########

# @examples
# library(Stat2Data)
# birdcalc <- BirdCalcium %>%select(Ca, Sex, Hormone)
# aov(Ca~factor(Hormone)*factor(Sex), data=birdcalc)
# APCSSA(Ca~factor(Hormone)*factor(Sex), data=birdcalc)

#' APCSSA
#'
#' @description This test uses all possible crossed comparisons (APC) based on aligned ranks with alignment done using the average (mean).
#'
#' @param dataFrame dataFrame should be in long format with the observed values in the first column, Factor A in the second, and Factor B in the third.
#' @param numSim Since the critical values are generated using simulations, this number lets the function know which critical values we are comparing the test statistics to.
#' @return The output is the APCSSA test statistics, which is the maximum of standardized statistics APCCRA and APCRCA.
#' @export
## Test statistic for APCSSA ----
APCSSA <- function(formula, data, numSim = 100000) {
  # Extract the data frame from the formula
  df <- model.frame(formula, data)

  # Ensure the correct column structure (values, Factor A, Factor B)
  if (ncol(df) != 3) {
    stop("Formula must include exactly one response and two factors.")
  }

  # Extract factor names from the formula
  response_name <- names(df)[1]
  factorA_name <- names(df)[2]
  factorB_name <- names(df)[3]

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

  # Retrieve the correct null distribution
  prev_name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  if (!exists(prev_name, envir = .GlobalEnv)) {
    stop("Error: The required null distribution does not exist. Run sim_nullAPCSSA() first.")
  }
  APCSSnullDist <- get(prev_name, envir = .GlobalEnv)

  ## Get the scaled, unstandardized test statistics
  APCCRAD <- .APCCRAD(dataFrame)
  APCRCAD <- .APCRCAD(dataFrame)

  ## Standardize and pick the max of the two
  APCCRADstar <- (APCCRAD - APCSSnullDist[1])/APCSSnullDist[2]
  APCRCADstar <- (APCRCAD - APCSSnullDist[3])/APCSSnullDist[4]

  APCSSA <- max(APCCRADstar, APCRCADstar)

  # Calculate p-value
  p_value <- .calc_p_value("APCSSA", i, j, k, APCSSA)

  # Display result
  result <- data.frame(
    Interaction = paste0(factorA_name, ":", factorB_name),
    Statistic = APCSSA,
    `P-value` = p_value
  )

  # Display result in aov-style format
  cat("\nAPCSSA Test Summary\n")
  print(result, row.names = FALSE)

  # Generate Viridis-like color palette with `hcl.colors()`
  num_levels <- nlevels(as.factor(dataFrame$B))
  color_palette <- hcl.colors(num_levels, "Viridis")

  # Generate interaction plot with gradient colors
  interaction.plot(
    x.factor = dataFrame$A,
    trace.factor = dataFrame$B,
    response = dataFrame$value,
    xlab = factorA_name,
    ylab = response_name,
    trace.label = factorB_name,
    col = color_palette,
    lwd = 2,
    main = "Interaction Plot"
  )
  # Return result invisibly for flexibility
  invisible(result)
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
