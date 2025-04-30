## First column contains the values of the observations;
## second contains the break-down of factor A (i of them);
## third contains the break-down of factor B (j of them);

# This file contains the codes for running different tests of interaction
# in a balanced two-way ANOVA

######## APCSSA/APCSSM ########
#' APCSSA
#'
#' @description
#' This function tests for interaction in a two-way layout using all possible crossed comparisons (APC) based on aligned ranks, with alignment performed using the mean.
#'
#' @param formula A formula specifying the model, with one response and two factors.
#' @param data A data frame in long format: the first column contains observed values, the second Factor A, and the third Factor B.
#' @param numSim An integer specifying the number of simulations used to estimate the null distribution. Defaults to 100000.
#'
#' @details
#' `APCSSA` performs a nonparametric test for interaction in a two-way layout with balanced replication. It computes two statistics—APCCRA and APCRCA—by aligning the data by means and ranking across rows or columns.
#'
#' These statistics use all possible crossed comparisons to detect interaction effects. The final test statistic is the maximum of the two standardized statistics. A p-value is estimated by comparing this statistic to a pre-simulated null distribution specific to the design dimensions.
#'
#' If the design is not covered by the package's pre-simulated settings or if higher precision is desired by increasing `numSim`, users must first generate the null distribution manually using \code{\link{sim_nullAPCSSA}},
#'
#' For more information, see the referenced article.
#'
#' @references
#' Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024). Nonparametric tests for interaction in two-way ANOVA with balanced replications. *arXiv preprint* arXiv:2410.04700.
#'
#' @return
#' A data frame with the following columns:
#' \itemize{
#'   \item \strong{Statistic}: The APCSSA test statistic (maximum of APCCRA and APCRCA).
#'   \item \strong{P-value}: The estimated p-value.
#' }
#' A summary table is printed to the console, and an interaction plot is generated.
#'
#' @seealso
#' \code{\link{nullAPCXXA}}, \code{\link{nullAPCSSA}}, \code{\link{sim_nullAPCSSA}}, \code{\link{APCSSM}}
#'
#' @examples
#' # Generate sample data with interaction and normal error
#' A <- factor(rep(1:3, each = 9))   # 3 levels
#' B <- factor(rep(rep(1:3, each = 3), times = 3))  # 3 levels, repeated
#' interaction <- as.numeric(as.character(A)) * as.numeric(as.character(B))
#' value <- interaction + rnorm(27)
#' data <- data.frame(value, A, B)
#'
#' # Run the APCSSA test
#' APCSSA(value ~ A + B, data = data)
#'
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
        nullAPCXXA(i, j, k)
        nullAPCSSA(i, j, k)
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
#' @description
#' This function tests for interaction in a two-way layout using all possible crossed comparisons (APC) based on aligned ranks, with alignment performed using the median.
#'
#' @param formula A formula specifying the model, with one response and two factors.
#' @param data A data frame in long format: the first column contains observed values, the second Factor A, and the third Factor B.
#' @param numSim An integer specifying the number of simulations used to estimate the null distribution. Defaults to 100000.
#'
#' @details
#' `APCSSM` performs a nonparametric test for interaction in a two-way layout with balanced replication. It computes two statistics—APCCRM and APCRCM—by aligning the data by means and ranking across rows or columns.
#'
#' These statistics use all possible crossed comparisons to detect interaction effects. The final test statistic is the maximum of the two standardized statistics. A p-value is estimated by comparing this statistic to a pre-simulated null distribution specific to the design dimensions.
#'
#' If the design is not covered by the package's pre-simulated settings or if higher precision is desired by increasing `numSim`, users must first generate the null distribution manually using \code{\link{sim_nullAPCSSM}}.
#'
#' For more information, see the referenced article.
#'
#' @references
#' Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024). Nonparametric tests for interaction in two-way ANOVA with balanced replications. *arXiv preprint* arXiv:2410.04700.
#'
#' @return
#' A data frame with the following columns:
#' \itemize{
#'   \item \strong{Statistic}: The APCSSM test statistic (maximum of APCCRM and APCRCM).
#'   \item \strong{P-value}: The estimated p-value.
#' }
#' A summary table is printed to the console, and an interaction plot is generated.
#'
#' @seealso
#' \code{\link{nullAPCXXM}}, \code{\link{nullAPCSSM}}, \code{\link{sim_nullAPCSSM}}, \code{\link{APCSSA}}
#'
#' @examples
#' # Generate sample data with interaction and normal error
#' A <- factor(rep(1:3, each = 9))   # 3 levels
#' B <- factor(rep(rep(1:3, each = 3), times = 3))  # 3 levels, repeated
#' interaction <- as.numeric(as.character(A)) * as.numeric(as.character(B))
#' value <- interaction + rnorm(27)
#' data <- data.frame(value, A, B)
#'
#' # Run the APCSSM test
#' APCSSM(value ~ A + B, data = data)
#'
#' @export
## Test statistic for APCSSM ----
APCSSM <- function(formula, data, numSim = 100000) {
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
  nullXXM_name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  nullSSM_name <- paste0("nullAPCSSM_", i, "x", j, "x", k)

  if (!exists(nullXXM_name, envir = .GlobalEnv) ||
      !exists(nullSSM_name, envir = .GlobalEnv)) {

    # If no stored null distributions, check if numSim == 100000
    if (numSim == 100000) {
      # Attempt to generate null distributions
      tryCatch({
        nullAPCXXM(i, j, k)
        nullAPCSSM(i, j, k)
      }, error = function(e) {
        stop("Error generating null distributions. Try running `sim_nullAPCSSM()` manually.")
      })
    } else {
      stop(paste("No valid null distribution found. Run `sim_nullAPCSSM()` to generate it."))
    }
  }

  # Retrieve the correct null distribution
  prev_name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  if (!exists(prev_name, envir = .GlobalEnv)) {
    stop("Error: The required null distribution does not exist. Run sim_nullAPCSSM() first.")
  }
  APCSSnullDist <- get(prev_name, envir = .GlobalEnv)

  ## Get the scaled, unstandardized test statistics
  APCCRMD <- .APCCRMD(dataFrame)
  APCRCMD <- .APCRCMD(dataFrame)

  ## Standardize and pick the max of the two
  APCCRMDstar <- (APCCRMD - APCSSnullDist[1])/APCSSnullDist[2]
  APCRCMDstar <- (APCRCMD - APCSSnullDist[3])/APCSSnullDist[4]

  APCSSM <- max(APCCRMDstar, APCRCMDstar)

  # Calculate p-value
  p_value <- .calc_p_value("APCSSM", i, j, k, APCSSM)

  # Display result
  result <- data.frame(
    Interaction = paste0(factorA_name, ":", factorB_name),
    Statistic = APCSSM,
    `P-value` = p_value
  )
  # Display result in aov-style format
  cat("\nAPCSSM Test Summary\n")
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

