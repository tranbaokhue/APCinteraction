## First column contains the values of the observations.
## Second column contains the break-down of factor A (i of them).
## Third contains the break-down of factor B (j of them).

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
#' @param numSim An integer specifying the number of simulations used to estimate the null distribution. Defaults to 100,000.
#'
#' @details
#' `APCSSA` performs a nonparametric test for interaction in a two-way layout with balanced replication. It computes two statistics - APCCRA and APCRCA - by aligning the data by means and ranking across rows or columns.
#'
#' These statistics use all possible crossed comparisons to detect interaction effects. The final test statistic is the maximum of the two standardized statistics. A *p*-value is estimated by comparing this statistic to a pre-simulated null distribution specific to the design dimensions.
#'
#' If the design is not covered by the package's existing collection of simulated settings or if a higher precision is desired, `numSim` > 100,000, users must first generate the null distribution manually using \code{\link{sim_nullAPCSSA}} before using this function to evaluate the interaction effect.
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
#'   \item \strong{*p*-value}: The estimated *p*-value.
#' }
#' A summary table is printed to the console, and an interaction plot is generated.
#'
#' @seealso
#' \code{\link{sim_nullAPCSSA}}, \code{\link{APCSSM}}
#'
#' @examples
#' # Set the seed for reproducibility
#' set.seed(206)
#'
#' # Parameters - number of levels for factors A and B, and replications
#' nA <- 3; nB <- 3; nrep <- 3
#'
#' # Generate levels and create the full design
#' A_vals <- 1:nA
#' B_vals <- 1:nB
#' design <- expand.grid(A = A_vals, B = B_vals)
#' design <- design[rep(seq_len(nrow(design)), each = nrep), ]
#'
#' # Compute response as interaction (product of numeric factor levels)
#' mu <- A_vals[design$A] * B_vals[design$B]
#' value <- mu + rnorm(length(mu))  # Add normal noise
#'
#' # Assemble the data
#' data <- data.frame(
#'   value = value,
#'   A = factor(design$A),
#'   B = factor(design$B)
#' )
#'
#' # Run the APCSSA test
#' APCSSA(value ~ A + B, data = data)
#'
#' @importFrom magrittr %>%
#' @export

## Test statistic for APCSSA ----
APCSSA <- function(formula, data, numSim = 100000) {
  # Extract the data frame from the formula
  df <- model.frame(formula, data)

  # Ensure the correct column structure (values, Factor A, Factor B)
  if (ncol(df) != 3) {
    stop("Formula must include exactly one response and two factors.")
  }

  # Pull out the response name
  response_name <- names(df)[1]

  # Turn the two RHS columns into factors
  fac <- setNames(lapply(df[2:3], factor), names(df)[2:3])

  # Order them by number of levels (smallest first)
  ord <- order(vapply(fac, nlevels, integer(1)))

  # Assign A/B names and values so that A always has ≤ levels than B
  factorA_name <- names(fac)[ord[1]]
  factorB_name <- names(fac)[ord[2]]
  dataFrame <- data.frame(
    value = df[[1]],
    A     = fac[[ord[1]]],
    B     = fac[[ord[2]]]
  )

  # Dimensions (guaranteed i ≤ j)
  i <- nlevels(dataFrame$A)
  j <- nlevels(dataFrame$B)
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
        .nullAPCXXA(i, j, k)
        .nullAPCSSA(i, j, k)
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
    `p-value` = p_value
  )

  # Display result in aov-style format
  cat(strrep("=", 40), "\n")
  cat("APCSSA Test Summary\n")
  cat(strrep("=", 40), "\n")

  print(result, row.names = FALSE)

  cat("\nNote: The p-value is estimated using the null distribution with", numSim, "simulations.\n\n")

  # Generate a color palette with `hcl.colors()`
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
#' @param numSim An integer specifying the number of simulations used to estimate the null distribution. Defaults to 100,000.
#'
#' @details
#' `APCSSM` performs a nonparametric test for interaction in a two-way layout with balanced replication. It computes two statistics - APCCRM and APCRCM - by aligning the data by means and ranking across rows or columns.
#'
#' These statistics use all possible crossed comparisons to detect interaction effects. The final test statistic is the maximum of the two standardized statistics. A *p*-value is estimated by comparing this statistic to a pre-simulated null distribution specific to the design dimensions.
#'
#' If the design is not covered by the package's existing collection of simulated settings or if a higher precision is desired, `numSim` > 100,000, users must first generate the null distribution manually using \code{\link{sim_nullAPCSSM}} before using this function to evaluate the interaction effect.
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
#'   \item \strong{*p*-value}: The estimated *p*-value.
#' }
#' A summary table is printed to the console, and an interaction plot is generated.
#'
#' @seealso
#' \code{\link{sim_nullAPCSSM}}, \code{\link{APCSSA}}
#'
#' @examples
#' # Set seed for reproducibility
#' set.seed(36)
#'
#' # Parameters - number of levels for factos A and B, replications, and interaction effect
#' nA <- 3; nB <- 4; nrep <- 2; c <- 1.25
#'
#'# Generate levels and create the full design
#' A_vals <- seq(-2, 2, length.out = nA)
#' B_vals <- seq(-1.5, 1.5, length.out = nB)
#' design <- expand.grid(A = 1:nA, B = 1:nB)
#' design <- design[rep(seq_len(nrow(design)), each = nrep), ]
#'
#' # Create the specific interaction matrix: alternate ±c in top rows
#' specInt <- matrix(0, nA, nB)
#' specInt[1:2, 1:nB] <- matrix(c(c, -c, -c, c), nrow = 2, byrow = TRUE)
#'
#' # Compute response
#' mu <- A_vals[design$A] + B_vals[design$B] + specInt[cbind(design$A, design$B)]
#' value <- mu + rt(length(mu), df = 1)
#'
#' # Assemble the data
#' data <- data.frame(
#'   value = value,
#'   A = factor(design$A, labels = round(A_vals, 2)),
#'   B = factor(design$B, labels = round(B_vals, 2))
#' )
#'
#' # With Cauchy errors, we opt for APCSSM to check for interaction
#' APCSSM(value ~ A + B, data = data)
#'
#' @importFrom magrittr %>%
#' @export
## Test statistic for APCSSM ----
APCSSM <- function(formula, data, numSim = 100000) {
  # Extract the data frame from the formula
  df <- model.frame(formula, data)

  # Ensure the correct column structure (values, Factor A, Factor B)
  if (ncol(df) != 3) {
    stop("Formula must include exactly one response and two factors.")
  }

  # Pull out the response name
  response_name <- names(df)[1]

  # Turn the two RHS columns into factors
  fac <- setNames(lapply(df[2:3], factor), names(df)[2:3])

  # Order them by number of levels (smallest first)
  ord <- order(vapply(fac, nlevels, integer(1)))

  # Assign A/B names and values so that A always has ≤ levels than B
  factorA_name <- names(fac)[ord[1]]
  factorB_name <- names(fac)[ord[2]]
  dataFrame <- data.frame(
    value = df[[1]],
    A     = fac[[ord[1]]],
    B     = fac[[ord[2]]]
  )

  # Dimensions (guaranteed i ≤ j)
  i <- nlevels(dataFrame$A)
  j <- nlevels(dataFrame$B)
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
        .nullAPCXXM(i, j, k)
        .nullAPCSSM(i, j, k)
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
    `p-value` = p_value
  )

  # Display result in aov-style format
  cat(strrep("=", 40), "\n")
  cat("APCSSM Test Summary\n")
  cat(strrep("=", 40), "\n")

  print(result, row.names = FALSE)

  cat("\nNote: The p-value is estimated using the null distribution with", numSim, "simulations.\n\n")

  # Generate a color palette with `hcl.colors()`
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

