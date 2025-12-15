# Package-internal environment for storing null distributions ----
.apc_cache <- new.env(parent = emptyenv())

# Functions for APCSSA/APCSSM statistics ----

## APCCRA ----
.APCCRAD <- function(dataFrame) {
  ## Number of levels of factors A (I) and B (J), and number of replications per cell (K)
  I <- nlevels(as.factor(dataFrame[, 2]))
  J <- nlevels(as.factor(dataFrame[, 3]))
  K <- nrow(dataFrame)/(I*J)

  ## Split dataframe into smaller dataframes by factor B (by columns) for aligning
  dataList <- split(dataFrame, as.factor(dataFrame[, 3]))

  ## Align by column average
  for (j in 1:J) { ## looping over the number of levels for factor B (J)
    dataList[[j]][, 1] <- dataList[[j]][, 1] - mean(dataList[[j]][, 1])
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor A (by rows) for ranking
  dataList <- split(dataFrame, as.factor(dataFrame[, 2]))

  ## Rank by row
  for (i in 1:I) { ## looping over the number of levels for factor A (I)
    dataList[[i]][, 1] <- rank(dataList[[i]][, 1], ties.method = "average")
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor B (by column) and then by factor A (by row)
  ## for calculation
  dataList <- split(dataFrame, as.factor(dataFrame[, 3]))
  for (j in 1:J) {
    dataList[[j]] <- split(dataList[[j]], as.factor(dataList[[j]][, 2]))
  }

  ## Calculate all the J*(J-1)/2 values of V_jj'
  APCCRA <- c()
  for (j1 in (1:(J - 1))) {
    for (j2 in ((j1 + 1):J)) {
      Vjjp <- 0
      ## Calculate each cross-comparison
      for (i1 in (1:(I - 1))) {
        for (i2 in ((i1 + 1):I)) {
          for (k1 in 1:K) {
            for (k2 in 1:K) {
              for (k3 in 1:K) {
                for (k4 in 1:K) {
                  Vjjp <- Vjjp + (dataList[[j1]][[i1]][, 1][k1] + dataList[[j2]][[i2]][, 1][k2] - dataList[[j1]][[i2]][, 1][k3] - dataList[[j2]][[i1]][, 1][k4])^2
                }
              }
            }
          }
        }
      }
      APCCRA <- c(APCCRA, Vjjp)
    }
  }
  APCCRA <- max(APCCRA)

  ## Scale
  APCCRAD <- 2*APCCRA/(K^4*I*(I - 1))
  return(APCCRAD)
}

## APCRCA ----
.APCRCAD <- function(dataFrame) {
  ## Number of levels of factors A (I) and B (J), and number of replication per cell (K)
  I <- nlevels(as.factor(dataFrame[, 2]))
  J <- nlevels(as.factor(dataFrame[, 3]))
  K <- nrow(dataFrame)/(I*J)

  ## Split dataframe into smaller dataframes by factor A (by rows) for aligning
  dataList <- split(dataFrame, as.factor(dataFrame[, 2]))

  ## Aligning by row average
  for (i in 1:I) { ## looping over the number of levels for factor A (I)
    dataList[[i]][, 1] <- dataList[[i]][, 1] - mean(dataList[[i]][, 1])
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor B (by columns) for ranking
  dataList <- split(dataFrame, as.factor(dataFrame[, 3]))

  ## Rank by column
  for (j in 1:J) { ## looping over the number of levels for factor B (J)
    dataList[[j]][, 1] <- rank(dataList[[j]][, 1], ties.method = "average")
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor A (by row) and then by factor B (by column)
  ## for calculation
  dataList <- split(dataFrame, as.factor(dataFrame[, 2]))
  for (i in 1:I) {
    dataList[[i]] <- split(dataList[[i]], as.factor(dataList[[i]][, 3]))
  }

  ## Calculate all the I*(I-1)/2 values of V_ii'
  APCRCA <- c()
  for (i1 in (1:(I - 1))) {
    for (i2 in ((i1 + 1):I)) {
      Viip <- 0
      ## Calculate each cross-comparison
      for (j1 in (1:(J - 1))) {
        for (j2 in ((j1 + 1):J)) {
          for (k1 in 1:K) {
            for (k2 in 1:K) {
              for (k3 in 1:K) {
                for (k4 in 1:K) {
                  Viip <- Viip + (dataList[[i1]][[j1]][, 1][k1] + dataList[[i2]][[j2]][, 1][k2] - dataList[[i1]][[j2]][, 1][k3] - dataList[[i2]][[j1]][, 1][k4])^2
                }
              }
            }
          }
        }
      }
      APCRCA <- c(APCRCA, Viip)
    }
  }
  APCRCA <- max(APCRCA)

  ## Scale
  APCRCAD <- 2*APCRCA/(K^4*J*(J - 1))
  return(APCRCAD)
}

## APCCRM ----
.APCCRMD <- function(dataFrame) {
  ## Number of levels of factors A (I) and B (J), and number of replication per cell (K)
  I <- nlevels(as.factor(dataFrame[, 2]))
  J <- nlevels(as.factor(dataFrame[, 3]))
  K <- nrow(dataFrame)/(I*J)

  ## Split dataframe into smaller dataframes by factor B (by columns) for aligning
  dataList <- split(dataFrame, as.factor(dataFrame[, 3]))

  ## Align by column median
  for (j in 1:J) { ## looping over the number of levels for factor B (J)
    dataList[[j]][, 1] <- dataList[[j]][, 1] - median(dataList[[j]][, 1])
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor A (by rows) for ranking
  dataList <- split(dataFrame, as.factor(dataFrame[, 2]))

  ## Rank by row
  for (i in 1:I) { ## looping over the number of levels for factor A (I)
    dataList[[i]][, 1] <- rank(dataList[[i]][, 1], ties.method = "average")
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor B (by column) and then by factor A (by row)
  ## for calculation
  dataList <- split(dataFrame, as.factor(dataFrame[, 3]))
  for (j in 1:J) {
    dataList[[j]] <- split(dataList[[j]], as.factor(dataList[[j]][, 2]))
  }

  ## Calculate all the J*(J-1)/2 values of V_jj'
  APCCRM <- c()
  for (j1 in (1:(J - 1))) {
    for (j2 in ((j1 + 1):J)) {
      Vjjp <- 0
      ## Calculate each cross-comparison
      for (i1 in (1:(I - 1))) {
        for (i2 in ((i1 + 1):I)) {
          for (k1 in 1:K) {
            for (k2 in 1:K) {
              for (k3 in 1:K) {
                for (k4 in 1:K) {
                  Vjjp <- Vjjp + (dataList[[j1]][[i1]][, 1][k1] + dataList[[j2]][[i2]][, 1][k2] - dataList[[j1]][[i2]][, 1][k3] - dataList[[j2]][[i1]][, 1][k4])^2
                }
              }
            }
          }
        }
      }
      APCCRM <- c(APCCRM, Vjjp)
    }
  }
  APCCRM <- max(APCCRM)
  ## Scale
  APCCRMD <- 2*APCCRM/(K^4*I*(I - 1))
  return(APCCRMD)
}


## APCRCM ----
.APCRCMD <- function(dataFrame) {
  ## Number of levels of factors A (I) and B (J), and number of replication per cell (K)
  I <- nlevels(as.factor(dataFrame[, 2]))
  J <- nlevels(as.factor(dataFrame[, 3]))
  K <- nrow(dataFrame)/(I*J)

  ## Split dataframe into smaller dataframes by factor A (by rows) for aligning
  dataList <- split(dataFrame, as.factor(dataFrame[, 2]))

  ## Aligning by row median
  for (i in 1:I) { ## looping over the number of levels for factor A (I)
    dataList[[i]][, 1] <- dataList[[i]][, 1] - median(dataList[[i]][, 1])
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor B (by columns) for ranking
  dataList <- split(dataFrame, as.factor(dataFrame[, 3]))

  ## Rank by column
  for (j in 1:J) { ## looping over the number of levels for factor B (J)
    dataList[[j]][, 1] <- rank(dataList[[j]][, 1], ties.method = "average")
  }

  ## Rebind the (new) dataframe by row
  dataFrame <- do.call("rbind", dataList)

  ## Split dataframe into smaller dataframes by factor A (by row) and then by factor B (by column)
  ## for calculation
  dataList <- split(dataFrame, as.factor(dataFrame[, 2]))
  for (i in 1:I) {
    dataList[[i]] <- split(dataList[[i]], as.factor(dataList[[i]][, 3]))
  }

  ## Calculate all the I*(I-1)/2 values of V_ii'
  APCRCM <- c()
  for (i1 in (1:(I - 1))) {
    for (i2 in ((i1 + 1):I)) {
      Viip <- 0
      ## Calculate each cross-comparison
      for (j1 in (1:(J - 1))) {
        for (j2 in ((j1 + 1):J)) {
          for (k1 in 1:K) {
            for (k2 in 1:K) {
              for (k3 in 1:K) {
                for (k4 in 1:K) {
                  Viip <- Viip + (dataList[[i1]][[j1]][, 1][k1] + dataList[[i2]][[j2]][, 1][k2] - dataList[[i1]][[j2]][, 1][k3] - dataList[[i2]][[j1]][, 1][k4])^2
                }
              }
            }
          }
        }
      }
      APCRCM <- c(APCRCM, Viip)
    }
  }
  APCRCM <- max(APCRCM)
  ## Scale
  APCRCMD <- 2*APCRCM/(K^4*J*(J - 1))
  return(APCRCMD)
}

## p-value calculations -----
#' This function gives the approximate p-value for the test statistics APCSSA or APCSSM
#'
#' @param type The test type, either "APCSSA" or "APCSSM"
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param stat The calculated statistics using one of the two tests on a data set
#'
#' @returns The p-value of given statistics
#' @keywords internal
#' @noRd
# Function to calculate the p-value given a test stats and test type (combined)
.calc_p_value <- function(type, i, j, k, stat) {

  # Validate input type
  if (!type %in% c("APCSSA", "APCSSM")) {
    stop("Invalid type. Use 'APCSSA' or 'APCSSM'.")
  }

  # Construct the expected object name
  nullDist <- paste0("null", type, "_", i, "x", j, "x", k)

  # Load data if not already in memory
  if (!exists(nullDist, envir = .apc_cache)) {
    if (type == "APCSSA") {
      nullAPCSSA(i, j, k)
    } else {
      nullAPCSSM(i, j, k)
    }
  }

  # Check if loading was successful
  if (!exists(nullDist, envir = .apc_cache)) {
    stop("No null distribution readily available. Try running sim_", type, "_null(i, j, k).")
  }

  # Retrieve the null distribution
  nullvals <- get(nullDist, envir = .apc_cache)

  # Compute the approximate p-value (right-tailed)
  p_value <- mean(nullvals >= stat)

  return(p_value)
}

# Functions relevant to loading existing null distributions ----
#' Filter the existing null table for APCCRA and APCRCA based on I, J, K values for null mean and standard deviations
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @keywords internal
#' @noRd
.nullAPCXXA <- function(i, j, k) {
  # Check if the data set exists
  if (!exists("nullAPCXXA_summary")) stop("Error: nullAPCXXA_summary not found. Try running sim_nullAPCSSA().")

  # Extract column names
  col_names <- colnames(nullAPCXXA_summary)

  # Filter null table
  result <- nullAPCXXA_summary %>%
    dplyr::filter(
      .data[[col_names[1]]] == i,
      .data[[col_names[2]]] == j,
      .data[[col_names[3]]] == k
    )  %>%
    dplyr::select(-c(1:3))  # Remove first three columns

  # Return result
  if (nrow(result) == 0) {
    warning("No matching rows found.")
    return(NULL)
  }

  name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  assign(name, result, envir = .apc_cache)

  return(result)
}

#' Filter the existing null table for APCCRM and APCRCM based on I, J, K values for null mean and standard deviations
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @keywords internal
#' @noRd
.nullAPCXXM <- function(i, j, k) {
  # Check if the data set exists
  if (!exists("nullAPCXXM_summary")) stop("Error: nullAPCXXM_summary is not available in the package.")

  # Extract column names
  col_names <- colnames(nullAPCXXM_summary)

  # Filter null table
  result <- nullAPCXXM_summary %>%
    dplyr::filter(
      .data[[col_names[1]]] == i,
      .data[[col_names[2]]] == j,
      .data[[col_names[3]]] == k
    )  %>%
    dplyr::select(-c(1:3))  # Remove first three columns

  # Return result
  if (nrow(result) == 0) {
    warning("No matching rows found.")
    return(NULL)
  }

  name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  assign(name, result, envir = .apc_cache)

  return(result)
}

#' Load standardized null for APCSSA to find critical values and approximate p-values later
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return A numeric vector of null values of length 100,000
#' @keywords internal
#' @noRd
# Function to load APCSSA null data
.nullAPCSSA <- function(i, j, k) {
  # Construct the file name
  file_name <- paste0("APCSSA_Null_Distribution_", i, "x", j, "x", k, "_100kSim.RData")

  # Get full file path
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Ensure file exists
  if (!file.exists(file_path)) {
    stop("Error: File not found: ", file_name, ". Try running sim_nullAPCSSA().")
  }

  # Load the data into a temporary environment
  env <- new.env()
  loaded_objects <- load(file_path, envir = env)

  # Ensure at least one object was loaded
  if (length(loaded_objects) == 0) {
    stop("No objects found in the loaded file: ", file_name)
  }

  # Rename the first loaded object and assign it to the global environment
  new_name <- paste0("nullAPCSSA_", i, "x", j, "x", k)
  assign(new_name, env[[loaded_objects[1]]], envir = .apc_cache)

  message("Data successfully loaded into the global environment.")
}


#' Load standardized null for APCSSM to find critical values and approximate p-values later
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return A numeric vector of null values of length 100,000
#' @keywords internal
#' @noRd
# Function to load APCSSM data
.nullAPCSSM <- function(i, j, k) {
  # Construct the file name
  file_name <- paste0("APCSSM_Null_Distribution_", i, "x", j, "x", k, "_100kSim.RData")

  # Get full file path
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Ensure file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_name, ". Try running simulate_APCSSM_null(i, j, k).")
  }

  # Load the data into a temporary environment
  env <- new.env()
  loaded_objects <- load(file_path, envir = env)

  # Ensure at least one object was loaded
  if (length(loaded_objects) == 0) {
    stop("No objects found in the loaded file: ", file_name)
  }

  # Rename the first loaded object and assign it to the global environment
  new_name <- paste0("nullAPCSSM_", i, "x", j, "x", k)
  assign(new_name, env[[loaded_objects[1]]], envir = .apc_cache)

  message("Data successfully loaded into the global environment.")
}

# Functions relevant to simulating new null distributions ----

#' This function helps simulate the first null for APCSSA test statistics.
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores parallelly or not
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRA and APCRCA) that will get standardized into APCSSA.
#' @keywords internal
#' @noRd
#' @importFrom pbapply pblapply
# First Null
.null1APCSSA <- function(i, j, k,
                         numSim   = 100000,
                         parallel = TRUE,
                         verbose  = TRUE) {
  I <- i; J <- j; K <- k

  # turn pbapply bar on or off
  pbapply::pboptions(type = if (verbose) "timer" else "none")

  if (verbose) {
    message("-> [.null1APCSSA] Generating ", numSim, " null data sets ...")
  }

  if (parallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("I", "J", "K"), envir = environment())

    APCnull <- pbapply::pblapply(seq_len(numSim), function(n) {
      data.frame(
        value = stats::rnorm(I * J * K),
        A     = rep(1:I, each = K, times = J),
        B     = rep(1:J, each = I * K)
      )
    }, cl = cl)

  } else {
    APCnull <- pbapply::pblapply(seq_len(numSim), function(n) {
      data.frame(
        value = stats::rnorm(I * J * K),
        A     = rep(1:I, each = K, times = J),
        B     = rep(1:J, each = I * K)
      )
    })
  }

  if (verbose) {
    message("-> [.null1APCSSA] Computing APCCRA for each null ...")
  }
  if (parallel) {
    nullDistCRA <- unlist(
      pbapply::pblapply(APCnull, .APCCRAD, cl = cl),
      use.names = FALSE
    )
  } else {
    nullDistCRA <- unlist(
      pbapply::pblapply(APCnull, .APCCRAD),
      use.names = FALSE
    )
  }

  if (verbose) {
    message("-> [.null1APCSSA] Computing APCRCA for each null ...")
  }
  if (parallel) {
    nullDistRCA <- unlist(
      pbapply::pblapply(APCnull, .APCRCAD, cl = cl),
      use.names = FALSE
    )
  } else {
    nullDistRCA <- unlist(
      pbapply::pblapply(APCnull, .APCRCAD),
      use.names = FALSE
    )
  }

  # summarize, assign, save
  nullAPCXXA_summary <- data.frame(
    E_CRA = mean(nullDistCRA),
    SD_CRA = sd(nullDistCRA),
    E_RCA = mean(nullDistRCA),
    SD_RCA = sd(nullDistRCA)
  )
  name      <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  assign(name, nullAPCXXA_summary, envir = .apc_cache)
  save_path <- file.path(getwd(), paste0(name, ".RData"))
  save(list = name, file = save_path)
  if (verbose) message("-> [.null1APCSSA] Saved summary to ", save_path)

  nullAPCXXA_summary
}

#' This function helps simulate the second null for APCSSA test statistics.
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores parallelly or not
#'
#' @returns A numeric vector with length equals to the numSim of all APCSSA statistics on the null data sets
#' @keywords internal
#' @noRd
#' @importFrom pbapply pblapply
# Second Null
.null2APCSSA <- function(i, j, k,
                         numSim   = 100000,
                         parallel = TRUE,
                         verbose  = TRUE) {
  I <- i; J <- j; K <- k

  prev_name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  if (!exists(prev_name, envir = .apc_cache)) {
    stop("Error: run .null1APCSSA() first for this (i,j,k).")
  }
  APCSSnullDist <- get(prev_name, envir = .apc_cache)

  pbapply::pboptions(type = if (verbose) "timer" else "none")

  if (verbose) {
    message("-> [.null2APCSSA] Generating ", numSim, " null data sets ...")
  }

  if (parallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("APCSSA", ".APCCRAD", ".APCRCAD", "APCSSnullDist", "I", "J", "K"),
      envir = environment()
    )

    APCnull <- pbapply::pblapply(seq_len(numSim), function(n) {
      data.frame(
        value = stats::rnorm(I * J * K),
        A     = rep(1:I, each = K, times = J),
        B     = rep(1:J, each = I * K)
      )
    }, cl = cl)

  } else {
    APCnull <- pbapply::pblapply(seq_len(numSim), function(n) {
      data.frame(
        value = stats::rnorm(I * J * K),
        A     = rep(1:I, each = K, times = J),
        B     = rep(1:J, each = I * K)
      )
    })
  }

  if (verbose) {
    message("-> [.null2APCSSA] Computing APCSSA for each null ...")
  }
  if (parallel) {
    nullDistSSA <- unlist(
      pbapply::pblapply(APCnull, APCSSA, cl = cl),
      use.names = FALSE
    )
  } else {
    nullDistSSA <- unlist(
      pbapply::pblapply(APCnull, APCSSA),
      use.names = FALSE
    )
  }

  nameA     <- paste0("nullAPCSSA_", i, "x", j, "x", k)
  assign(nameA, nullDistSSA, envir = .apc_cache)
  save_path <- file.path(getwd(), paste0(nameA, ".RData"))
  save(list = nameA, file = save_path)
  if (verbose) message("-> [.null2APCSSA] Saved result to ", save_path)

  nullDistSSA
}

#' This function helps simulate the first null for APCSSM test statistics.
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores parallelly or not
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRM and APCRCM) that will get standardized into APCSSM.
#' @keywords internal
#' @noRd
#' @importFrom pbapply pblapply
# First Null
.null1APCSSM <- function(i, j, k,
                         numSim   = 100000,
                         parallel = TRUE,
                         verbose  = TRUE) {
  I <- i; J <- j; K <- k

  # toggle pbapply style
  pbapply::pboptions(type = if (verbose) "timer" else "none")

  if (verbose) {
    message("-> [.null1APCSSM] Generating ", numSim, " null data sets ...")
  }

  if (parallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = c("I", "J", "K"), envir = environment())

    APCnull <- pbapply::pblapply(
      seq_len(numSim),
      function(n) {
        data.frame(
          value = stats::rnorm(I * J * K),
          A     = rep(1:I, each = K, times = J),
          B     = rep(1:J, each = I * K)
        )
      },
      cl = cl
    )
  } else {
    APCnull <- pbapply::pblapply(
      seq_len(numSim),
      function(n) {
        data.frame(
          value = stats::rnorm(I * J * K),
          A     = rep(1:I, each = K, times = J),
          B     = rep(1:J, each = I * K)
        )
      }
    )
  }

  if (verbose) {
    message("-> [.null1APCSSM] Computing CRM for each null ...")
  }
  nullDistCRM <- unlist(
    if (parallel) {
      pbapply::pblapply(APCnull, .APCCRMD, cl = cl)
    } else {
      pbapply::pblapply(APCnull, .APCCRMD)
    },
    use.names = FALSE
  )

  if (verbose) {
    message("-> [.null1APCSSM] Computing RCM for each null ...")
  }
  nullDistRCM <- unlist(
    if (parallel) {
      pbapply::pblapply(APCnull, .APCRCMD, cl = cl)
    } else {
      pbapply::pblapply(APCnull, .APCRCMD)
    },
    use.names = FALSE
  )

  # summarize & save
  nullAPCXXM_summary <- data.frame(
    E_CRM = mean(nullDistCRM),
    SD_CRM = sd(nullDistCRM),
    E_RCM = mean(nullDistRCM),
    SD_RCM = sd(nullDistRCM)
  )
  name      <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  assign(name, nullAPCXXM_summary, envir = .apc_cache)
  save_path <- file.path(getwd(), paste0(name, ".RData"))
  save(list = name, file = save_path)

  if (verbose) {
    message("-> [.null1APCSSM] Saved result to: ", save_path)
  }

  nullAPCXXM_summary
}

#' This function helps simulate the second null for APCSSM test statistics.
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations required to generate this null distribution (the number of cases/data sets the tests are applied to)
#' @param parallel This is a Boolean option to run this simulation using multiple cores parallelly or not
#'
#' @returns A numeric vector with length equals to the numSim of all APCSSM statistics on the null data sets
#' @keywords internal
#' @noRd
#' @importFrom pbapply pblapply
# Second Null

.null2APCSSM <- function(i, j, k,
                         numSim   = 100000,
                         parallel = TRUE,
                         verbose  = TRUE) {
  I <- i; J <- j; K <- k

  prev_name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  if (!exists(prev_name, envir = .apc_cache)) {
    stop("Error: .null1APCSSM() must be run first for this (i,j,k).")
  }
  APCSSnullDist <- get(prev_name, envir = .apc_cache)

  pbapply::pboptions(type = if (verbose) "timer" else "none")

  if (verbose) {
    message("-> [.null2APCSSM] Generating ", numSim, " null data sets ...")
  }

  if (parallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("APCSSM", ".APCCRMD", ".APCRCMD", "APCSSnullDist", "I", "J", "K"),
      envir = environment()
    )

    APCnull <- pbapply::pblapply(
      seq_len(numSim),
      function(n) {
        data.frame(
          value = stats::rnorm(I * J * K),
          A     = rep(1:I, each = K, times = J),
          B     = rep(1:J, each = I * K)
        )
      },
      cl = cl
    )
  } else {
    APCnull <- pbapply::pblapply(
      seq_len(numSim),
      function(n) {
        data.frame(
          value = stats::rnorm(I * J * K),
          A     = rep(1:I, each = K, times = J),
          B     = rep(1:J, each = I * K)
        )
      }
    )
  }

  if (verbose) {
    message("-> [.null2APCSSM] Computing APCSSM for each null ...")
  }
  nullDistSSM <- unlist(
    if (parallel) {
      pbapply::pblapply(APCnull, APCSSM, cl = cl)
    } else {
      pbapply::pblapply(APCnull, APCSSM)
    },
    use.names = FALSE
  )

  nameA     <- paste0("nullAPCSSM_", i, "x", j, "x", k)
  assign(nameA, nullDistSSM, envir = .apc_cache)
  save_path <- file.path(getwd(), paste0(nameA, ".RData"))
  save(list = nameA, file = save_path)

  if (verbose) {
    message("-> [.null2APCSSM] Saved result to: ", save_path)
  }

  nullDistSSM
}

# Making sure the functions are global = seen everywhere
utils::globalVariables(c("nullAPCSSA", "nullAPCSSM"))
