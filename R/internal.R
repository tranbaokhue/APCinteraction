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


# Functions relevant to loading existing null distributions ----





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
.null1APCSSA <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  I <- i
  J <- j
  K <- k

  if (parallel) {
    # Create cluster
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(cl, varlist = c("I", "J", "K"), envir = environment())

    # Parallelized null matrix generation
    APCnull <- parallel::parLapply(cl, 1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distributions in parallel
    nullDistCRA <- unlist(parallel::parLapply(cl, APCnull, .APCCRAD), use.names = FALSE)
    nullDistRCA <- unlist(parallel::parLapply(cl, APCnull, .APCRCAD), use.names = FALSE)

    # Stop the cluster
    parallel::stopCluster(cl)

  } else {
    # Non-parallel version
    APCnull <- lapply(1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distributions sequentially
    nullDistCRA <- unlist(lapply(APCnull, .APCCRAD), use.names = FALSE)
    nullDistRCA <- unlist(lapply(APCnull, .APCRCAD), use.names = FALSE)
  }

  # Summarize results
  nullAPCXXA_summary <- data.frame(
    E_CRA = mean(nullDistCRA),
    SD_CRA = sd(nullDistCRA),
    E_RCA = mean(nullDistRCA),
    SD_RCA = sd(nullDistRCA)
  )

  # Assign to global environment with a formatted name
  name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  assign(name, nullAPCXXA_summary, envir = .GlobalEnv)

  # Save to the working directory
  save_path <- file.path(getwd(), paste0(name, ".RData"))
  save(list = name, file = save_path)

  message("Saved result to: ", save_path)

  return(nullAPCXXA_summary)
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
# Second Null
.null2APCSSA <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  I <- i
  J <- j
  K <- k

  # Retrieve the correct null distribution from null1APCSSA
  prev_name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  if (!exists(prev_name, envir = .GlobalEnv)) {
    stop("Error: The required null distribution from null1APCSSA does not exist. Run null1APCSSA() first.")
  }
  APCSSnullDist <- get(prev_name, envir = .GlobalEnv)

  if (parallel) {
    # Create and register cluster
    cl <- parallel::makeCluster(parallel::detectCores() - 1)

    # Export required functions and data
    parallel::clusterExport(cl, varlist = c("APCSSA", ".APCCRAD", ".APCRCAD", "APCSSnullDist", "I", "J", "K"), envir = environment())

    # Parallelized null matrix generation
    APCnull <- parallel::parLapply(cl, 1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distribution in parallel
    nullDistSSA <- unlist(parallel::parLapply(cl, APCnull, APCSSA), use.names = FALSE)

    # Stop the cluster
    parallel::stopCluster(cl)

  } else {
    # Non-parallel version
    APCnull <- lapply(1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distribution sequentially
    nullDistSSA <- unlist(lapply(APCnull, APCSSA), use.names = FALSE)
  }

  # Assign to global environment
  nameA <- paste0("nullAPCSSA_", i, "x", j, "x", k)
  assign(nameA, nullDistSSA, envir = .GlobalEnv)

  # Save to working directory
  save_path <- file.path(getwd(), paste0(nameA, ".RData"))
  save(list = nameA, file = save_path)

  message("Saved result to: ", save_path)

  return(nullDistSSA)
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
# First Null
.null1APCSSM <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  I <- i
  J <- j
  K <- k

  if (parallel) {
    # Create cluster
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    parallel::clusterExport(cl, varlist = c("I", "J", "K"), envir = environment())

    # Parallelized null matrix generation
    APCnull <- parallel::parLapply(cl, 1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distributions in parallel
    nullDistCRM <- unlist(parallel::parLapply(cl, APCnull, .APCCRMD), use.names = FALSE)
    nullDistRCM <- unlist(parallel::parLapply(cl, APCnull, .APCRCMD), use.names = FALSE)

    # Stop the cluster
    parallel::stopCluster(cl)

  } else {
    # Non-parallel version
    APCnull <- lapply(1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distributions sequentially
    nullDistCRM <- unlist(lapply(APCnull, .APCCRMD), use.names = FALSE)
    nullDistRCM <- unlist(lapply(APCnull, .APCRCMD), use.names = FALSE)
  }

  # Summarize results
  nullAPCXXM_summary <- data.frame(
    E_CRM = mean(nullDistCRM),
    SD_CRM = sd(nullDistCRM),
    E_RCM = mean(nullDistRCM),
    SD_RCM = sd(nullDistRCM)
  )

  # Assign to global environment with a formatted name
  name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  assign(name, nullAPCXXM_summary, envir = .GlobalEnv)

  # Save to the working directory
  save_path <- file.path(getwd(), paste0(name, ".RData"))
  save(list = name, file = save_path)

  message("Saved result to: ", save_path)

  return(nullAPCXXM_summary)
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
# Second Null
.null2APCSSM <- function(i, j, k, numSim = 100000, parallel = TRUE) {
  I <- i
  J <- j
  K <- k

  # Retrieve the correct null distribution from null1APCSSA
  prev_name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  if (!exists(prev_name, envir = .GlobalEnv)) {
    stop("Error: The required null distribution from null1APCSSM does not exist. Run null1APCSSM() first.")
  }
  APCSSnullDist <- get(prev_name, envir = .GlobalEnv)

  if (parallel) {
    # Create and register cluster
    cl <- parallel::makeCluster(parallel::detectCores() - 1)

    # Export required functions and data
    parallel::clusterExport(cl, varlist = c("APCSSM", ".APCCRMD", ".APCRCMD", "APCSSnullDist", "I", "J", "K"), envir = environment())

    # Parallelized null matrix generation
    APCnull <- parallel::parLapply(cl, 1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distribution in parallel
    nullDistSSM <- unlist(parallel::parLapply(cl, APCnull, APCSSM), use.names = FALSE)

    # Stop the cluster
    parallel::stopCluster(cl)

  } else {
    # Non-parallel version
    APCnull <- lapply(1:numSim, function(n) {
      data.frame(value = stats::rnorm(I * J * K),
                 A = rep(1:I, each = K, times = J),
                 B = rep(1:J, each = I * K))
    })

    # Compute null distribution sequentially
    nullDistSSM <- unlist(lapply(APCnull, APCSSM), use.names = FALSE)
  }

  # Assign to global environment
  nameA <- paste0("nullAPCSSM_", i, "x", j, "x", k)
  assign(nameA, nullDistSSM, envir = .GlobalEnv)

  # Save to working directory
  save_path <- file.path(getwd(), paste0(nameA, ".RData"))
  save(list = nameA, file = save_path)

  message("Saved result to: ", save_path)

  return(nullDistSSM)
}
