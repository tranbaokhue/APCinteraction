## APC CRA ----
.APCCRADts <- function(dataFrame) {
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

## APC RCA ----
.APCRCADts <- function(dataFrame) {
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

## APC CRM ----
.APCCRMDts <- function(dataFrame) {
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


## APC RCM ----
.APCRCMDts <- function(dataFrame) {
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
