#' Filter the existing null table for APCCRA and APCRCA based on I, J, K values for null mean and standard deviations
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @export
nullAPCXXA <- function(i, j, k) {
  # Check if the data set exists
  if (!exists("nullAPCXXA_summary")) stop("Error: nullAPCXXA_summary is not available in the package.")

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

  obj_name <- paste0("nullAPCXXA_", i, "x", j, "x", k)
  assign(obj_name, result, envir = .GlobalEnv)

  return(result)
}

#' Filter the existing null table for APCCRM and APCRCM based on I, J, K values for null mean and standard deviations
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @export
nullAPCXXM <- function(i, j, k) {
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

  obj_name <- paste0("nullAPCXXM_", i, "x", j, "x", k)
  assign(obj_name, result, envir = .GlobalEnv)

  return(result)
}

#' Load standardized null for APCSSA to find critical values and approximate p-values later
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return A numeric vector of null values of length 100,000
#' @export
# Function to load APCSSA null data
nullAPCSSA <- function(i, j, k) {
  # Construct the file name
  file_name <- paste0("APCSSA Null Distribution ", i, "x", j, "x", k, "_100kSim.RData")

  # Get full file path
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Ensure file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_name, ". Try running simulate_APCSSA_null(i, j, k).")
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
  assign(new_name, env[[loaded_objects[1]]], envir = .GlobalEnv)

  message("Data successfully loaded into the global environment.")
}


#' Load standardized null for APCSSM to find critical values and approximate p-values later
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return A numeric vector of null values of length 100,000
#' @export
# Function to load APCSSM data
nullAPCSSM <- function(i, j, k) {
  # Construct the file name
  file_name <- paste0("APCSSM Null Distribution ", i, "x", j, "x", k, "_100kSim.RData")

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
  assign(new_name, env[[loaded_objects[1]]], envir = .GlobalEnv)

  message("Data successfully loaded into the global environment.")
}

#' This function gives the approximate p-value for the test statistics APCSSA or APCSSM
#'
#' @param type The test type, either "APCSSA" or "APCSSM"
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param stat The calculated statistics using one of the two tests on a data set
#'
#' @returns The p-value of given statistics
#' @export
# Function to calculate the p-value given a test stats and test type (combined)
calc_p_value <- function(type, i, j, k, stat) {

  # Validate input type
  if (!type %in% c("APCSSA", "APCSSM")) {
    stop("Invalid type. Use 'APCSSA' or 'APCSSM'.")
  }

  # Construct the expected object name
  nullDist <- paste0("null", type, "_", i, "x", j, "x", k)

  # Load data if not already in memory
  if (!exists(nullDist, envir = .GlobalEnv)) {
    if (type == "APCSSA") {
      nullAPCSSA(i, j, k)
    } else {
      nullAPCSSM(i, j, k)
    }
  }

  # Check if loading was successful
  if (!exists(nullDist, envir = .GlobalEnv)) {
    stop("No null distribution readily available. Try running sim_", type, "_null(i, j, k).")
  }

  # Retrieve the null distribution
  nullvals <- get(nullDist, envir = .GlobalEnv)

  # Compute the approximate p-value (right-tailed)
  p_value <- mean(nullvals >= stat)

  return(p_value)
}


#



sim_nullAPCSSA <- function(i, j, k){

}

# First Null

#### Null Matrices ----
cl <- parallel::makeCluster(detectCores()-1)
registerDoParallel(cl)
# Parallelized foreach loop
I <- Ai
J <- Bj
K <- NumReps
numTrial <- 100000
APCnull<- NULL
APCnull <- foreach(i = 1:numTrial) %dopar% {
  nullData <- data.frame(value = rnorm(I * J * K),
                         A = rep(1:I, each = K, times = J),
                         B = rep(1:J, each = I * K)
  )
  APCnull[[i]] <- nullData
}
# Stop the cluster
stopCluster(cl)

#### Null Mean & Variance ----
cl <- parallel::makeCluster(detectCores()-1)
registerDoParallel(cl)
nullDistCRA <- unlist(parLapply(cl, APCnull,APCCRADts), use.names = FALSE)
nullDistRCA <- unlist(parLapply(cl, APCnull,APCRCADts), use.names = FALSE)
nullDistCRM <- unlist(parLapply(cl, APCnull,APCCRMDts), use.names = FALSE)
nullDistRCM <- unlist(parLapply(cl, APCnull,APCRCMDts), use.names = FALSE)
stopCluster(cl)

summary<- c(mean(nullDistCRA), sd(nullDistCRA), mean(nullDistRCA), sd(nullDistRCA), mean(nullDistCRM), sd(nullDistCRM), mean(nullDistRCM), sd(nullDistRCM))

## Means and standard deviations of CRA, RCA, CRM, RCM
APCSSnullDist_AixBjxK<-list(nullDistCRA, nullDistRCA, nullDistCRM, nullDistRCM, summary)
save(APCSSnullDist_AixBjxK,file="APCSSnullDist_AixBjxK_100kSim.RData")

APCSSNull_results_AixBjxK<-tibble(Stats=c("CRA Mean", "CRA sd","RCA Mean", "RCA sd","CRM Mean","CRM sd","RCM Mean","RCM sd"))%>%
  mutate(Value=APCSSnullDist_AixBjxK[[5]])

## 100k Sim ----
nullDistAPCSSA <-NULL
nullDistAPCSSM <-NULL
APCSSnullDist <- NULL
APCSSnullDist_AixBjxK[[5]]
APCSSnullDist<- c()

cl <- parallel::makeCluster(detectCores()-1)
registerDoParallel(cl)
I <- Ai
J <- Bj
K <- NumReps
numTrial <- 100000
APCSSAnull <- NULL
APCSSMnull <- NULL
clusterExport(cl,list("APCCRADts","APCRCADts","APCSSAts","APCCRMDts",
                      "APCRCMDts","APCSSMts","APCSSnullDist"))

# Second null for Average
APCSSAnull <- foreach(i = 1:numTrial) %dopar% {
  nullData <- data.frame(value = rnorm(I * J * K),
                         A = rep(1:I, each = K, times = J),
                         B = rep(1:J, each = I * K)
  )
  APCSSAnull[[i]] <- nullData
}

nullDistAPCSSA <- unlist(parLapply(cl, APCSSAnull,APCSSAts), use.names = FALSE)
save(nullDistAPCSSA,file="APCSSA Null Distribution AixBjxK_100kSim.RData")
stopCluster(cl)

# Second null for Median
cl <- parallel::makeCluster(detectCores()-1)
registerDoParallel(cl)
I <- Ai
J <- Bj
K <- NumReps
numTrial <- 100000
APCSSAnull <- NULL
APCSSMnull <- NULL
clusterExport(cl,list("APCCRADts","APCRCADts","APCSSAts","APCCRMDts",
                      "APCRCMDts","APCSSMts","APCSSnullDist"))

APCSSMnull <- foreach(i = 1:numTrial) %dopar% {
  nullData <- data.frame(value = rnorm(I * J * K),
                         A = rep(1:I, each = K, times = J),
                         B = rep(1:J, each = I * K)
  )
  APCSSMnull[[i]] <- nullData
}

nullDistAPCSSM <- unlist(parLapply(cl, APCSSMnull,APCSSMts), use.names = FALSE)
save(nullDistAPCSSM,file="APCSSM Null Distribution AixBjxK_100kSim.RData")
stopCluster(cl)



