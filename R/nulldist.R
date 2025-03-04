#' Filter the existing null table for APCCRA and APCRCA based on I, J, K values for null mean and standard deviations
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @export
firstNullAPC_A <- function(i, j, k) {
  # Check if the data set exists
  if (!exists("FirstNullAPCSSA")) stop("Error: FirstNullAPCSSA is not available in the package.")

  # Extract column names
  col_names <- colnames(FirstNullAPCSSA)

  # Filter null table
  result <- FirstNullAPCSSA %>%
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

  return(result)
}

#' Filter the existing null table for APCCRM and APCRCM based on I, J, K values for null mean and standard deviations
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @export
firstNullAPC_M <- function(i, j, k) {
  # Check if the data set exists
  if (!exists("FirstNullAPCSSM")) stop("Error: FirstNullAPCSSM is not available in the package.")

  # Extract column names
  col_names <- colnames(FirstNullAPCSSM)

  # Filter null table
  result <- FirstNullAPCSSM %>%
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
secondNullAPCSSA <- function(i, j, k) {
  # Construct the file name
  file_name <- paste0("APCSSA Null Distribution ", i, "x", j, "x", k, "_100kSim.RData")

  # Get full file path
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Ensure file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_name)
  }

  # Load the data directly into the global environment
  load(file_path, envir = .GlobalEnv)

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
secondNullAPCSSM <- function(i, j, k) {
  # Construct the file name
  file_name <- paste0("APCSSM Null Distribution ", i, "x", j, "x", k, "_100kSim.RData")

  # Get full file path
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Ensure file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_name)
  }

  # Load the data directly into the global environment
  load(file_path, envir = .GlobalEnv)

  message("Data successfully loaded into the global environment.")
}

