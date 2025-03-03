#' Filter the existing null table for APCCRA and APCRCA based on I, J, K values for null mean and standard deviations
#'
#' @param I The number of levels in factor A
#' @param J The number of levels in factor B
#' @param K The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @export
filter_nullAPCSSA <- function(I1, J1, K1) {
  # Check if the data set exists
  if (!exists("FirstNullAPCSSA")) stop("Error: FirstNullAPCSSA is not available in the package.")

  # Extract column names
  col_names <- colnames(FirstNullAPCSSA)

  # Filter null table
  result <- FirstNullAPCSSA %>%
    dplyr::filter(
      .data[[col_names[1]]] == I1,
      .data[[col_names[2]]] == J1,
      .data[[col_names[3]]] == K1
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
#' @param I The number of levels in factor A
#' @param J The number of levels in factor B
#' @param K The number of observations at each level of factor A and B
#' @return The null mean and standard deviation of tests APCCRA and APCRCA
#' @export
filter_nullAPCSSM <- function(I1, J1, K1) {
  # Check if the data set exists
  if (!exists("FirstNullAPCSSM")) stop("Error: FirstNullAPCSSM is not available in the package.")

  # Extract column names
  col_names <- colnames(FirstNullAPCSSM)

  # Filter null table
  result <- FirstNullAPCSSM %>%
    dplyr::filter(
      .data[[col_names[1]]] == I1,
      .data[[col_names[2]]] == J1,
      .data[[col_names[3]]] == K1
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
#' @param I The number of levels in factor A
#' @param J The number of levels in factor B
#' @param K The number of observations at each level of factor A and B
#' @return A numeric vector of null values of length 100,000
#' @export
# Function to load APCSSA null data
load_APCSSAnull <- function(I, J, K) {
  # Construct the file name based on I, J, K
  file_name <- paste0("APCSSA Null Distribution ", I, "x", J, "x", K, "_100kSim.RData")

  # Construct the full path to the file in extdata directory
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Check if the file exists before trying to load
  if (file.exists(file_path)) {
    # Load the data
    load(file_path)
    message("Loaded data from ", file_name)
  } else {
    stop("File not found: ", file_name)
  }
}

#' Load standardized null for APCSSM to find critical values and approximate p-values later
#'
#' @param I The number of levels in factor A
#' @param J The number of levels in factor B
#' @param K The number of observations at each level of factor A and B
#' @return A numeric vector of null values of length 100,000
#' @export
# Function to load APCSSM data
load_APCSSMnull <- function(I, J, K) {
  # Construct the file name based on I, J, K
  file_name <- paste0("APCSSM Null Distribution ", I, "x", J, "x", K, "_100kSim.RData")

  # Construct the full path to the file in extdata directory
  file_path <- system.file("extdata", file_name, package = "APCinteraction")

  # Check if the file exists before trying to load
  if (file.exists(file_path)) {
    # Load the data
    load(file_path)
    message("Loaded data from ", file_name)
  } else {
    stop("File not found: ", file_name)
  }
}



