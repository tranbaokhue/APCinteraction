#' Filter the existing null table for APCSSA based on I, J, K values for null mean & sd
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

#' Filter the existing null table for APCSSM based on I, J, K values for null mean & sd
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



