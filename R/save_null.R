#' Save null distribution to a directory
#'
#' @description
#' Saves a null distribution from the package's cache to RData files in the
#' specified directory.
#'
#' @param type Character string: either "APCSSA" or "APCSSM."
#' @param i The number of levels for factor A.
#' @param j The number of levels for factor B.
#' @param k The number of replications at each combination of Factor A and Factor B.
#' @param path Character string specifying the directory to save the files to.
#'
#' @return Invisibly returns the file paths to where the data were saved.
#'
#' @examples
#' \donttest{
#' # First generate the desired null distribution
#' sim_nullAPCSSA(2, 2, 2, 50000)
#'
#' # Save to a temporary directory
#' save_null("APCSSA", 2, 2, 2, path = tempdir())
#'
#' # Or save to a custom location
#' save_null("APCSSA", 2, 2, 2, path = "path/to/directory")
#' }
#'
#' @seealso
#' \code{\link{sim_nullAPCSSA}}, \code{\link{sim_nullAPCSSM}}
#'
#' @export
save_null <- function(type, i, j, k, path) {
  if (!type %in% c("APCSSA", "APCSSM")) {
    stop("type must be 'APCSSA' or 'APCSSM'")
  }

  # Check for both null distribution
  name_XXA <- paste0("nullAPCXX", substr(type, 6, 6), "_", i, "x", j, "x", k)
  name_SS <- paste0("null", type, "_", i, "x", j, "x", k)

  if (!exists(name_XXA, envir = .apc_cache)) {
    stop("Null distribution ", name_XXA, " not found in cache. Run sim_null", type, "() first.")
  }
  if (!exists(name_SS, envir = .apc_cache)) {
    stop("Null distribution ", name_SS, " not found in cache. Run sim_null", type, "() first.")
  }

  # Validate save directory
  if (!dir.exists(path)) {
    stop("Directory does not exist: ", path)
  }

  # Save both files
  file_XXA <- file.path(path, paste0(name_XXA, ".RData"))
  file_SS <- file.path(path, paste0(name_SS, ".RData"))

  save(list = name_XXA, file = file_XXA, envir = .apc_cache)
  save(list = name_SS, file = file_SS, envir = .apc_cache)

  message("Saved:\n  ", file_XXA, "\n  ", file_SS)

  invisible(c(file_XXA, file_SS))
}
