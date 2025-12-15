#' Save null distribution to working directory
#'
#' @description
#' Saves a null distribution from the package cache to an RData file in the working directory.
#'
#' @param type Character string: either "APCSSA" or "APCSSM"
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param path Optional custom file path. If NULL (default), saves to working directory.
#'
#' @return Invisibly returns the file path where the data was saved
#'
#' @examples
#' \dontrun{
#' # After running simulation
#' sim_nullAPCSSA(2, 2, 2)
#'
#' # Save to working directory
#' save_null("APCSSA", 2, 2, 2)
#'
#' # Or save to custom location
#' save_null("APCSSA", 2, 2, 2, path = "~/my_nulls/")
#' }
#'
#' @seealso
#' \code{\link{sim_nullAPCSSA}}, \code{\link{sim_nullAPCSSM}}
#'
#' @export
save_null <- function(type, i, j, k, path = NULL) {
  if (!type %in% c("APCSSA", "APCSSM")) {
    stop("type must be 'APCSSA' or 'APCSSM'")
  }

  # Check for both null distribution types
  name_XXA <- paste0("nullAPCXX", substr(type, 6, 6), "_", i, "x", j, "x", k)
  name_SS <- paste0("null", type, "_", i, "x", j, "x", k)

  if (!exists(name_XXA, envir = .apc_cache)) {
    stop("Null distribution ", name_XXA, " not found in cache. Run sim_null", type, "() first.")
  }
  if (!exists(name_SS, envir = .apc_cache)) {
    stop("Null distribution ", name_SS, " not found in cache. Run sim_null", type, "() first.")
  }

  # Determine save directory
  save_dir <- if (is.null(path)) getwd() else path
  if (!dir.exists(save_dir)) {
    stop("Directory does not exist: ", save_dir)
  }

  # Save both files
  file_XXA <- file.path(save_dir, paste0(name_XXA, ".RData"))
  file_SS <- file.path(save_dir, paste0(name_SS, ".RData"))

  save(list = name_XXA, file = file_XXA, envir = .apc_cache)
  save(list = name_SS, file = file_SS, envir = .apc_cache)

  message("Saved:\n  ", file_XXA, "\n  ", file_SS)

  invisible(c(file_XXA, file_SS))
}
