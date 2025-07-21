#' APCinteraction: All possible crossed comparisons (APC) nonparametric tests of interaction in balanced two-way ANOVA models
#'
#' Provides nonparametric tests for interaction in two-way ANOVA designs with balanced
#' replications using all possible comparisons. The `APCSSA` and `APCSSM` statistics
#' extend previous methods and allow greater flexibility and power in detecting interactions
#' when the data have non-normal errors. The package includes tools for computing these
#' test statistics, generating interaction plots, and—when necessary—simulating their null
#' distributions.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{APCSSA}}}{Compute the APCSSA interaction test statistic.}
#'   \item{\code{\link{sim_nullAPCSSA}}}{Simulate and save the null distribution for the APCSSA statistic.}
#'   \item{\code{\link{APCSSM}}}{Compute the APCSSM interaction test statistic.}
#'   \item{\code{\link{sim_nullAPCSSM}}}{Simulate and save the null distribution for the APCSSM statistic.}
#' }
#'
#' @docType _PACKAGE
#' @name APCinteraction
#' @author Bao Khue Tran
#' @keywords internal
NULL
