#' APCinteraction: All Possible Comparisons (APC) Nonparametric Tests for Interaction in Balanced Two-way ANOVA Models
#'
#' Provides nonparametric tests for interaction in two-way ANOVA designs with balanced
#' replications using All Possible Comparisons (APC). The `APCSSA` and `APCSSM` statistics
#' extend previous methods, allow greater flexibility, and demonstrate higher power in detecting interactions
#' for non-normal data. The package includes optimized functions for computing these
#' test statistics, generating interaction plots, and simulating their null
#' distributions.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{APCSSA}}}{Compute the APCSSA interaction test statistic.}
#'   \item{\code{\link{sim_nullAPCSSA}}}{Simulate the null distribution for the APCSSA statistic.}
#'   \item{\code{\link{APCSSM}}}{Compute the APCSSM interaction test statistic.}
#'   \item{\code{\link{sim_nullAPCSSM}}}{Simulate the null distribution for the APCSSM statistic.}
#'
#'   \item{\code{\link{save_null}}}{Save the simulated null distributions locally for future usage.}
#' }
#'
#' @name APCinteraction
#' @author
#' Bao Khue Tran (maintainer) \email{baokhuetran@outlook.com}, Andrew Nguyen, Amy Wagaman, and Bradley Hartlaub
#'
#' @keywords interaction nonparametric anova

"_PACKAGE"

#' @importFrom stats median sd model.frame setNames interaction.plot
#' @importFrom grDevices hcl.colors
#' @importFrom rlang .data
NULL
