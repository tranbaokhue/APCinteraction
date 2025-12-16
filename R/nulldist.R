# The following functions help user simulate their own null distributions
# to estimate the p-value depending on the 2-way ANOVA design they have.

#' Simulate the null distribution for APCSSA
#'
#' @description
#' A complete function that simulate the null distribution for APCSSA and automatically save it to the working directory
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations (data sets) to generate
#' @param parallel Boolean; if TRUE (the default), uses one fewer than the total number of cores for parallel computation.
#' @param verbose Logical; if TRUE (the default), prints a startup notice and progress bar.
#'
#'
#' @details
#' While we have provided null distributions for various settings in this package, there are still two-way settings that we don't readily have the null distribution to evaluate the significance of the test statistics (see **References** for the complete list of settings available). For example, one might wish to estimate the *p*-value with a higher `numSim`, ie. 250,000.
#'
#' @note
#' Benchmarks for 100,000 simulations (no bias toward any chip type):
#' - **3×3×3** (243 summands): ~ 5 hours on an Intel i5
#' - **2×5×5** (625 summands): < 5 hours on Apple M1 (7 cores, parallel)
#' - **6×6×5** (9375 summands): ~ 108 hours on Apple M1 (7 cores, parallel)
#' Factor levels (i,j) have some impact, but total summands scale as *k^4* (see **References** for more information),
#' so replication depth (k) dominates run time.
#'
#' **Tip:** leave `parallel = TRUE` to shorten run times dramatically.
#'
#' @references
#' Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024). Nonparametric tests for interaction in two-way ANOVA with balanced replications. *arXiv preprint* arXiv:2410.04700.
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRA and APCRCA) that will get standardized into APCSSA and a numeric vector with length equals to the `numSim` of all APCSSA statistics on the null data sets.
#'
#' @examples
#' \donttest{
#' sim_nullAPCSSA(2, 2, 2, 5000)  # This should take only a few seconds
#' }
#'
#' @seealso
#' \code{\link{APCSSA}}, \code{\link{sim_nullAPCSSM}}
#'
#' @importFrom magrittr %>%
#' @export
sim_nullAPCSSA <- function(i, j, k,
                           numSim   = 100000,
                           parallel = TRUE,
                           verbose  = TRUE) {
  if (verbose) {
    message(
      "Warning: sim_nullAPCSSA can take hours to run (depending on CPU & `parallel`).\n",
      "    See `?sim_nullAPCSSA` for details.\n"
    )
  }

  # 1) generate and summarize null #1
  .null1APCSSA(i, j, k,
               numSim   = numSim,
               parallel = parallel,
               verbose  = verbose)

  # 2) generate and save null #2
  .null2APCSSA(i, j, k,
               numSim   = numSim,
               parallel = parallel,
               verbose  = verbose)

  invisible(NULL)
}


#' Simulate the null distribution for APCSSM
#'
#' @description
#' A complete function that simulate the null distribution for APCSSM and automatically save it to the working directory
#'
#' @param i The number of levels in factor A
#' @param j The number of levels in factor B
#' @param k The number of observations at each level of factor A and B
#' @param numSim The number of simulations (data sets) to generate
#' @param parallel Boolean; if TRUE (the default), uses one fewer than the total number of cores for parallel computation.
#' @param verbose Logical; if TRUE (the default), prints a startup notice and progress bar.
#'
#' @details
#' While we have provided null distributions for various settings in this package, there are still two-way settings that we don't readily have the null distribution to evaluate the significance of the test statistics (see **References** for the complete list of settings available). For example, one might wish to estimate the *p*-value with a higher `numSim`, ie. 250,000.
#'
#' @note
#' Benchmarks for 100,000 simulations (no bias toward any chip type):
#' - **3×3×3** (243 summands): ~5 hours on an Intel i5
#' - **2×5×5** (625 summands): < 5 hours on Apple M1 (7 cores, parallel)
#' - **6×6×5** (9375 summands): ~108 hours on Apple M1 (7 cores, parallel)
#' Factor levels (i,j) have some impact, but total summands scale as *k^4* (see **References** for more information),
#' so replication depth (k) dominates run time.
#'
#' **Tip:** leave `parallel = TRUE` to shorten run times dramatically.
#'
#' @references
#' Tran, B. K., Wagaman, A. S., Nguyen, A., Jacobson, D., & Hartlaub, B. (2024). Nonparametric tests for interaction in two-way ANOVA with balanced replications. *arXiv preprint* arXiv:2410.04700.
#'
#' @returns A data frame with the null mean and standard deviation for the two test statistics (APCCRM and APCRCM) that will get standardized into APCSSM and a numeric vector with length equals to the `numSim` of all APCSSM statistics on the null data sets.
#'
#' @examples
#' \donttest{
#' sim_nullAPCSSA(2, 2, 2, 5000)  # This should take only a few seconds
#' }
#'
#' @seealso
#' \code{\link{APCSSM}}, \code{\link{sim_nullAPCSSA}}
#'
#' @importFrom magrittr %>%
#' @export
sim_nullAPCSSM <- function(i, j, k,
                           numSim   = 100000,
                           parallel = TRUE,
                           verbose  = TRUE) {
  if (verbose) {
    message(
      "Warning: sim_nullAPCSSM can take hours to run (depending on CPU & `parallel`).\n",
      "    See `?sim_nullAPCSSM` for details.\n"
    )
  }

  # 1) generate and summarize null #1
  .null1APCSSM(i, j, k,
               numSim   = numSim,
               parallel = parallel,
               verbose  = verbose)

  # 2) generate and save null #2
  .null2APCSSM(i, j, k,
               numSim   = numSim,
               parallel = parallel,
               verbose  = verbose)

  invisible(NULL)
}

