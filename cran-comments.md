## R CMD check results

0 errors \| 0 warnings \| 0 notes

-   This is a resubmission. Changes based on reviewer feedback:
    -   Replaced direct `cat()`/`print()` calls in `APCSSA()` and `APCSSM()` with a proper S3 `summary.APC()` method. The test functions now return objects of class `"APC"` silently; users call `summary.APC()` to view results and plots.
    -   Removed default file-writing to working directory in `save_null()`. The `path` argument is now required (no default). Examples use `tempdir()`.
    -   Examples for `sim_nullAPCSSA()`, `sim_nullAPCSSM()`, and `save_null()` use `\dontrun{}` because these functions spawn parallel processes (via `parallel::makeCluster`), which is incompatible with CRAN's check environment that limits the number of cores.

## Note on `R/sysdata.rda`

`R/sysdata.rda` (5.6 KB) contains internal standardization parameters (means and standard deviations) for all 90 supported design configurations. These are required to compute the APCSSA and APCSSM test statistics without any external dependency. The companion package `APCinteractionData` is only needed for *p*-value estimation.

## Companion data package

Pre-computed null distributions (90 per test statistic, 100,000 simulations each) are provided by the companion package `APCinteractionData`, available on GitHub at <https://github.com/tranbaokhue/APCinteractionData>.

This package is listed in `Suggests` and is not required for core functionality. Two null distribution files are bundled in this package to support the examples. Users can also generate their own null distributions with `sim_nullAPCSSA()` / `sim_nullAPCSSM()`.

## Test environments

-   local macOS (aarch64-apple-darwin24.4.0), R 4.5.1
