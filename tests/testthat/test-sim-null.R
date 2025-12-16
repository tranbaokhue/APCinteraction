test_that("sim_nullAPCSSA generates expected null distributions for 2x2x2", {
  skip_on_cran()  # The simulation can take long

  set.seed(12345)

  # Run with small numSim for testing - 8 seconds to run on a Mac M4 core
  result <- sim_nullAPCSSA(2, 2, 2, numSim = 5000, parallel = FALSE, verbose = FALSE)

  # Check that null distributions were created in cache
  expect_true(exists("nullAPCXXA_2x2x2", envir = APCinteraction:::.apc_cache))
  expect_true(exists("nullAPCSSA_2x2x2", envir = APCinteraction:::.apc_cache))

  # Check structure of nullAPCXXA
  nullXXA <- get("nullAPCXXA_2x2x2", envir = APCinteraction:::.apc_cache)
  expect_s3_class(nullXXA, "data.frame")
  expect_named(nullXXA, c("E_CRA", "SD_CRA", "E_RCA", "SD_RCA"))
  expect_equal(nrow(nullXXA), 1)

  # Check structure of nullAPCSSA
  nullSSA <- get("nullAPCSSA_2x2x2", envir = APCinteraction:::.apc_cache)
  expect_type(nullSSA, "double")
  expect_length(nullSSA, 5000)

  # All values should be numeric and finite
  expect_true(all(is.finite(nullSSA)))
})

test_that("sim_nullAPCSSM generates expected null distributions for 2x2x2", {
  skip_on_cran()  # Too slow for CRAN checks

  set.seed(12345)

  # Run with small numSim for testing
  result <- sim_nullAPCSSM(2, 2, 2, numSim = 5000, parallel = FALSE, verbose = FALSE)

  # Check that null distributions were created in cache
  expect_true(exists("nullAPCXXM_2x2x2", envir = APCinteraction:::.apc_cache))
  expect_true(exists("nullAPCSSM_2x2x2", envir = APCinteraction:::.apc_cache))

  # Check structure of nullAPCXXM
  nullXXM <- get("nullAPCXXM_2x2x2", envir = APCinteraction:::.apc_cache)
  expect_s3_class(nullXXM, "data.frame")
  expect_named(nullXXM, c("E_CRM", "SD_CRM", "E_RCM", "SD_RCM"))
  expect_equal(nrow(nullXXM), 1)

  # Check structure of nullAPCSSM
  nullSSM <- get("nullAPCSSM_2x2x2", envir = APCinteraction:::.apc_cache)
  expect_type(nullSSM, "double")
  expect_length(nullSSM, 5000)

  # All values should be numeric and finite
  expect_true(all(is.finite(nullSSM)))
})

test_that("sim_nullAPCSSA produces reproducible results with same seed", {
  skip_on_cran()

  # Clear cache
  if (exists("nullAPCXXA_2x2x2", envir = APCinteraction:::.apc_cache)) {
    rm(list = "nullAPCXXA_2x2x2", envir = APCinteraction:::.apc_cache)
  }
  if (exists("nullAPCSSA_2x2x2", envir = APCinteraction:::.apc_cache)) {
    rm(list = "nullAPCSSA_2x2x2", envir = APCinteraction:::.apc_cache)
  }

  set.seed(999)
  sim_nullAPCSSA(2, 2, 2, numSim = 100, parallel = FALSE, verbose = FALSE)
  result1 <- get("nullAPCSSA_2x2x2", envir = APCinteraction:::.apc_cache)

  # Clear cache again
  rm(list = c("nullAPCXXA_2x2x2", "nullAPCSSA_2x2x2"), envir = APCinteraction:::.apc_cache)

  set.seed(999)
  sim_nullAPCSSA(2, 2, 2, numSim = 100, parallel = FALSE, verbose = FALSE)
  result2 <- get("nullAPCSSA_2x2x2", envir = APCinteraction:::.apc_cache)

  # Should be identical
  expect_equal(result1, result2)
})
