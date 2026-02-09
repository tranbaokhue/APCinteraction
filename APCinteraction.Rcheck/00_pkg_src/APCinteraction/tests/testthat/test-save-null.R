test_that("save_null creates RData files", {
  skip_on_cran()

  # Generate null distribution
  set.seed(777)
  sim_nullAPCSSA(2, 2, 2, numSim = 100, parallel = FALSE, verbose = FALSE)

  # Create temp directory for test
  temp_dir <- tempdir()

  # Save
  result <- save_null("APCSSA", 2, 2, 2, path = temp_dir)

  # Check files were created
  expect_true(file.exists(file.path(temp_dir, "nullAPCXXA_2x2x2.RData")))
  expect_true(file.exists(file.path(temp_dir, "nullAPCSSA_2x2x2.RData")))

  # Check return value
  expect_length(result, 2)
  expect_true(all(grepl("\\.RData$", result)))

  # Clean up
  unlink(file.path(temp_dir, "nullAPCXXA_2x2x2.RData"))
  unlink(file.path(temp_dir, "nullAPCSSA_2x2x2.RData"))
})

test_that("save_null errors when null distribution not in cache", {
  # Try to save non-existent distribution
  expect_error(
    save_null("APCSSA", 99, 99, 99),
    "not found in cache"
  )
})

test_that("save_null errors with invalid type", {
  expect_error(
    save_null("INVALID", 2, 2, 2),
    "must be 'APCSSA' or 'APCSSM'"
  )
})

test_that("save_null errors with non-existent directory", {
  skip_on_cran()

  # Ensure distribution exists
  if (!exists("nullAPCXXA_2x2x2", envir = APCinteraction:::.apc_cache)) {
    set.seed(777)
    sim_nullAPCSSA(2, 2, 2, numSim = 100, parallel = FALSE, verbose = FALSE)
  }

  expect_error(
    save_null("APCSSA", 2, 2, 2, path = "/this/path/does/not/exist"),
    "Directory does not exist"
  )
})
