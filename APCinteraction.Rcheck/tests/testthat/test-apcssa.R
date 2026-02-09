test_that("APCSSA returns correct structure", {
  skip_on_cran()  # Requires null distribution

  # Generate null distribution first
  set.seed(111)
  sim_nullAPCSSA(2, 2, 2, numSim = 1000, parallel = FALSE, verbose = FALSE)

  # Create test data with interaction
  set.seed(222)
  nA <- 2; nB <- 2; nrep <- 2
  A_vals <- 1:nA
  B_vals <- 1:nB
  design <- expand.grid(A = A_vals, B = B_vals)
  design <- design[rep(seq_len(nrow(design)), each = nrep), ]

  mu <- A_vals[design$A] * B_vals[design$B]
  value <- mu + rnorm(length(mu))

  data <- data.frame(
    value = value,
    A = factor(design$A),
    B = factor(design$B)
  )

  # Run test
  result <- APCSSA(value ~ A + B, data = data)

  # Check structure
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Interaction", "Statistic", "p.value"))
  expect_equal(nrow(result), 1)

  # Check types
  expect_type(result$Interaction, "character")
  expect_type(result$Statistic, "double")
  expect_type(result$p.value, "double")

  # Check p-value is in [0,1]
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("APCSSA handles formula interface correctly", {
  skip_on_cran()

  # Ensure null distribution exists
  if (!exists("nullAPCXXA_2x2x2", envir = APCinteraction:::.apc_cache)) {
    set.seed(111)
    sim_nullAPCSSA(2, 2, 2, numSim = 1000, parallel = FALSE, verbose = FALSE)
  }

  set.seed(333)
  data <- data.frame(
    response = rnorm(8),
    factorA = factor(rep(1:2, each = 4)),
    factorB = factor(rep(1:2, times = 4))
  )

  result <- APCSSA(response ~ factorA + factorB, data = data)

  # Should work and return interaction name
  expect_match(result$Interaction, "factorA:factorB")
})

test_that("APCSSA produces exact expected values from documentation example", {
  skip_on_cran()  # Requires pre-computed null distribution

  # Example from documentation
  set.seed(206)

  # Parameters - number of levels for factors A and B, and replications
  nA <- 3; nB <- 3; nrep <- 3

  # Generate levels and create the full design
  A_vals <- 1:nA
  B_vals <- 1:nB
  design <- expand.grid(A = A_vals, B = B_vals)
  design <- design[rep(seq_len(nrow(design)), each = nrep), ]

  # Compute response as interaction (product of numeric factor levels)
  mu <- A_vals[design$A] * B_vals[design$B]
  value <- mu + rnorm(length(mu))  # Add normal noise

  data <- data.frame(
    value = value,
    A = factor(design$A),
    B = factor(design$B)
  )

  # Run the APCSSA test
  result <- APCSSA(value ~ A + B, data = data)

  # Check exact expected values
  expect_equal(result$Statistic, 4.084583, tolerance = 1e-6)
  expect_equal(result$p.value, 0.00263, tolerance = 1e-5)
})

test_that("APCSSA errors with incorrect formula", {
  data <- data.frame(
    value = rnorm(8),
    A = factor(rep(1:2, each = 4)),
    B = factor(rep(1:2, times = 4))
  )

  # Wrong formula (only one factor)
  expect_error(APCSSA(value ~ A, data = data))

  # Wrong formula (three factors)
  data$C <- factor(rep(1:2, 4))
  expect_error(APCSSA(value ~ A + B + C, data = data))
})
