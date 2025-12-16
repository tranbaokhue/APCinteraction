test_that("APCSSM returns correct structure", {
  skip_on_cran()  # Requires null distribution

  # Generate null distribution first
  set.seed(444)
  sim_nullAPCSSM(2, 2, 2, numSim = 1000, parallel = FALSE, verbose = FALSE)

  # Create test data with interaction
  set.seed(555)
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
  result <- APCSSM(value ~ A + B, data = data)

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

test_that("APCSSM handles formula interface correctly", {
  skip_on_cran()

  # Ensure null distribution exists
  if (!exists("nullAPCXXM_2x2x2", envir = APCinteraction:::.apc_cache)) {
    set.seed(444)
    sim_nullAPCSSM(2, 2, 2, numSim = 1000, parallel = FALSE, verbose = FALSE)
  }

  set.seed(666)
  data <- data.frame(
    response = rnorm(8),
    factorA = factor(rep(1:2, each = 4)),
    factorB = factor(rep(1:2, times = 4))
  )

  result <- APCSSM(response ~ factorA + factorB, data = data)

  # Should work and return interaction name
  expect_match(result$Interaction, "factorA:factorB")
})

test_that("APCSSM produces exact expected values from documentation example", {
  skip_on_cran()  # Requires pre-computed null distribution

  # Example from documentation
  set.seed(36)

  # Parameters - number of levels for factors A and B, replications, and interaction effect
  nA <- 3; nB <- 4; nrep <- 2; c <- 1.25

  # Generate levels and create the full design
  A_vals <- seq(-2, 2, length.out = nA)
  B_vals <- seq(-1.5, 1.5, length.out = nB)
  design <- expand.grid(A = 1:nA, B = 1:nB)
  design <- design[rep(seq_len(nrow(design)), each = nrep), ]

  # Create the specific interaction matrix: alternate ±c in top rows
  specInt <- matrix(0, nA, nB)
  specInt[1:2, 1:nB] <- matrix(c(c, -c, -c, c), nrow = 2, byrow = TRUE)

  # Compute response
  mu <- A_vals[design$A] + B_vals[design$B] + specInt[cbind(design$A, design$B)]
  value <- mu + rt(length(mu), df = 1)

  data <- data.frame(
    value = value,
    A = factor(design$A, labels = round(A_vals, 2)),
    B = factor(design$B, labels = round(B_vals, 2))
  )

  # Run the APCSSM test
  result <- APCSSM(value ~ A + B, data = data)

  # Check exact expected values
  expect_equal(result$Statistic, 2.483915, tolerance = 1e-6)
  expect_equal(result$p.value, 0.02707, tolerance = 1e-5)
})

test_that("APCSSM errors with incorrect formula", {
  data <- data.frame(
    value = rnorm(8),
    A = factor(rep(1:2, each = 4)),
    B = factor(rep(1:2, times = 4))
  )

  # Wrong formula (only one factor)
  expect_error(APCSSM(value ~ A, data = data))

  # Wrong formula (three factors)
  data$C <- factor(rep(1:2, 4))
  expect_error(APCSSM(value ~ A + B + C, data = data))
})
