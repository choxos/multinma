# Test file for Issue #53: Return Knot Locations from Spline Models

# Helper function to create minimal survival network for testing
create_test_survival_network <- function() {
  # Create minimal IPD survival data
  ipd_data <- data.frame(
    .study = rep(c("Study1", "Study2"), each = 50),
    .trt = rep(c("A", "B"), times = 50),
    .y = rep(NA, 100),  # Will be filled with Surv object
    time = c(rexp(50, 0.1), rexp(50, 0.08)),  # Event times
    status = rep(1, 100)  # All events for simplicity
  )
  
  # Create Surv object
  ipd_data$.y <- survival::Surv(ipd_data$time, ipd_data$status)
  
  # Create AgD survival data
  agd_data <- data.frame(
    .study = "Study3",
    .trt = "C",
    .y = list(survival::Surv(rexp(30, 0.12), rep(1, 30))),
    .sample_size = 30
  )
  
  # Create network
  network <- set_ipd(ipd_data, study = .study, trt = .trt, y = .y)
  network <- set_agd_surv(network, agd_data, study = .study, trt = .trt, y = .y, sample_size = .sample_size)
  
  return(network)
}

test_that("nma returns knot locations for mspline models", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Fit mspline model
  fit <- nma(network, 
             likelihood = "mspline",
             iter = 100,  # Small number for testing
             warmup = 50,
             chains = 2)
  
  # Check that knots component is present
  expect_true("knots" %in% names(fit))
  expect_true("basis" %in% names(fit))
  
  # Check knots structure
  expect_true(is.list(fit$knots))
  expect_true(all(sapply(fit$knots, is.numeric)))
  
  # Should have knots for each study
  expected_studies <- levels(network$ipd$.study)
  expect_equal(names(fit$knots), expected_studies)
  
  # All studies should have same number of knots
  knot_lengths <- sapply(fit$knots, length)
  expect_true(all(knot_lengths == knot_lengths[1]))
  
  # Knots should be in sorted order
  for (study_knots in fit$knots) {
    expect_equal(study_knots, sort(study_knots))
  }
})

test_that("nma returns knot locations for pexp models", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Fit piecewise exponential model
  fit <- nma(network,
             likelihood = "pexp", 
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Check that knots component is present
  expect_true("knots" %in% names(fit))
  expect_true("basis" %in% names(fit))
  
  # Check knots structure
  expect_true(is.list(fit$knots))
  expect_true(all(sapply(fit$knots, is.numeric)))
  
  # Should have knots for each study
  expected_studies <- levels(network$ipd$.study)
  expect_equal(names(fit$knots), expected_studies)
})

test_that("nma preserves user-specified knot configurations", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Define custom knots for each study
  custom_knots <- list(
    "Study1" = c(0, 2, 5, 8, 15),
    "Study2" = c(0, 1.5, 4, 7, 12),
    "Study3" = c(0, 3, 6, 10, 18)
  )
  
  # Fit model with custom knots
  fit <- nma(network,
             likelihood = "mspline",
             knots = custom_knots,
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Check that custom knots are preserved (but sorted)
  expect_equal(names(fit$knots), names(custom_knots))
  
  for (study in names(custom_knots)) {
    expect_equal(fit$knots[[study]], sort(custom_knots[[study]]))
  }
})

test_that("nma handles single vector knots specification", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Define single vector of knots to be shared
  shared_knots <- c(0, 2, 5, 8, 15)
  
  # Fit model with shared knots
  fit <- nma(network,
             likelihood = "mspline",
             knots = shared_knots,
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Check that knots are shared across all studies
  expected_studies <- levels(network$ipd$.study)
  expect_equal(names(fit$knots), expected_studies)
  
  # All studies should have the same knots (sorted)
  sorted_shared_knots <- sort(shared_knots)
  for (study in expected_studies) {
    expect_equal(fit$knots[[study]], sorted_shared_knots)
  }
})

test_that("nma knots work with aux_regression", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Add a covariate for regression
  network$ipd$age <- rnorm(nrow(network$ipd), 60, 10)
  network$agd_surv$age_mean <- 62
  network$agd_surv$age_sd <- 8
  
  # Define shared knots for aux_regression (required when aux_regression is used)
  shared_knots <- c(0, 3, 7, 12, 20)
  
  # Fit model with aux_regression
  fit <- nma(network,
             likelihood = "mspline",
             aux_regression = ~age,
             knots = shared_knots,
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Check that knots are preserved
  expect_true("knots" %in% names(fit))
  
  # Should have same knots for all studies when aux_regression is used
  sorted_shared_knots <- sort(shared_knots)
  for (study in names(fit$knots)) {
    expect_equal(fit$knots[[study]], sorted_shared_knots)
  }
})

test_that("nma does NOT return knots for non-spline models", {
  network <- create_test_survival_network()
  
  # Fit non-spline survival model
  fit <- nma(network,
             likelihood = "weibull",
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Should not have knots component
  expect_false("knots" %in% names(fit))
  expect_false("basis" %in% names(fit))
})

test_that("knots are consistent with basis component", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Fit mspline model
  fit <- nma(network,
             likelihood = "mspline",
             n_knots = 5,
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Check that knots and basis have same study names
  expect_equal(names(fit$knots), names(fit$basis))
  
  # Check that each study has appropriate number of knots
  # n_knots = 5 means 5 internal knots + 2 boundary knots = 7 total
  for (study in names(fit$knots)) {
    expect_equal(length(fit$knots[[study]]), 7)
  }
})

test_that("knots work with different n_knots settings", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Test with different numbers of knots
  for (n_knots in c(3, 5, 7)) {
    fit <- nma(network,
               likelihood = "mspline",
               n_knots = n_knots,
               iter = 100,
               warmup = 50,
               chains = 2)
    
    # Check that each study has correct number of total knots
    # n_knots internal + 2 boundary = n_knots + 2 total
    expected_total <- n_knots + 2
    
    for (study in names(fit$knots)) {
      expect_equal(length(fit$knots[[study]]), expected_total,
                   info = paste("Failed for n_knots =", n_knots))
    }
  }
})

test_that("knots functionality integrates with existing workflow", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Fit model and extract knots
  fit1 <- nma(network,
              likelihood = "mspline",
              iter = 100,
              warmup = 50,
              chains = 2)
  
  extracted_knots <- fit1$knots
  
  # Use extracted knots to fit another model
  fit2 <- nma(network,
              likelihood = "mspline",
              knots = extracted_knots,
              iter = 100,
              warmup = 50,
              chains = 2)
  
  # Knots should be identical (after sorting)
  expect_equal(fit1$knots, fit2$knots)
  
  # Should be able to use the model normally
  expect_true(inherits(fit2, "stan_nma"))
  expect_true("knots" %in% names(fit2))
})

test_that("class documentation includes knots component", {
  # This is more of a documentation test
  # Check that the stan_nma class has knots documented
  
  # Get the documentation for stan_nma class
  help_text <- capture.output(help("stan_nma-class", package = "multinma"))
  
  # Should mention knots in the documentation
  # This is a basic check - in practice, manual verification is better
  expect_true(any(grepl("knots", help_text, ignore.case = TRUE)))
})

test_that("knots handle edge cases gracefully", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Test with minimal knots (3 total: 1 internal + 2 boundary)
  fit <- nma(network,
             likelihood = "mspline", 
             n_knots = 1,
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Should still return knots
  expect_true("knots" %in% names(fit))
  
  # Each study should have 3 knots total
  for (study in names(fit$knots)) {
    expect_equal(length(fit$knots[[study]]), 3)
  }
})

test_that("knots are properly ordered by study factor levels", {
  skip_if_not_installed("splines2")
  
  network <- create_test_survival_network()
  
  # Ensure study factor has specific order
  network$ipd$.study <- factor(network$ipd$.study, levels = c("Study2", "Study1"))
  
  fit <- nma(network,
             likelihood = "mspline",
             iter = 100,
             warmup = 50,
             chains = 2)
  
  # Knots should be ordered according to factor levels
  expect_equal(names(fit$knots), levels(network$ipd$.study))
}) 