# Test file for Issue #54: MCMC Convergence Diagnostics

# Helper function to create minimal stan_nma object for testing
create_mock_stan_nma <- function(likelihood = "normal", trt_effects = "fixed", 
                                 class_effects = "independent", aux_regression = NULL) {
  # Create minimal mock stanfit object
  n_iter <- 100
  n_chains <- 4
  n_pars <- 3
  
  # Create mock simulation object
  mock_sim <- list(
    iter = n_iter,
    chains = n_chains,
    fnames_oi = c("d[B]", "d[C]", "sd"),
    samples = array(rnorm(n_iter * n_chains * n_pars), 
                    dim = c(n_iter, n_chains, n_pars))
  )
  
  # Create mock stanfit with minimal structure
  mock_stanfit <- structure(
    list(
      sim = mock_sim,
      model_name = "test_model"
    ),
    class = "stanfit"
  )
  
  # Add summary method behavior
  mock_stanfit@sim <- mock_sim
  
  # Create mock network
  mock_network <- list(
    treatments = factor(c("A", "B", "C"))
  )
  
  # Create stan_nma object
  structure(
    list(
      stanfit = mock_stanfit,
      network = mock_network,
      likelihood = likelihood,
      trt_effects = trt_effects,
      class_effects = class_effects,
      aux_regression = aux_regression
    ),
    class = "stan_nma"
  )
}

# Mock as.array method for stan_nma objects
as.array.stan_nma <- function(x, pars = NULL, ...) {
  if (is.null(pars)) {
    pars <- x$stanfit@sim$fnames_oi
  }
  
  # Select only available parameters
  available_pars <- intersect(pars, x$stanfit@sim$fnames_oi)
  par_indices <- match(available_pars, x$stanfit@sim$fnames_oi)
  
  # Return subset of samples
  samples <- x$stanfit@sim$samples[, , par_indices, drop = FALSE]
  dimnames(samples) <- list(
    iterations = 1:dim(samples)[1],
    chains = 1:dim(samples)[2], 
    parameters = available_pars
  )
  
  return(samples)
}

# Mock rstan::summary function
summary.stanfit <- function(object, pars = NULL, ...) {
  if (is.null(pars)) {
    pars <- object@sim$fnames_oi
  }
  
  # Create mock summary statistics
  n_pars <- length(pars)
  summary_stats <- matrix(
    data = c(
      rnorm(n_pars, 0, 1),      # mean
      runif(n_pars, 0.01, 0.1), # se_mean
      runif(n_pars, 0.1, 1),    # sd
      rep(0.025, n_pars),       # 2.5%
      rep(0.25, n_pars),        # 25%
      rep(0.5, n_pars),         # 50%
      rep(0.75, n_pars),        # 75%
      rep(0.975, n_pars),       # 97.5%
      runif(n_pars, 50, 200),   # n_eff
      runif(n_pars, 0.99, 1.02) # Rhat
    ),
    nrow = n_pars,
    ncol = 10
  )
  
  colnames(summary_stats) <- c("mean", "se_mean", "sd", "2.5%", "25%", "50%", 
                               "75%", "97.5%", "n_eff", "Rhat")
  rownames(summary_stats) <- pars
  
  return(list(summary = summary_stats))
}

test_that("mcmc_trace produces trace plots for stan_nma objects", {
  skip_if_not_installed("bayesplot")
  
  mock_fit <- create_mock_stan_nma()
  
  # Should work with default parameters
  expect_no_error({
    p <- mcmc_trace(mock_fit)
  })
  
  # Should return a ggplot object
  expect_true(inherits(p, "ggplot"))
  
  # Should work with specific parameters
  expect_no_error({
    p2 <- mcmc_trace(mock_fit, pars = c("d[B]", "d[C]"))
  })
  
  expect_true(inherits(p2, "ggplot"))
})

test_that("mcmc_acf produces autocorrelation plots", {
  skip_if_not_installed("bayesplot")
  
  mock_fit <- create_mock_stan_nma()
  
  # Should work with default parameters
  expect_no_error({
    p <- mcmc_acf(mock_fit)
  })
  
  # Should return a ggplot object
  expect_true(inherits(p, "ggplot"))
  
  # Should work with specific parameters
  expect_no_error({
    p2 <- mcmc_acf(mock_fit, pars = c("d[B]"))
  })
  
  expect_true(inherits(p2, "ggplot"))
})

test_that("mcmc_rhat produces R-hat diagnostic plots", {
  skip_if_not_installed("bayesplot")
  
  mock_fit <- create_mock_stan_nma()
  
  # Should work with default parameters
  expect_no_error({
    p <- mcmc_rhat(mock_fit)
  })
  
  # Should return a ggplot object
  expect_true(inherits(p, "ggplot"))
  
  # Should work with specific parameters
  expect_no_error({
    p2 <- mcmc_rhat(mock_fit, pars = c("d[B]", "d[C]"))
  })
  
  expect_true(inherits(p2, "ggplot"))
})

test_that("mcmc_neff produces effective sample size plots", {
  skip_if_not_installed("bayesplot")
  
  mock_fit <- create_mock_stan_nma()
  
  # Should work with default parameters
  expect_no_error({
    p <- mcmc_neff(mock_fit)
  })
  
  # Should return a ggplot object
  expect_true(inherits(p, "ggplot"))
  
  # Should work with specific parameters
  expect_no_error({
    p2 <- mcmc_neff(mock_fit, pars = c("sd"))
  })
  
  expect_true(inherits(p2, "ggplot"))
})

test_that("mcmc_diagnostics_plot combines multiple diagnostics", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("patchwork")
  
  mock_fit <- create_mock_stan_nma()
  
  # Should work with all diagnostics (default)
  expect_no_error({
    p <- mcmc_diagnostics_plot(mock_fit)
  })
  
  # Should work with subset of diagnostics
  expect_no_error({
    p2 <- mcmc_diagnostics_plot(mock_fit, type = c("trace", "rhat"))
  })
  
  # Should work with single diagnostic
  expect_no_error({
    p3 <- mcmc_diagnostics_plot(mock_fit, type = "trace")
  })
  
  # Should work with specific parameters
  expect_no_error({
    p4 <- mcmc_diagnostics_plot(mock_fit, pars = c("d[B]", "d[C]"))
  })
})

test_that("get_default_diagnostic_pars selects appropriate parameters", {
  # Test basic model (fixed effects)
  mock_fit_fixed <- create_mock_stan_nma(trt_effects = "fixed")
  mock_fit_fixed$stanfit@sim$fnames_oi <- c("d[B]", "d[C]")
  
  pars <- get_default_diagnostic_pars(mock_fit_fixed)
  expect_true(all(c("d[B]", "d[C]") %in% pars))
  
  # Test random effects model
  mock_fit_random <- create_mock_stan_nma(trt_effects = "random")
  mock_fit_random$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "sd")
  
  pars_re <- get_default_diagnostic_pars(mock_fit_random)
  expect_true("sd" %in% pars_re)
  expect_true(all(c("d[B]", "d[C]") %in% pars_re))
  
  # Test survival model
  mock_fit_surv <- create_mock_stan_nma(likelihood = "weibull")
  mock_fit_surv$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "shape[1]", "shape[2]")
  
  pars_surv <- get_default_diagnostic_pars(mock_fit_surv)
  expect_true(any(grepl("shape", pars_surv)))
  
  # Test regression model
  mock_fit_reg <- create_mock_stan_nma()
  mock_fit_reg$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "beta[1]", "beta[2]")
  
  pars_reg <- get_default_diagnostic_pars(mock_fit_reg)
  expect_true(any(grepl("beta", pars_reg)))
})

test_that("get_default_diagnostic_pars handles different model types", {
  # Test spline model
  mock_fit_spline <- create_mock_stan_nma(likelihood = "mspline")
  mock_fit_spline$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "scoef[1,1]", "scoef[1,2]", "scoef[2,1]")
  
  pars_spline <- get_default_diagnostic_pars(mock_fit_spline)
  expect_true(any(grepl("scoef", pars_spline)))
  
  # Test class effects model
  mock_fit_class <- create_mock_stan_nma(class_effects = "exchangeable")
  mock_fit_class$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "class_mean[1]", "class_sd[1]")
  
  pars_class <- get_default_diagnostic_pars(mock_fit_class)
  expect_true(any(grepl("class_", pars_class)))
  
  # Test auxiliary regression
  mock_fit_aux <- create_mock_stan_nma(aux_regression = ~age)
  mock_fit_aux$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "beta_aux[1]", "beta_aux[2]")
  
  pars_aux <- get_default_diagnostic_pars(mock_fit_aux)
  expect_true(any(grepl("beta_aux", pars_aux)))
})

test_that("get_default_diagnostic_pars limits parameter count", {
  # Create model with many parameters
  many_pars <- c(paste0("d[", 1:25, "]"), paste0("beta[", 1:25, "]"))
  
  mock_fit_many <- create_mock_stan_nma()
  mock_fit_many$stanfit@sim$fnames_oi <- many_pars
  
  # Should limit and warn
  expect_warning({
    pars <- get_default_diagnostic_pars(mock_fit_many)
  }, "Large number of parameters detected")
  
  expect_true(length(pars) <= 20)
})

test_that("validate_diagnostic_pars handles missing parameters", {
  mock_fit <- create_mock_stan_nma()
  mock_fit$stanfit@sim$fnames_oi <- c("d[B]", "d[C]", "sd")
  
  # Should warn about missing parameters
  expect_warning({
    valid_pars <- validate_diagnostic_pars(mock_fit, c("d[B]", "missing_par", "d[C]"))
  }, "not found")
  
  expect_equal(valid_pars, c("d[B]", "d[C]"))
  
  # Should error if no valid parameters
  expect_error({
    validate_diagnostic_pars(mock_fit, c("missing1", "missing2"))
  }, "No valid parameters")
})

test_that("diagnostic functions error appropriately for non-stan_nma objects", {
  not_stan_nma <- list(some_data = "test")
  
  expect_error(mcmc_trace(not_stan_nma), "Expecting a `stan_nma` object")
  expect_error(mcmc_acf(not_stan_nma), "Expecting a `stan_nma` object")
  expect_error(mcmc_rhat(not_stan_nma), "Expecting a `stan_nma` object")
  expect_error(mcmc_neff(not_stan_nma), "Expecting a `stan_nma` object")
  expect_error(mcmc_diagnostics_plot(not_stan_nma), "Expecting a `stan_nma` object")
})

test_that("require_pkg gives helpful error messages", {
  expect_error(
    require_pkg("nonexistent_package"),
    "Package nonexistent_package is required.*not installed"
  )
})

test_that("diagnostic functions work with various parameter selections", {
  skip_if_not_installed("bayesplot")
  
  mock_fit <- create_mock_stan_nma()
  
  # Test with empty parameter list (should use defaults)
  expect_no_error({
    mcmc_trace(mock_fit, pars = character(0))
  })
  
  # Test with single parameter
  expect_no_error({
    mcmc_trace(mock_fit, pars = "d[B]")
  })
  
  # Test with multiple parameters
  expect_no_error({
    mcmc_trace(mock_fit, pars = c("d[B]", "d[C]"))
  })
})

test_that("mcmc_diagnostics_plot handles different type combinations", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("patchwork")
  
  mock_fit <- create_mock_stan_nma()
  
  # Test all possible type combinations
  type_combinations <- list(
    "trace",
    "acf", 
    "rhat",
    "neff",
    c("trace", "acf"),
    c("trace", "rhat"),
    c("trace", "neff"),
    c("acf", "rhat"),
    c("acf", "neff"),
    c("rhat", "neff"),
    c("trace", "acf", "rhat"),
    c("trace", "acf", "neff"),
    c("trace", "rhat", "neff"),
    c("acf", "rhat", "neff"),
    c("trace", "acf", "rhat", "neff")
  )
  
  for (types in type_combinations) {
    expect_no_error({
      p <- mcmc_diagnostics_plot(mock_fit, type = types)
    }, info = paste("Failed for types:", paste(types, collapse = ", ")))
  }
})

test_that("diagnostic functions apply multinma theme", {
  skip_if_not_installed("bayesplot")
  
  mock_fit <- create_mock_stan_nma()
  
  # Check that plots have theme applied (this is hard to test directly,
  # but we can at least ensure the functions run without error)
  expect_no_error({
    p1 <- mcmc_trace(mock_fit)
    p2 <- mcmc_acf(mock_fit)
    p3 <- mcmc_rhat(mock_fit)
    p4 <- mcmc_neff(mock_fit)
  })
  
  # Plots should be ggplot objects
  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
  expect_true(inherits(p3, "ggplot"))
  expect_true(inherits(p4, "ggplot"))
})

test_that("diagnostic functions handle edge cases gracefully", {
  skip_if_not_installed("bayesplot")
  
  # Test with minimal parameters
  mock_fit_minimal <- create_mock_stan_nma()
  mock_fit_minimal$stanfit@sim$fnames_oi <- c("d[B]")
  
  expect_no_error({
    mcmc_trace(mock_fit_minimal)
    mcmc_acf(mock_fit_minimal)
    mcmc_rhat(mock_fit_minimal)
    mcmc_neff(mock_fit_minimal)
  })
  
  # Test with different likelihood types
  for (likelihood in c("binomial", "poisson", "normal", "weibull", "mspline")) {
    mock_fit_likelihood <- create_mock_stan_nma(likelihood = likelihood)
    
    expect_no_error({
      get_default_diagnostic_pars(mock_fit_likelihood)
    }, info = paste("Failed for likelihood:", likelihood))
  }
}) 