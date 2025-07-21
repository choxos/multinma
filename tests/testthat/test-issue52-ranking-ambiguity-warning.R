# Test file for Issue #52: TTE Parametric Distribution Ranking Warnings

# Mock stan_nma object creation helper
create_mock_stan_nma <- function(likelihood = "weibull", aux_by = ".study", aux_regression = NULL, consistency = "consistency") {
  # Create minimal mock stan_nma object
  mock_network <- list(
    treatments = factor(c("A", "B", "C"))
  )
  
  # Create mock aux_regression with proper terms structure if provided
  if (!is.null(aux_regression)) {
    # Create a minimal formula object with proper terms
    aux_formula <- as.formula(aux_regression)
    aux_terms <- terms(aux_formula)
    attr(aux_terms, "factor") <- matrix(0, nrow = 0, ncol = length(attr(aux_terms, "term.labels")))
    colnames(attr(aux_terms, "factor")) <- attr(aux_terms, "term.labels")
  } else {
    aux_terms <- NULL
  }
  
  structure(
    list(
      likelihood = likelihood,
      aux_by = aux_by,
      aux_regression = aux_terms,
      consistency = consistency,
      network = mock_network
    ),
    class = "stan_nma"
  )
}

# Mock relative_effects function for testing
relative_effects <- function(x, newdata = NULL, study = NULL, all_contrasts = FALSE, summary = FALSE) {
  # Return minimal structure for testing
  ntrt <- nlevels(x$network$treatments)
  n_iter <- 100
  n_chains <- 4
  
  # Create mock MCMC array
  sim_array <- array(rnorm(n_iter * n_chains * (ntrt - 1)), 
                     dim = c(n_iter, n_chains, ntrt - 1),
                     dimnames = list(
                       iterations = 1:n_iter,
                       chains = 1:n_chains,
                       parameters = paste0("d[", levels(x$network$treatments)[-1], "]")
                     ))
  
  return(list(sim = sim_array, studies = NULL))
}

test_that("posterior_ranks warns for TTE parametric models with aux_by including .trt", {
  # Test Weibull model with .trt in aux_by
  weibull_model <- create_mock_stan_nma(
    likelihood = "weibull", 
    aux_by = c(".study", ".trt")
  )
  
  expect_warning(
    posterior_ranks(weibull_model),
    "Treatment rankings may be ambiguous.*time-to-event parametric models.*aux_by"
  )
  
  # Test Gompertz model with .trt in aux_by
  gompertz_model <- create_mock_stan_nma(
    likelihood = "gompertz",
    aux_by = c(".study", ".trt")
  )
  
  expect_warning(
    posterior_ranks(gompertz_model),
    "Treatment rankings may be ambiguous"
  )
  
  # Test multiple TTE parametric distributions
  tte_likelihoods <- c("weibull-aft", "lognormal", "loglogistic", "gamma", "gengamma")
  
  for (likelihood in tte_likelihoods) {
    model <- create_mock_stan_nma(
      likelihood = likelihood,
      aux_by = c(".study", ".trt")
    )
    
    expect_warning(
      posterior_ranks(model),
      "Treatment rankings may be ambiguous",
      info = paste("Failed for likelihood:", likelihood)
    )
  }
})

test_that("posterior_ranks warns for TTE parametric models with aux_regression including .trt", {
  # Test with simple .trt term
  weibull_model_reg <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_regression = "~.trt"
  )
  
  expect_warning(
    posterior_ranks(weibull_model_reg),
    "Treatment rankings may be ambiguous.*aux_regression including .trt"
  )
  
  # Test with .trt interaction
  weibull_model_interaction <- create_mock_stan_nma(
    likelihood = "weibull", 
    aux_regression = "~.trt + covariate + .trt:covariate"
  )
  
  expect_warning(
    posterior_ranks(weibull_model_interaction),
    "Treatment rankings may be ambiguous"
  )
  
  # Test Gamma model with .trt regression
  gamma_model_reg <- create_mock_stan_nma(
    likelihood = "gamma",
    aux_regression = "~.trt + age"
  )
  
  expect_warning(
    posterior_ranks(gamma_model_reg),
    "Treatment rankings may be ambiguous"
  )
})

test_that("posterior_ranks does NOT warn for TTE models without auxiliary treatment effects", {
  # Test with default aux_by (no .trt)
  weibull_model_no_trt <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_by = ".study"  # Default - no .trt
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(weibull_model_no_trt))
  )
  
  # Test with aux_regression not including .trt
  weibull_model_reg_no_trt <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_regression = "~age + gender"  # No .trt
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(weibull_model_reg_no_trt))
  )
})

test_that("posterior_ranks does NOT warn for non-TTE parametric models", {
  # Test binomial model with .trt in aux_by (should not warn)
  binomial_model <- create_mock_stan_nma(
    likelihood = "binomial",
    aux_by = c(".study", ".trt")
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(binomial_model))
  )
  
  # Test normal model with .trt regression (should not warn) 
  normal_model <- create_mock_stan_nma(
    likelihood = "normal",
    aux_regression = "~.trt + covariate"
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(normal_model))
  )
  
  # Test Poisson model
  poisson_model <- create_mock_stan_nma(
    likelihood = "poisson",
    aux_by = c(".study", ".trt")
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(poisson_model))
  )
})

test_that("posterior_ranks does NOT warn for exponential models (no auxiliary parameters)", {
  # Exponential models don't have auxiliary parameters, so should not warn
  exp_model <- create_mock_stan_nma(
    likelihood = "exponential",
    aux_by = c(".study", ".trt")  # This should be ignored for exponential
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(exp_model))
  )
  
  # Exponential AFT model
  exp_aft_model <- create_mock_stan_nma(
    likelihood = "exponential-aft",
    aux_regression = "~.trt"
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(exp_aft_model))
  )
})

test_that("posterior_ranks does NOT warn for spline-based survival models", {
  # M-spline models have different ambiguity context
  mspline_model <- create_mock_stan_nma(
    likelihood = "mspline",
    aux_by = c(".study", ".trt")
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(mspline_model))
  )
  
  # Piecewise exponential model
  pexp_model <- create_mock_stan_nma(
    likelihood = "pexp", 
    aux_regression = "~.trt"
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(pexp_model))
  )
})

test_that("posterior_rank_probs triggers warnings via posterior_ranks", {
  # Since posterior_rank_probs calls posterior_ranks internally, it should also trigger warnings
  weibull_model <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_by = c(".study", ".trt")
  )
  
  expect_warning(
    posterior_rank_probs(weibull_model),
    "Treatment rankings may be ambiguous"
  )
})

test_that("warning message contains all required information", {
  weibull_model <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_by = c(".study", ".trt")
  )
  
  # Capture the warning message
  warning_msg <- ""
  withCallingHandlers(
    posterior_ranks(weibull_model),
    warning = function(w) {
      warning_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  
  # Check that warning contains key information
  expect_true(grepl("time-to-event parametric models", warning_msg))
  expect_true(grepl("auxiliary treatment effects", warning_msg))
  expect_true(grepl("main treatment effects \\(d\\) only", warning_msg))
  expect_true(grepl("d_aux", warning_msg))
  expect_true(grepl("Consider examining both sets of effects", warning_msg))
})

test_that("ranking calculations proceed normally after warning", {
  weibull_model <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_by = c(".study", ".trt")
  )
  
  # Should warn but still return valid results
  result <- expect_warning(
    posterior_ranks(weibull_model),
    "Treatment rankings may be ambiguous"
  )
  
  # Check that results are still generated
  expect_true(inherits(result, "nma_ranks"))
  expect_true(inherits(result, "nma_summary"))
  expect_true("summary" %in% names(result))
  expect_true("sims" %in% names(result))
})

test_that("warning detection handles edge cases correctly", {
  # Test NULL aux_by
  model_null_aux_by <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_by = NULL
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(model_null_aux_by))
  )
  
  # Test empty aux_by
  model_empty_aux_by <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_by = character(0)
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(model_empty_aux_by))
  )
  
  # Test NULL aux_regression
  model_null_aux_reg <- create_mock_stan_nma(
    likelihood = "weibull",
    aux_regression = NULL
  )
  
  expect_silent(
    suppressMessages(posterior_ranks(model_null_aux_reg))
  )
}) 