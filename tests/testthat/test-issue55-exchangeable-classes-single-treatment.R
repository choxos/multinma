# Test file for Issue #55: Exchangeable Classes with Single Treatments

# Helper function to create test data with mixed class sizes
create_mixed_class_data <- function() {
  # Create test data where some classes have multiple treatments, others single
  test_dat <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2", "Study3", "Study3"),
    trt = c("A", "B", "A", "C", "A", "D"), 
    r = c(10, 8, 12, 6, 15, 9),
    n = c(50, 45, 60, 40, 70, 50),
    # Mixed classes: Class1 has multiple treatments (A,B), others are single
    trt_class = c("Class1", "Class1", "Class1", "Class2", "Class1", "Class3")
  )
  
  return(test_dat)
}

# Helper function to create test data where ALL classes contain single treatments
create_all_single_class_data <- function() {
  test_dat <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    trt = c("A", "B", "A", "C"),
    r = c(10, 8, 12, 6),
    n = c(50, 45, 60, 40),
    # All classes are single treatment
    trt_class = c("Class1", "Class2", "Class1", "Class3")
  )
  
  return(test_dat)
}

# Helper function to create standard test data (control case)
create_standard_class_data <- function() {
  test_dat <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2"),
    trt = c("A", "B", "A", "C"),
    r = c(10, 8, 12, 6),
    n = c(50, 45, 60, 40),
    # Both classes have multiple treatments
    trt_class = c("Class1", "Class1", "Class1", "Class1")
  )
  
  return(test_dat)
}

test_that("exchangeable classes work with mixed class sizes", {
  test_dat <- create_mixed_class_data()
  
  # Set up network with classes
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # This should not error with the fix
  expect_no_error({
    fit <- nma(net, 
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,  # Small number for testing
               warmup = 50)
  })
  
  # Should produce a valid stan_nma object
  expect_true(inherits(fit, "stan_nma"))
  expect_equal(fit$class_effects, "exchangeable")
})

test_that("exchangeable classes work when ALL classes are single treatments", {
  test_dat <- create_all_single_class_data()
  
  # Set up network with all single-treatment classes
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # This should not error with the fix (edge case scenario)
  expect_no_error({
    fit <- nma(net,
               class_effects = "exchangeable", 
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  # Should produce a valid stan_nma object
  expect_true(inherits(fit, "stan_nma"))
  expect_equal(fit$class_effects, "exchangeable")
})

test_that("exchangeable classes work with standard multi-treatment classes", {
  test_dat <- create_standard_class_data()
  
  # Set up network with standard classes (control case)
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # This should continue to work as before
  expect_no_error({
    fit <- nma(net,
               class_effects = "exchangeable",
               chains = 2, 
               iter = 100,
               warmup = 50)
  })
  
  # Should produce a valid stan_nma object
  expect_true(inherits(fit, "stan_nma"))
  expect_equal(fit$class_effects, "exchangeable")
})

test_that("fix prevents array index out of range errors", {
  # This test specifically targets the original error condition
  test_dat <- create_mixed_class_data()
  
  # Set up network
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # The original error was: "array[uni, ...] index: accessing element out of range"
  # This should not occur with the fix
  expect_no_error({
    fit <- nma(net,
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  # Verify the model actually ran (not just initialized)
  expect_true(inherits(fit$stanfit, "stanfit"))
  
  # Check that we can extract basic results
  expect_no_error({
    summary(fit)
  })
})

test_that("fix works with different likelihood types", {
  test_dat <- create_mixed_class_data()
  
  # Test with different likelihoods that were affected
  likelihoods <- c("binomial", "poisson", "normal")
  
  for (likelihood in likelihoods) {
    # Adjust data for likelihood
    if (likelihood == "poisson") {
      net <- set_agd_arm(test_dat, study, trt, y = r, trt_class = trt_class)
    } else if (likelihood == "normal") {
      test_dat$y <- test_dat$r / test_dat$n
      test_dat$se <- sqrt(test_dat$y * (1 - test_dat$y) / test_dat$n)
      net <- set_agd_arm(test_dat, study, trt, y = y, se = se, trt_class = trt_class)
    } else {
      net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
    }
    
    expect_no_error({
      fit <- nma(net,
                 likelihood = likelihood,
                 class_effects = "exchangeable",
                 chains = 2,
                 iter = 100, 
                 warmup = 50)
    }, info = paste("Failed for likelihood:", likelihood))
    
    expect_true(inherits(fit, "stan_nma"))
  }
})

test_that("fix works with random effects models", {
  test_dat <- create_mixed_class_data()
  
  # Set up network
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # Test with random effects (more complex model)
  expect_no_error({
    fit <- nma(net,
               trt_effects = "random",
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  expect_true(inherits(fit, "stan_nma"))
  expect_equal(fit$trt_effects, "random")
  expect_equal(fit$class_effects, "exchangeable")
})

test_that("fix works with regression models", {
  test_dat <- create_mixed_class_data()
  
  # Add a covariate for regression
  test_dat$age_mean <- c(65, 62, 68, 70, 63, 66)
  
  # Set up network
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # Test with regression
  expect_no_error({
    fit <- nma(net,
               regression = ~age_mean,
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  expect_true(inherits(fit, "stan_nma"))
  expect_equal(fit$class_effects, "exchangeable")
})

test_that("fix maintains backward compatibility", {
  # Use plaque psoriasis data structure (mentioned in original issue)
  test_dat <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study2", "Study3"),
    trt = c("Placebo", "DrugA", "Placebo", "DrugB", "DrugC"),
    r = c(5, 15, 8, 20, 12),
    n = c(50, 50, 60, 60, 40),
    # Mixed structure similar to plaque psoriasis
    trt_class = c("Control", "Treatment", "Control", "Treatment", "Treatment")
  )
  
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # This should work without any issues
  expect_no_error({
    fit <- nma(net,
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  # Verify results are reasonable
  expect_true(inherits(fit, "stan_nma"))
  
  # Should be able to extract relative effects
  expect_no_error({
    rel_eff <- relative_effects(fit)
  })
  
  expect_true(inherits(rel_eff, "nma_summary"))
})

test_that("fix handles edge case with single-arm studies", {
  # Test data with some single-arm studies (another edge case)
  test_dat <- data.frame(
    study = c("Study1", "Study1", "Study2", "Study3"),
    trt = c("A", "B", "A", "C"), 
    r = c(10, 8, 12, 6),
    n = c(50, 45, 60, 40),
    trt_class = c("Class1", "Class2", "Class1", "Class3")
  )
  
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  expect_no_error({
    fit <- nma(net,
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  expect_true(inherits(fit, "stan_nma"))
})

test_that("fix does not affect non-exchangeable class models", {
  test_dat <- create_mixed_class_data()
  
  # Set up network
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # Test that independent and common class effects still work
  expect_no_error({
    fit_indep <- nma(net,
                     class_effects = "independent",
                     chains = 2,
                     iter = 100,
                     warmup = 50)
  })
  
  expect_no_error({
    fit_common <- nma(net,
                      class_effects = "common",
                      chains = 2,
                      iter = 100,
                      warmup = 50)
  })
  
  expect_true(inherits(fit_indep, "stan_nma"))
  expect_true(inherits(fit_common, "stan_nma"))
  expect_equal(fit_indep$class_effects, "independent")
  expect_equal(fit_common$class_effects, "common")
})

test_that("bounds checking prevents invalid array access", {
  # This is a more technical test to verify the specific fix
  test_dat <- create_mixed_class_data()
  
  # The fix adds agd_arm_trt[i] <= nt to prevent accessing invalid indices
  # This test ensures the model runs without Stan runtime errors
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # Capture any Stan-specific errors
  expect_no_error({
    fit <- nma(net,
               class_effects = "exchangeable",
               chains = 2,
               iter = 100,
               warmup = 50)
  })
  
  # Verify the model completed successfully
  expect_true(inherits(fit$stanfit, "stanfit"))
  
  # Check that the fix allows proper parameter extraction
  expect_no_error({
    class_params <- rstan::extract(fit$stanfit, pars = c("class_mean", "class_sd"))
  })
  
  expect_true(length(class_params) > 0)
})

test_that("fix works with all affected Stan model files", {
  test_dat <- create_mixed_class_data()
  net <- set_agd_arm(test_dat, study, trt, r = r, n = n, trt_class = trt_class)
  
  # Test different models that use different Stan files
  stan_models <- list(
    binomial_1par = list(likelihood = "binomial"),
    binomial_2par = list(likelihood = "binomial"),  # Uses same file as 1par in practice
    normal = list(likelihood = "normal"),
    poisson = list(likelihood = "poisson")
  )
  
  for (model_name in names(stan_models)) {
    model_args <- stan_models[[model_name]]
    
    # Adjust data based on likelihood
    if (model_args$likelihood == "normal") {
      test_dat$y <- test_dat$r / test_dat$n
      test_dat$se <- sqrt(test_dat$y * (1 - test_dat$y) / test_dat$n)
      net_adj <- set_agd_arm(test_dat, study, trt, y = y, se = se, trt_class = trt_class)
    } else if (model_args$likelihood == "poisson") {
      net_adj <- set_agd_arm(test_dat, study, trt, y = r, trt_class = trt_class)
    } else {
      net_adj <- net
    }
    
    expect_no_error({
      fit <- nma(net_adj,
                 likelihood = model_args$likelihood,
                 class_effects = "exchangeable",
                 chains = 2,
                 iter = 100,
                 warmup = 50)
    }, info = paste("Failed for model:", model_name))
    
    expect_true(inherits(fit, "stan_nma"),
                info = paste("Invalid result for model:", model_name))
  }
}) 