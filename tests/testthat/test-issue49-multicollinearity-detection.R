# Test file for Issue #49: Enhanced Multicollinearity Detection in add_integration()

# Helper function to suppress all warnings and messages for testing
quiet <- function(expr) {
  suppressWarnings(suppressMessages(expr))
}

test_that("detect_multicollinearity detects high pairwise correlations", {
  # High correlation matrix
  high_cor_matrix <- matrix(c(1, 0.97, 0.97, 1), nrow = 2)
  var_names <- c("age", "age_squared")
  
  expect_warning(
    detect_multicollinearity(high_cor_matrix, var_names),
    "High correlations detected.*age & age_squared.*r = 0.97"
  )
  
  # Multiple high correlations
  multi_high_cor <- matrix(c(1, 0.96, 0.98, 
                             0.96, 1, 0.97,
                             0.98, 0.97, 1), nrow = 3)
  var_names_multi <- c("x1", "x2", "x3")
  
  expect_warning(
    detect_multicollinearity(multi_high_cor, var_names_multi),
    "High correlations detected"
  )
})

test_that("detect_multicollinearity detects near-singular matrices", {
  # Near-singular matrix (very small eigenvalue)
  near_singular <- matrix(c(1, 0.999999, 0.999999, 1), nrow = 2)
  var_names <- c("var1", "var2")
  
  expect_warning(
    detect_multicollinearity(near_singular, var_names),
    "near-singular.*smallest eigenvalue"
  )
})

test_that("detect_multicollinearity detects high condition number", {
  # Matrix with high condition number but not near-singular
  high_condition <- matrix(c(1, 0.9, 0.85,
                             0.9, 1, 0.8,
                             0.85, 0.8, 1), nrow = 3)
  var_names <- c("x1", "x2", "x3")
  
  # Force high condition number by making it more ill-conditioned
  eigenvals <- eigen(high_condition)$values
  # Modify to create high condition number
  high_condition[1,2] <- high_condition[2,1] <- 0.95
  high_condition[1,3] <- high_condition[3,1] <- 0.94
  high_condition[2,3] <- high_condition[3,2] <- 0.93
  
  expect_warning(
    detect_multicollinearity(high_condition, var_names, threshold_condition = 10),
    "high condition number"
  )
})

test_that("detect_multicollinearity provides positive feedback for well-conditioned matrices", {
  # Well-conditioned matrix
  well_conditioned <- matrix(c(1, 0.3, 0.2,
                               0.3, 1, 0.1,
                               0.2, 0.1, 1), nrow = 3)
  var_names <- c("x1", "x2", "x3")
  
  expect_message(
    detect_multicollinearity(well_conditioned, var_names),
    "well-conditioned"
  )
})

test_that("detect_multicollinearity handles single variable case", {
  single_var_matrix <- matrix(1, nrow = 1)
  var_names <- c("single_var")
  
  expect_message(
    detect_multicollinearity(single_var_matrix, var_names),
    "Single covariate.*no multicollinearity assessment needed"
  )
})

test_that("add_integration triggers multicollinearity detection for user-provided correlations", {
  # Create test data
  test_data <- data.frame(study = 1, x1_mean = 2, x1_sd = 0.5, x2_mean = 3, x2_sd = 0.8)
  
  # High correlation matrix
  high_cor <- matrix(c(1, 0.97, 0.97, 1), nrow = 2)
  
  expect_warning(
    add_integration(test_data,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qnorm, mean = x2_mean, sd = x2_sd),
                    cor = high_cor),
    "High correlations detected"
  )
})

test_that("add_integration triggers multicollinearity detection for IPD-derived correlations", {
  skip_if_not_installed("dplyr")
  
  # Create test network with IPD that has high correlations
  ipd_data <- data.frame(
    .study = rep("Study1", 100),
    .trt = rep("A", 100),
    .y = rnorm(100),
    x1 = rnorm(100),
    x2 = NA  # Will be set to create high correlation
  )
  
  # Create high correlation between x1 and x2
  ipd_data$x2 <- ipd_data$x1 + rnorm(100, 0, 0.1)  # Very high correlation
  
  # Create a second study to avoid warnings about single study
  ipd_data2 <- ipd_data
  ipd_data2$.study <- "Study2"
  ipd_data2$.trt <- "B"
  ipd_data_combined <- rbind(ipd_data, ipd_data2)
  
  agd_data <- data.frame(
    .study = "Study3",
    .trt = "C",
    .y = 0,
    .se = 1,
    x1_mean = 2,
    x1_sd = 0.5,
    x2_mean = 3,
    x2_sd = 0.8
  )
  
  # Create network
  network <- set_ipd(ipd_data_combined, study = .study, trt = .trt, y = .y)
  network <- set_agd_arm(network, agd_data, study = .study, trt = .trt, y = .y, se = .se)
  
  expect_warning(
    add_integration(network,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qnorm, mean = x2_mean, sd = x2_sd)),
    "High correlations detected"
  )
})

test_that("add_integration triggers multicollinearity detection for adjusted correlations", {
  # Create test data
  test_data <- data.frame(study = 1, x1_mean = 2, x1_sd = 0.5, x2_mean = 3, x2_sd = 0.8)
  
  # Correlation matrix that will become problematic after adjustment
  problematic_cor <- matrix(c(1, 0.9, 0.9, 1), nrow = 2)
  
  # This should trigger detection after correlation adjustment
  expect_warning(
    add_integration(test_data,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qbern, prob = 0.5),  # Binary variable to trigger adjustment
                    cor = problematic_cor,
                    cor_adjust = "spearman"),
    "multicollinearity|condition number",
    fixed = FALSE
  )
})

test_that("multicollinearity detection works with custom thresholds", {
  # Matrix with moderate correlation
  moderate_cor <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
  var_names <- c("x1", "x2")
  
  # Should not warn with default threshold (0.95)
  expect_silent(
    quiet(detect_multicollinearity(moderate_cor, var_names))
  )
  
  # Should warn with lower threshold
  expect_warning(
    detect_multicollinearity(moderate_cor, var_names, threshold_high = 0.7),
    "High correlations detected"
  )
  
  # Test custom condition number threshold
  expect_warning(
    detect_multicollinearity(moderate_cor, var_names, threshold_condition = 2),
    "high condition number"
  )
})

test_that("multicollinearity detection handles edge cases gracefully", {
  # Perfect correlation (singular matrix)
  perfect_cor <- matrix(c(1, 1, 1, 1), nrow = 2)
  var_names <- c("x1", "x2")
  
  expect_warning(
    detect_multicollinearity(perfect_cor, var_names),
    "near-singular"
  )
  
  # Identity matrix (no correlations)
  identity_cor <- diag(3)
  var_names_identity <- c("x1", "x2", "x3")
  
  expect_message(
    detect_multicollinearity(identity_cor, var_names_identity),
    "well-conditioned"
  )
})

test_that("multicollinearity detection integrates properly with existing workflow", {
  # Test that detection doesn't break normal functionality
  test_data <- data.frame(study = 1, x1_mean = 2, x1_sd = 0.5, x2_mean = 3, x2_sd = 0.8)
  
  # Well-conditioned matrix should work normally
  good_cor <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  
  result <- expect_message(
    add_integration(test_data,
                    x1 = distr(qnorm, mean = x1_mean, sd = x1_sd),
                    x2 = distr(qnorm, mean = x2_mean, sd = x2_sd),
                    cor = good_cor),
    "well-conditioned"
  )
  
  # Check that integration points were actually generated
  expect_true(inherits(result, "integration_tbl"))
  expect_true(all(c(".int_x1", ".int_x2") %in% colnames(result)))
})

test_that("plot_integration_error placeholder function works", {
  # Test the placeholder function
  expect_message(
    plot_integration_error(NULL),
    "Integration error plotting functionality will be available in future versions"
  )
  
  expect_null(plot_integration_error(NULL))
}) 