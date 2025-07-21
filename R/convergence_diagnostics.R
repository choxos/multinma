#' MCMC Convergence Diagnostics
#'
#' Functions for assessing MCMC convergence and sampling efficiency from fitted
#' multinma models. These functions provide easy access to standard diagnostic
#' plots and statistics to help evaluate whether MCMC chains have converged.
#'
#' @param x A [stan_nma] object created by [nma()]
#' @param pars Character vector of parameter names to include in diagnostics.
#'   If `NULL` (default), a sensible subset of parameters is automatically
#'   selected based on the model type.
#' @param ... Additional arguments passed to bayesplot functions
#'
#' @return For plotting functions, a ggplot object. For diagnostic functions,
#'   the return depends on the specific diagnostic.
#'
#' @details These functions integrate with the bayesplot package to provide
#'   high-quality visualizations of MCMC diagnostics. The parameter selection
#'   when `pars = NULL` includes:
#'   
#'   - Treatment effects (`d[...]`)
#'   - Regression coefficients (`beta[...]`) if regression is used
#'   - Heterogeneity parameters (`sd`) if random effects are used
#'   - Auxiliary parameters for survival models
#'   - Class effects parameters if applicable
#'
#' @name mcmc_diagnostics
#' @examples
#' \donttest{
#' # Smoking cessation example
#' library(multinma)
#' smk_net <- set_agd_arm(smoking, study, trt, r = r, n = n)
#' smk_fit <- nma(smk_net, trt_effects = "random")
#' 
#' # Trace plots
#' mcmc_trace(smk_fit)
#' 
#' # Autocorrelation plots  
#' mcmc_acf(smk_fit)
#' 
#' # R-hat diagnostics
#' mcmc_rhat(smk_fit)
#' 
#' # Effective sample size
#' mcmc_neff(smk_fit)
#' 
#' # Combined diagnostics plot
#' mcmc_diagnostics_plot(smk_fit)
#' }
NULL

#' @rdname mcmc_diagnostics
#' @export
mcmc_trace <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("Expecting a `stan_nma` object.")
  }
  
  # Get default parameters if not specified
  if (is.null(pars)) {
    pars <- get_default_diagnostic_pars(x)
  }
  
  # Validate parameters
  pars <- validate_diagnostic_pars(x, pars)
  
  # Extract MCMC array
  mcmc_array <- as.array(x, pars = pars)
  
  # Create trace plot using bayesplot
  require_pkg("bayesplot")
  p <- bayesplot::mcmc_trace(mcmc_array, ...)
  
  # Apply multinma theme
  p <- p + theme_multinma()
  
  return(p)
}

#' @rdname mcmc_diagnostics
#' @export  
mcmc_acf <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("Expecting a `stan_nma` object.")
  }
  
  # Get default parameters if not specified
  if (is.null(pars)) {
    pars <- get_default_diagnostic_pars(x)
  }
  
  # Validate parameters
  pars <- validate_diagnostic_pars(x, pars)
  
  # Extract MCMC array
  mcmc_array <- as.array(x, pars = pars)
  
  # Create autocorrelation plot using bayesplot
  require_pkg("bayesplot")
  p <- bayesplot::mcmc_acf(mcmc_array, ...)
  
  # Apply multinma theme
  p <- p + theme_multinma()
  
  return(p)
}

#' @rdname mcmc_diagnostics
#' @export
mcmc_rhat <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("Expecting a `stan_nma` object.")
  }
  
  # Get default parameters if not specified
  if (is.null(pars)) {
    pars <- get_default_diagnostic_pars(x)
  }
  
  # Validate parameters
  pars <- validate_diagnostic_pars(x, pars)
  
  # Get R-hat values
  rhat_vals <- rstan::summary(x$stanfit, pars = pars)$summary[, "Rhat"]
  
  # Create R-hat plot using bayesplot
  require_pkg("bayesplot")
  p <- bayesplot::mcmc_rhat(rhat_vals, ...)
  
  # Apply multinma theme
  p <- p + theme_multinma()
  
  return(p)
}

#' @rdname mcmc_diagnostics
#' @export
mcmc_neff <- function(x, pars = NULL, ...) {
  if (!inherits(x, "stan_nma")) {
    abort("Expecting a `stan_nma` object.")
  }
  
  # Get default parameters if not specified
  if (is.null(pars)) {
    pars <- get_default_diagnostic_pars(x)
  }
  
  # Validate parameters
  pars <- validate_diagnostic_pars(x, pars)
  
  # Get effective sample size values
  stan_summary <- rstan::summary(x$stanfit, pars = pars)$summary
  neff_vals <- stan_summary[, "n_eff"]
  
  # Total number of iterations (across all chains)
  total_iter <- x$stanfit@sim$iter * x$stanfit@sim$chains
  neff_ratio <- neff_vals / total_iter
  
  # Create effective sample size plot using bayesplot
  require_pkg("bayesplot")
  p <- bayesplot::mcmc_neff(neff_ratio, ...)
  
  # Apply multinma theme
  p <- p + theme_multinma()
  
  return(p)
}

#' @rdname mcmc_diagnostics
#' @param type Character vector specifying which diagnostics to include.
#'   Options are "trace", "acf", "rhat", and "neff". Default includes all.
#' @export
mcmc_diagnostics_plot <- function(x, pars = NULL, type = c("trace", "acf", "rhat", "neff"), ...) {
  if (!inherits(x, "stan_nma")) {
    abort("Expecting a `stan_nma` object.")
  }
  
  require_pkg("patchwork")
  
  # Validate type argument
  type <- match.arg(type, several.ok = TRUE)
  
  # Get default parameters if not specified
  if (is.null(pars)) {
    pars <- get_default_diagnostic_pars(x)
  }
  
  # Create individual plots
  plots <- list()
  
  if ("trace" %in% type) {
    plots$trace <- mcmc_trace(x, pars = pars, ...)
  }
  
  if ("acf" %in% type) {
    plots$acf <- mcmc_acf(x, pars = pars, ...)
  }
  
  if ("rhat" %in% type) {
    plots$rhat <- mcmc_rhat(x, pars = pars, ...)
  }
  
  if ("neff" %in% type) {
    plots$neff <- mcmc_neff(x, pars = pars, ...)
  }
  
  # Combine plots using patchwork
  if (length(plots) == 1) {
    combined_plot <- plots[[1]]
  } else if (length(plots) == 2) {
    combined_plot <- plots[[1]] / plots[[2]]
  } else if (length(plots) == 3) {
    combined_plot <- plots[[1]] / (plots[[2]] | plots[[3]])
  } else {
    combined_plot <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
  }
  
  return(combined_plot)
}

#' Get default diagnostic parameters
#' 
#' Internal function to select sensible default parameters for MCMC diagnostics
#' based on the model structure.
#' 
#' @param x A stan_nma object
#' @return Character vector of parameter names
#' @noRd
get_default_diagnostic_pars <- function(x) {
  # Get all parameter names
  all_pars <- x$stanfit@sim$fnames_oi
  
  # Priority parameters to include
  pars <- character(0)
  
  # Always include treatment effects
  d_pars <- all_pars[grepl("^d\\[", all_pars)]
  pars <- c(pars, d_pars)
  
  # Include regression coefficients if present
  beta_pars <- all_pars[grepl("^beta\\[", all_pars)]
  pars <- c(pars, beta_pars)
  
  # Include heterogeneity parameter if random effects
  if (x$trt_effects == "random") {
    sd_pars <- all_pars[grepl("^sd$", all_pars)]
    pars <- c(pars, sd_pars)
  }
  
  # Include auxiliary parameters for survival models
  if (x$likelihood %in% c("weibull", "gompertz", "weibull-aft", "lognormal", 
                          "loglogistic", "gamma", "gengamma")) {
    aux_pars <- all_pars[grepl("^(shape|sigma|k|sdlog)\\[", all_pars)]
    pars <- c(pars, aux_pars)
  }
  
  # Include spline coefficients for flexible survival models
  if (x$likelihood %in% c("mspline", "pexp")) {
    scoef_pars <- all_pars[grepl("^scoef\\[", all_pars)]
    # Limit to avoid too many parameters
    pars <- c(pars, head(scoef_pars, 10))
  }
  
  # Include class effects if present
  if (!is.null(x$class_effects) && x$class_effects == "exchangeable") {
    class_pars <- all_pars[grepl("^(class_mean|class_sd)\\[", all_pars)]
    pars <- c(pars, class_pars)
  }
  
  # Include auxiliary regression parameters
  if (!is.null(x$aux_regression)) {
    aux_reg_pars <- all_pars[grepl("^beta_aux\\[", all_pars)]
    pars <- c(pars, aux_reg_pars)
  }
  
  # Remove any duplicates and limit total number
  pars <- unique(pars)
  
  # Limit to reasonable number for visualization
  if (length(pars) > 20) {
    warn("Large number of parameters detected. Limiting to first 20 for diagnostics. Use `pars` argument to specify custom selection.")
    pars <- pars[1:20]
  }
  
  return(pars)
}

#' Validate diagnostic parameters
#' 
#' Internal function to validate that requested parameters exist in the model.
#' 
#' @param x A stan_nma object
#' @param pars Character vector of parameter names
#' @return Validated character vector of parameter names
#' @noRd
validate_diagnostic_pars <- function(x, pars) {
  # Get all available parameter names
  all_pars <- x$stanfit@sim$fnames_oi
  
  # Check which requested parameters are available
  missing_pars <- setdiff(pars, all_pars)
  
  if (length(missing_pars) > 0) {
    warn(paste("The following parameters were not found and will be ignored:",
               paste(missing_pars, collapse = ", ")))
  }
  
  # Return only available parameters
  available_pars <- intersect(pars, all_pars)
  
  if (length(available_pars) == 0) {
    abort("No valid parameters specified for diagnostics.")
  }
  
  return(available_pars)
}

#' Require package with informative error
#' 
#' Internal function to check for required packages and give helpful error messages.
#' 
#' @param pkg Package name
#' @noRd
require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    abort(paste("Package", pkg, "is required for MCMC diagnostics but is not installed.",
                "Install it with: install.packages(\"", pkg, "\")", sep = ""))
  }
} 