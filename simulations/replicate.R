# Single-replicate logic: draw one population, apply every estimator under
# study, and return the coefficient estimates as a matrix.

ESTIMATOR_NAMES <- c(
  "census",          # oracle OLS on the fully observed population
  "sample_ols",      # complete-case OLS on donors only
  "pseudo_knn",      # pseudocensus OLS, uniform k-NN imputation weights
  "pseudo_local",    # pseudocensus OLS, local debiasing weights
  "pseudo_global",   # pseudocensus OLS, global debiasing weights
  "debiased_knn",    # bias-corrected pseudocensus OLS, uniform weights
  "debiased_local",  # bias-corrected pseudocensus OLS, local weights
  "debiased_global"  # bias-corrected pseudocensus OLS, global weights
)

#' Run one Monte Carlo replicate for a design cell.
#'
#' The seed is a deterministic function of (cell, rep_id) alone, so results
#' are reproducible regardless of how replicates are scheduled across workers.
#'
#' The design matrix cbind(1, x) is passed to every routine, including the
#' weight derivations: reproducing the intercept column imposes the
#' sum-to-one constraint on each recipient's weights.
#'
#' The response vector handed to the pseudocensus estimators has recipient
#' entries set to NA, guaranteeing no leakage of unobserved responses.
#'
#' @returns A length(ESTIMATOR_NAMES)-by-(p + 1) matrix of estimates.
run_replicate <- function(cfg, rep_id, base_seed = 894651L) {
  set.seed(base_seed + 1000003L * cfg$cell + rep_id)

  pop <- generate_population(
    n = cfg$n,
    p = cfg$p,
    mean_spec = cfg$mean_spec,
    error_spec = cfg$error_spec,
    missing_spec = cfg$missing_spec,
    covariate_design = cfg$covariate_design,
    rho = cfg$rho,
    response_rate = cfg$response_rate
  )

  y_observed <- replace(pop$y, !pop$observed, NA)
  donors <- identify_donors(y_observed)
  recipients <- identify_recipients(y_observed)

  X <- cbind(1, pop$x)
  k <- min(cfg$k, length(donors))
  r <- length(recipients)

  w_uniform <- matrix(1 / k, nrow = r, ncol = k)

  w_local <- derive_local_debiasing_weights(
    x = X,
    recipients = recipients,
    donors = donors,
    n_weights = k
  )

  w_global <- derive_global_debiasing_weights(
    x = X,
    recipients = recipients,
    donors = donors,
    n_weights = k
  ) |>
    matrix(nrow = r, ncol = k, byrow = TRUE)

  pseudo <- function(weights) {
    fit_pseudocensus_ols(X, y_observed, donors, recipients, weights)
  }

  debiased <- function(weights) {
    fit_debiased_pseudocensus_ols(X, y_observed, donors, recipients, weights)
  }

  estimates <- rbind(
    census = stats::.lm.fit(X, pop$y)$coefficients,
    sample_ols = fit_sample_ols(X, pop$y, donors),
    pseudo_knn = pseudo(w_uniform),
    pseudo_local = pseudo(w_local),
    pseudo_global = pseudo(w_global),
    debiased_knn = debiased(w_uniform),
    debiased_local = debiased(w_local),
    debiased_global = debiased(w_global)
  )

  colnames(estimates) <- c("(Intercept)", paste0("x", seq_len(cfg$p)))
  estimates
}
