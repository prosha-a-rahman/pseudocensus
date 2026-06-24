#' Derive weights that recover the recipient covariates
#'
#' Identify vectors of weights for each recipient which, when applied to
#' weighted nearest neighbours estimates, recovers the recipient covariate.
#'
#' @inheritParams order_donor_indices
#' @param n_weights A positive integer defining the number of weights used in
#'   in the weighted nearest neighbours estimates, uniformly for all `recipients`.
#'   The number of weights cannot exceed the number of valid donors, i.e.,
#'   `n_weights <= length(donors)`. The default argument sets
#'   `n_weights = length(donors)`.
#'
#' @returns A `length(recipients)`-by-`n_weights` numeric matrix of weights
#' corresponding rowwise to `recipients` and columnwise to proximity ordering
#' `1:n_weights`.
#'
#' @export
#'
derive_local_debiasing_weights <- function(covariates,
                                           recipients,
                                           donors,
                                           n_weights = NULL) {
  if (is.null(n_weights)) n_weights <- length(donors)

  ordered_donor_covariates <- order_donor_covariates(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    k = n_weights
  )

  # Transpose matrices so that each column is a point in the covariate space
  # to appropriately define linear system
  ordered_donor_covariates <- lapply(ordered_donor_covariates, t)

  recipient_covariates <- covariates[recipients, , drop = FALSE]

  solve_linear_system_wrapper <- function(idx) {
    solve_linear_system(
      A = ordered_donor_covariates[[idx]],
      b = recipient_covariates[idx, ]
    )
  }

  K <- ncol(ordered_donor_covariates[[1]])

  recipient_idx <- seq_along(ordered_donor_covariates)

  vapply(recipient_idx, solve_linear_system_wrapper, numeric(K)) |> t()
}




#' Compute global bias minimising weights
#'
#' Compute a single vector of weights that, when applied uniformly for all
#' `recipients`, minimises the bias of the pseudocensus OLS estimator.
#'
#' @details
#' It is implied that the `covariates` is a usual design matrix prepended by
#' a vector of ones, i.e., `covariates[1, ] = rep(1, nrow(covariates))`.
#'
#'
#' @inheritParams derive_local_debiasing_weights
#'
#' @returns A numeric vector of length `n_weights`.
#'
#' @export
derive_global_debiasing_weights <- function(covariates,
                                            recipients,
                                            donors,
                                            n_weights) {
  ordered_donor_covariates <- order_donor_covariates(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    k = n_weights
  )

  ordered_design_matrices <- transpose_list(L = ordered_donor_covariates)

  recipient_design_matrix <- covariates[recipients, , drop = FALSE]

  block_matrix <- do.call(cbind, ordered_design_matrices)

  A <- crossprod(recipient_design_matrix, block_matrix)

  n_parameters <- ncol(covariates)

  A <- matrix(A, nrow = n_parameters^2, ncol = n_weights)

  b <- crossprod(recipient_design_matrix) |> as.numeric()

  solve_linear_system(A, b)
}
