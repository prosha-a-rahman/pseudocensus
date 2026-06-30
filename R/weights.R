#' Derive nearest neighbours weights that reproduce the recipient covariates
#'
#' Identify vectors of weights for each recipient which, when applied to
#' weighted nearest neighbours, reproduces the original recipient covariate.
#'
#' @inheritParams order_donor
#' @param n_weights An integer in `1:length(donors)` defining the number of
#'   donors to be weighted per recipient. Defaults to `length(donors)` which
#'   utilises the entire donor pool.
#'
#' @returns A `length(recipients)`-by-`n_weights` numeric matrix of weights
#' corresponding rowwise to `recipients` and columnwise to proximity ordering
#' `1:n_weights`.
#'
#' @export
derive_local_debiasing_weights <- function(x,
                                           recipients,
                                           donors,
                                           n_weights = length(donors)) {
  ordered_donor_x <- order_donor_covariates(
    x = x,
    recipients = recipients,
    donors = donors,
    k = n_weights
  )

  # Transpose matrices so that each column is a point in the covariate space
  # to appropriately define linear system
  ordered_donor_x <- lapply(ordered_donor_x, t)

  x_recipient <- x[recipients, , drop = FALSE]

  solve_linear_system_wrapper <- function(idx) {
    solve_linear_system(
      A = ordered_donor_x[[idx]],
      b = x_recipient[idx, ]
    )
  }

  K <- ncol(ordered_donor_x[[1]])

  recipient_idx <- seq_along(ordered_donor_x)

  vapply(recipient_idx, solve_linear_system_wrapper, numeric(K)) |> t()
}




#' Derive a global vector of weights that minimises aggregate covariate
#' imputation error
#'
#' Compute a single vector of weights which, when applied uniformly for all
#' `recipients`, minimises the aggregate Euclidean distance between imputed and
#' actual covariates.
#'
#'
#'
#' @inheritParams derive_local_debiasing_weights
#'
#' @returns A numeric vector of length `n_weights`.
#'
#' @export
derive_global_debiasing_weights <- function(x,
                                            recipients,
                                            donors,
                                            n_weights = length(donors)) {
  ordered_donor_x <- order_donor_covariates(
    x = x,
    recipients = recipients,
    donors = donors,
    k = n_weights
  )

  ordered_donor_x <- transpose_list(L = ordered_donor_x)

  x_recipients <- x[recipients, , drop = FALSE]

  block_matrix <- do.call(cbind, ordered_donor_x)

  A <- crossprod(x_recipients, block_matrix)

  p <- ncol(x)

  A <- matrix(A, nrow = p^2, ncol = n_weights)

  b <- crossprod(x_recipients) |> as.numeric()

  solve_linear_system(A, b)
}
