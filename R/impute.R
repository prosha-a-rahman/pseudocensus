#' Impute covariates or responses using weighted nearest neighbours
#'
#' Respectively construct synthetic covariates or responses for a set of
#' recipients from weighted combinations of their nearest donor covariates or
#' responses.
#'
#' @inheritParams order_donor
#' @param weights A numeric matrix where `nrow(weights) = length(recipients)`
#'   and `ncol(weights) <= length(donors)`. Each row in `weights` supplies the
#'   imputation weights used in weighted nearest neighbours computation.
#'   `weights` corresponds rowwise to `recipients` and columnwise to the rank of
#'   the proximally sorted donors.
#' @param y A length `nrow(x)` numeric vector of responses indexed consistently
#'   with the rows of `x`. Only used by `impute_responses()`.
#'
#' @returns
#' Writing `n_recipients = length(recipients)`, `p = ncol(x)`, and
#' `k = ncol(weights)`:
#'
#' * `impute_covariates()`: an `n_recipients`-by-`p` numeric matrix, where each
#'   row is the weighted combination of the `k` nearest donor covariates
#'   corresponding to `recipients`.
#'
#' * `impute_responses()`: a length `n_recipients` numeric vector where each
#'   element is the weighted combination of responses for the `k` nearest donors
#'   corresponding to `recipients`.
#'
#' @seealso [order_donor_indices()] to see how donors are ranked;
#'   [derive_local_debiasing_weights()] and [derive_global_debiasing_weights()]
#'   for evaluating bias minimising weights.
#'
#' @name impute
NULL


#' @rdname impute
#' @export
impute_covariates <- function(x, recipients, donors, weights) {
  ordered_donor_x <- order_donor_covariates(
    x = x,
    recipients = recipients,
    donors = donors,
    k = ncol(weights)
  )

  imputation_wrapper <- function(idx) {
    crossprod(weights[idx, ], ordered_donor_x[[idx]]) |> as.numeric()
  }

  recipient_idx <- seq_along(ordered_donor_x)

  p <- ncol(x)

  vapply(recipient_idx, imputation_wrapper, numeric(p)) |> t()
}


#' @rdname impute
#' @export
impute_responses <- function(x, recipients, donors, weights, y) {
  ordered_donor_y <- order_donor_responses(
    x = x,
    recipients = recipients,
    donors = donors,
    y = y,
    k = ncol(weights)
  )

  rowSums(ordered_donor_y * weights)
}
