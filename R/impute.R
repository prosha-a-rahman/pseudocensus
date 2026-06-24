#' Impute covariates or responses using weighted nearest neighbours
#'
#' Construct synthetic covariates or responses for a set of recipients from
#' weighted combinations of their nearest donor values.
#'
#' @inheritParams order_donor_indices
#' @param weights A numeric matrix of imputation weights used for weighted
#'   nearest neighbours imputation. `weights` corresponds rowwise to `recipients`
#'   and columnwise to the rank of the proximally sorted donors. The number of
#'   weights cannot exceed the number of donors, i.e.,
#'   `ncol(weights) <= length(donors)`.
#'
#' @returns
#' Writing `n_recipients = length(recipients)` and `p = ncol(covariates)`:
#'
#' * `impute_covariates()`: an `n_recipients`-by-`p` numeric matrix, where each
#'   row is the the weighted combination of nearest donor covariates
#'   corresponding to `recipients`.
#'
#' * `impute_responses()`: a numeric vector of length `n_recipients`, where each
#'   element is the weighted combination of the nearest donor responses
#'   corresponding to `recipients`.
#'
#' @seealso [order_donor_indices()] for how donors are ranked;
#'   [derive_local_debiasing_weights()] and [derive_global_debiasing_weights()]
#'   for evaluating bias minimising weights.
#'
#' @export
impute_covariates <- function(covariates, recipients, donors, weights) {
  ordered_donor_covariates <- order_donor_covariates(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    k = ncol(weights)
  )

  imputation_wrapper <- function(idx) {
    crossprod(weights[idx, ], ordered_donor_covariates[[idx]]) |> as.numeric()
  }

  recipient_idx <- seq_along(ordered_donor_covariates)

  p <- ncol(covariates)

  vapply(recipient_idx, imputation_wrapper, numeric(p)) |> t()
}


#' @rdname impute_covariates
#' @export
impute_responses <- function(covariates, recipients, donors, weights, responses) {
  ordered_donor_responses <- order_donor_responses(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    responses = responses,
    k = ncol(weights)
  )

  rowSums(ordered_donor_responses * weights)
}
