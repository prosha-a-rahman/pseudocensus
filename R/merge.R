#' Impute recipient values in place by weighted nearest neighbours
#'
#' Impute the covariates or responses of `recipients` from weighted combinations
#' of their proximally ranked donors covariates or responses, and splice the
#' estimates back into the full data object in place of the original recipient
#' entries.
#'
#' @inheritParams impute
#'
#' @returns
#' * `merge_covariates()`: A numeric matrix with the same dimensions as
#'   `x`, in which rows indexed by `recipients` are replaced with their weighted
#'   nearest neighbours estimate and all other rows are unchanged.
#'
#' * `merge_responses()`: A numeric vector with the same length as `y`,
#'   where elements indexed by `recipients` are replaced with their weighted
#'   nearest neighbours estimate and all other elements are unchanged.
#'
#' @seealso [order_donor_indices()] to see how donors are ranked;
#'   [derive_local_debiasing_weights()] and [derive_global_debiasing_weights()]
#'   for constructing bias minimising weights.
#'
#' @name merge
NULL

#' @rdname merge
#' @export
merge_covariates <- function(x, recipients, donors, weights) {
  imputed_recipient_x <- impute_covariates(
    x = x,
    recipients = recipients,
    donors = donors,
    weights = weights
  )

  x[recipients, ] <- imputed_recipient_x

  x
}


#' @rdname merge
#' @export
merge_responses <- function(x, recipients, donors, weights, y) {
  imputed_recipient_y <- impute_responses(
    x = x,
    recipients = recipients,
    donors = donors,
    weights = weights,
    y = y
  )

  y[recipients] <- imputed_recipient_y

  y
}
