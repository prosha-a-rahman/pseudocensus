#' Impute a response based on weighted k nearest neighbours
#'
#' `impute_response()` constructs a synthetic response for the `recipient` by
#'   calculating a weighted average of proximal `donor` responses—using \eqn{k}
#'   nearest neighbours. `impute_covariate()` supplies the equivalent for the
#'   covariate vector of the `recipient`.
#'
#' @inheritParams order_donors
#' @param weights Numeric vector of weights used for calculating the weighted
#'   average.
#' @param responses Numeric vector of responses, pontentitally containing
#'   missing `NA` values
#'
#' @return
#' * `impute_response()`: A numeric value of the weighted average of the
#'   responses of the `length(weights)` nearest donors.
#' * `impute_covariate()`: A numeric vector of the weighted average of the
#'   covariates of the `length(weights)` nearest donors.
#'
#' @export
impute_response <- function(
  recipient,
  covariates,
  responses,
  weights,
  donors = NULL,
  p = 2
) {

  if (is.null(donors)) { donors <- identify_donors(responses) }

  ordered_donors <- order_donors(
    recipient = recipient,
    donors = donors,
    covariates = covariates,
    p = p
  )

  ordered_donor_responses <- responses[ordered_donors[seq_along(weights)]]

  sum(ordered_donor_responses * weights)

}

#' @rdname impute_response
#' @export
impute_covariate <- function(recipient, donors, covariates, weights, p = 2) {

  ordered_donors <- order_donors(
    recipient = recipient,
    donors = donors,
    covariates = covariates,
    p = p
  )

  ordered_donor_covariates <- covariates[ordered_donors[seq_along(weights)], ]

  colSums(ordered_donor_covariates * weights)

}
