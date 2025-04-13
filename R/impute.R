#' Impute a location or response based on weighted k nearest neighbours
#'
#' Supply synthetic approximations for the locations or responses of a
#' population subset using weighted \eqn{k} nearest neighbours.
#'
#' @inheritParams order_donors
#' @param all_recipients Index vector that identifies the recipients whose
#'   responses are to be imputed.
#' @param weights Numeric vector of weights used for calculating the weighted
#'   average response or location.
#' @param responses Numeric vector of responses, potentially containing
#'   missing `NA` values.
#' @param donors Index vector of donors. If `donors = NULL`, then the donors are
#'  identified by the indexes of non-`NA` values in `responses`.
#'
#' @details
#' `impute_all_recipient_responses()` generates responses for each element in
#' `all_recipients` using a weighted average of proximal donor responses.
#' `impute_all_recipient_location()` provides an analogue for the location of
#' `all_recipients`.
#'
#' @return
#' * `impute_all_recipient_responses()`: A length `length(all_recipients)`
#'   numeric vector of the weighted averages of the `length(weights)` nearest
#'   donors for each recipient.
#'
#' * `impute_recipient_location()`: A numeric vector that approximates the
#'   location of the recipient using the weighted average of the location of the
#'   `length(weights)` nearest donors. The length of the output is the
#'   dimension of the location—`ncol(locations)`.
#'
#' * `impute_all_recipient_locations()`: A
#'   `length(all_recipients) × ncol(locations)` matrix where each row contains
#'   the approximated location for a recipient.
#'
#' @seealso [order_donors()] for how distances are calculated.
#'
#' @export
impute_all_recipient_responses <- function(
    all_recipients,
    locations,
    responses,
    weights,
    donors = NULL,
    p = 2
) {

  if (is.null(donors)) { donors <- identify_donors(responses) }

  n_weights <- length(weights)

  order_donors_wrapper <- function(recipient) {

    order_donors(
      recipient = recipient,
      donors = donors,
      locations = locations,
      p = p,
      k = n_weights
    )

  }

  ordered_donors <- vapply(
    all_recipients,
    order_donors_wrapper,
    integer(n_weights)
  )

  donor_responses <- responses[ordered_donors]

  dim(donor_responses) <- dim(ordered_donors)

  crossprod(weights, donor_responses)[1, ]

}




#' @rdname impute_all_recipient_responses
#' @export
impute_recipient_location <- function(
    recipient,
    donors,
    locations,
    weights,
    p = 2
) {

  n_weights <- length(weights)

  ordered_donors <- order_donors(
    recipient = recipient,
    donors = donors,
    locations = locations,
    p = p,
    k = n_weights
  )

  ordered_donor_locations <- locations[ordered_donors, ]

  crossprod(weights, ordered_donor_locations)[1, ]

}




#' @rdname impute_all_recipient_responses
#' @export
impute_all_recipient_locations <- function(
    all_recipients,
    donors,
    locations,
    weights,
    p = 2
) {

  impute_recipient_location_wrapper <- function(recipient) {

    impute_recipient_location(
      recipient = recipient,
      donors = donors,
      locations = locations,
      weights = weights,
      p = p
    )

  }

  dim_locations <- ncol(locations)

  vapply(
    all_recipients,
    impute_recipient_location_wrapper,
    numeric(dim_locations)
  ) |> t()

}
