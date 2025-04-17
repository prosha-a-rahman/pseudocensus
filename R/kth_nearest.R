#' Extract proximally ordered locations
#'
#' Find the k-th nearest location in a set to a specific point.
#'
#' @inheritParams impute_response
#' @inheritParams order_donors
#'
#' @details
#' The function extracts the location of the `k`-th nearest donor for each
#' `recipient`. This is functionally equivalent to [impute_recipient_location()]
#' with `weights = c(rep(0, k), 1)`.
#'
#'
#' @returns
#' * `fetch_kth_nearest_donor_location()` : A row-element of `locations` whose
#'   index is in `donors`.
#'
#' * `fetch_all_kth_nearest_donor_locations()`: A matrix with
#'   `length(all_recipients)` rows, where each row maps to a row in `locations`
#'   whose indexes correspond to elements of `donors`.
#'
#' @seealso [order_donors()] for how distances are calculated.
#' @seealso [impute_location()] for functional equivalent.
#'
#' @export
fetch_kth_nearest_location <- function(
    recipients,
    locations,
    donors,
    p = 2,
    k
) {

  if (length(recipients) == 1) {

    return(
      fetch_kth_nearest_location_single_recipient(
        recipient = recipients,
        locations = locations,
        donors = donors,
        p = p,
        k = k
      )
    )

  }


  fetch_kth_nearest_location_single_recipient_wrapper <- function(recipient) {

    fetch_kth_nearest_location_single_recipient(
      recipient = recipient,
      locations = locations,
      donors = donors,
      p = 2,
      k = k
    )

  }

  dim_locations <- ncol(locations)

  vapply(
    all_recipients,
    fetch_kth_nearest_location_single_recipient_wrapper,
    numeric(dim_locations)
  ) |> t()

}




fetch_kth_nearest_location_single_recipient <- function(
    recipient,
    locations,
    donors,
    p = 2,
    k
) {

  ordered_donors <- order_donors(
    recipient = recipient,
    donors = donors,
    locations = locations,
    p = p,
    k = k
  )

  kth_nearest_donor <- ordered_donors[k]

  locations[kth_nearest_donor, ]

}




fetch_kth_nearest_location_matrix <- function(
    recipients,
    donors,
    locations,
    p = 2,
    k
) {

  n_population <- nrow(locations)

  dim_locations <- ncol(locations)

  kth_imputed_location_matrix <- matrix(
    0, ncol = dim_locations, nrow = n_population
  )

  kth_imputed_locations[recipients, ] <- fetch_kth_nearest_location(
    recipients = recipients,
    locations = locations,
    donors = donors,
    p = p,
    k = k
  )

  kth_imputed_locations

}


