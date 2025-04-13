
#' @export
fetch_kth_nearest_donor_location <- function(
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




#' @export
fetch_all_kth_nearest_donor_locations <- function(
    all_recipients,
    locations,
    donors,
    p = 2,
    k
) {

  fetch_kth_nearest_donor_location_wrapper <- function(recipient) {

    fetch_kth_nearest_donor_location(
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
    fetch_kth_nearest_donor_location_wrapper,
    numeric(dim_locations)
  ) |> t()

}


