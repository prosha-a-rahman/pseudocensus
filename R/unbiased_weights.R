calc_debiasing_weights_actual_sample <- function(
    locations,
    donors,
    recipients,
    n_weights,
    p = 2
) {

  n_population <- nrow(locations)

  dim_locations <- ncol(locations)

  A <- matrix(nrow = dim_locations^2, ncol = n_weights)

  for (i in 1:n_weights) {

    kth_nearest_location_matrix <- fetch_kth_nearest_location_matrix(
      recipients = recipients,
      donors = donors,
      locations = locations,
      p = p,
      k = k
    )

    A[, i] <- c(crossprod(locations, kth_nearest_location_matrix))

  }

  recipient_locations <- matrix(0, nrow = n_population, ncol = dim_locations)

  recipient_locations[recipients, ] <- locations[recipients, ]

  b <- c(crossprod(recipient_locations))

  qr.solve(A, b)

}
