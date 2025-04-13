calc_debiasing_weights_actual_sample <- function(
    locations,
    donors,
    all_recipients,
    n_weights,
    p = 2
) {

  n_population <- nrow(locations)

  dim_locations <- ncol(locations)

  matrix_to_be_inverted <- matrix(nrow = dim_locations^2, ncol = n_weights)

  for (i in 1:n_weights) {

    temp_matrix <- matrix(0, nrow = n_population, ncol = dim_locations)

    temp_matrix[all_recipients, ] <- fetch_all_kth_nearest_donor_locations(
      all_recipients = all_recipients,
      locations = locations,
      donors = donors,
      p = p,
      k = i
    )

    matrix_to_be_inverted[, i] <- c(crossprod(locations, temp_matrix))

  }

  recipient_locations <- matrix(0, nrow = n_population, ncol = dim_locations)

  recipient_locations[all_recipients, ] <- locations[all_recipients, ]

  target_matrix <- c(crossprod(recipient_locations, recipient_locations))

  inverted_matrix <- MASS::ginv(matrix_to_be_inverted)

  inverted_matrix %*% target_matrix

}
