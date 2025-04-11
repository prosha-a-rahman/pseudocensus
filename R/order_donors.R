#' Rank the donors by their distance from a recipient
#'
#' Order donors according to their proximity to a recipient using the `p`-norm
#' distance between their respective covariates.
#'
#' @param recipient Integer index identifying the recipient in `covariates`.
#' @param donors Integer vector of indices identifying the donors in `covariates`.
#' @param covariates Numeric matrix where each row is an individual's covariates.
#' @param p A double >= 1 (or `Inf`) specifying which `p`-norm to use for
#'   distance calculations. The deafult is `p = 2` (Euclidean distances).
#'
#' @details
#' `order_donors()` ascendingly ranks the donors according to their proximity to
#' a recipient based on their covariates with respect to the vector `p`-norm.
#'
#'
#' @return A permutation of the integer donor indices according to their
#'   ascendingly ranked proximity to the recipient.
#'
#' @export
order_donors <- function(recipient, donors, covariates, p = 2) {

  recipient_covariate <- covariates[recipient, , drop = FALSE]

  donor_covariates <- covariates[donors, , drop = FALSE]

  # donor_distances <- lp_distance(
  #   reference_point = recipient_covariate,
  #   point_matrix =  donor_covariates,
  #   p = p
  # )

  ranked_donors <- Rfast::dista(
    xnew = recipient_covariate,
    x = donor_covariates,
    type = "minkowski",
    p = p,
    k = length(donors),
    index = TRUE
  ) |> as.integer()

  donors[ranked_donors]

}


