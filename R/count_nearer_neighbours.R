#' Count the number of points closer to the recipient
#'
#' `count_nearer_neighbours()` counts the number of points whose covariates are
#'   relatively closer to the `recipient` covariate than the `donor` covariate
#'   is to the  `recipient` covariate.
#
#' @param recipient The index of the recipient
#' @param donor The index of the prospective donor
#' @param covariates A design matrix of covariates listed row-wise
#' @param p The order of the \eqn{\ell_{p}}-norm. The default value `p = 2`
#'            yields the Euclidean metric.
#'
#' @return A non-negative integer.``li
#' @export
#'
count_nearer_neighbours <- function(donor, recipient, covariates, p = 2) {

  if (recipient == donor) {

    return(0)

  }

  recipient_covariate <- covariates[recipient, ]

  donor_covariate <- covariates[donor, , drop = FALSE]

  donor_distance <- lp_distance(recipient_covariate, donor_covariate, p)

  neighbour_distances <- lp_distance(recipient_covariate, covariates, p)

  sum(neighbour_distances <= donor_distance) - 2

}
