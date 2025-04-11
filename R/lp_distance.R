#' Compute distances between points using p-norm metrics
#'
#' Calculates the distance of a set of points to a reference point with respect
#'   to the vector `p`-norm.
#'
#' @param reference_point A numeric vector representing the reference point from
#'   which distances are measures
#' @param point_matrix A numeric matrix with `length(reference_point)` rows,
#'   where each column represents a point in the same vector space as
#'   `reference_point`.
#' @param p A double >= 1 (or `Inf`) specifying which `p`-norm to use for
#'   distance calculations. The deafult is `p = 2` (Euclidean distances).
#'
#' @details
#' `lp_distance()` computes the distance between the specified `reference_point`
#' and each point represented by the columns of `point_matrix` using the vector
#' `p`-norm.
#'
#' For finite `p`, the `p`-norm distance between vectors \eqn{x} and \eqn{y} is
#' \deqn{||x - y||_p = \left(\sum_{i=1}^{n} |x_i - y_i|^p\right)^{1/p}.}
#'
#' For `p = Inf`, the distance is the maximum of the absolute differences:
#' \deqn{||x - y||_{\infty} = \max_{i} |x_i - y_i|}
#'
#' Common p-norm metrics:
#' * `p = 1`: Manhattan/Taxicab distance (sum of absolute differences)
#' * `p = 2`: Euclidean distance (straight-line distance)
#' * `p = Inf`: Chebyshev distance (maximum absolute difference)
#'
#' When `reference_point` is the zero vector, `lp_distance()` calculates the
#' `p`-norm of each point in the `point_matrix`.
#'
#' @return A length `ncol(point_matrix)` numeric vector of non-negative distances.
#'
#' @export
lp_distance <- function(reference_point, point_matrix, p = 2) {

  absolute_difference_matrix <- abs(sweep(point_matrix, 1, reference_point))

  if (is.finite(p)) { colSums(absolute_difference_matrix^p)^(1 / p) }

  else if (p == Inf) { apply(absolute_difference_matrix, 2, max) }

}
