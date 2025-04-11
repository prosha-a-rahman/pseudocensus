#' Compute distances of a set of points from a translated origin
#'
#' Computes the distance between a vector (`origin`) and columns of a matrix
#'   (`design_matrix`) with respect to the metric induced by the vector `p`-norm.
#'
#' @param origin A numeric vector.
#' @param design_matrix A numeric matrix with `length(origin)` rows.
#' @param p A number greater than, or equal to, one (1) (possibly `Inf`).
#'
#' @details
#' Each column in `design_matrix` represents a point in the same
#'   `length(origin)`-dimensional Euclidean space in which `origin` is defined.
#'   The function `lp_distance()` computes the distance of each of these points
#'   from `origin` with respect to the usual vector `p`-norm. If `origin` is the
#'   usual Euclidean origin \eqn{\mathbf{0}}, then `lp_distance()` simply
#'   computes the usual vector `p`-norm. Some commonly used norms include:
#'
#'   * `p = 1` the Manhattan norm.
#'   * `p = 2` the usual Euclidean norm.
#'   * `p = Inf` the Chebyshev/uniform/max norm.
#'
#' @return A non-negative numeric vector of distances.
#' @export
lp_distance <- function(origin, design_matrix, p = 2) {

  absolute_difference_matrix <- abs(design_matrix - origin)

  if (is.finite(p)) {

    colSums(absolute_difference_matrix^p)^(1 / p)

  } else if (p == Inf) {

    apply(absolute_difference_matrix, 2, max)

  }

}
