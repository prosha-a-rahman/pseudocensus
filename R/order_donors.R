#' Rank donors by distance
#'
#' Order donors according to their proximity to a recipient.
#'
#' @param recipient Index of the recipient.
#' @param donors Index vector of donors.
#' @param locations Numeric matrix identifying rowwise location of indexes.
#' @param p Number >= 1 (or `Inf`) specifying which `p`-norm to use for
#'   distance calculations. The default is `p = 2` (Euclidean distances).
#' @param k Positive integer indicating if only the `k` nearest donors should
#' be returned. If `k = NULL` (default), then all `donors` are ranked.
#'
#' @details
#' `order_donors()` ascendingly ranks the set of `donors` by their proximity to
#' the `recipient`. Distances are calculated with respect to the rowwise values
#' in `locations`. If `k` is specified, then only the `k` nearest donors to the
#' `recipient` are returned, otherwise the entire set of donors are ranked.
#'
#' For finite `p`, the `p`-norm distance between vectors \eqn{x} and \eqn{y} is
#' \deqn{||x - y||_p = \left(\sum_{i=1}^{n} |x_i - y_i|^p\right)^{1/p}.}
#'
#' For `p = Inf`, the distance is the maximum of the absolute differences:
#' \deqn{||x - y||_{\infty} = \max_{i} |x_i - y_i|.}
#'
#' Common `p`-norm metrics:
#' * `p = 1`: Manhattan/taxicab distance (sum of absolute differences)
#' * `p = 2`: Euclidean distance (straight-line distance)
#' * `p = Inf`: Chebyshev distance (maximum absolute difference)
#'
#' `order_donors()` uses `Rfast::dista()` implementation for fast distance
#' calculations.
#'
#' @return An integer vector containing a permuted subset of `donors`, ordered
#' according to their ascendingly ranked distance from the `recipient`. If `k`
#' is specified, only the `k` closest donors are returned, otherwise
#' `order_donors()` returns a permutation of `donors`.
#'
#'
#' @seealso [Rfast::dista()] for the underlying distance calculation.
#'
#' @export
order_donors <- function(recipient, donors, locations, p = 2, k = NULL) {

  if (is.null(k)) { k <- length(donors) }

  recipient_location <- locations[recipient, , drop = FALSE]

  donor_locations <- locations[donors, , drop = FALSE]

  ranked_donors <- Rfast::dista(
    xnew = recipient_location,
    x = donor_locations,
    type = "minkowski",
    k = k,
    index = TRUE,
    p = p
  ) |> as.integer()

  donors[ranked_donors]

}


