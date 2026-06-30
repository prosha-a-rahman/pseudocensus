#' Rank donors by their distance to each recipient
#'
#' Rank a candidate set of points by their proximity to a set of reference
#' points and return their associated indices, covariates, or responses.
#'
#' @details
#' Distances are computed between the rows of `x`, each of which is a
#' `ncol(x)`-dimensional covariate vector defining the matching variables. The
#' row-index vectors `donors` and `recipients` respectively identify the
#' candidate and reference points in `x`. Donors are independently ranked by
#' their proximity to each recipient using the Euclidean metric, and the nearest
#' `k` are returned. The associated set of indices, covariates, or responses for
#' the ranked donors are returned depending on which of `order_donor_*()` is
#' called.
#'
#' It is possible for `donors` and `recipients` to overlap, in which case the
#' nearest donor for overlapping recipients will be the recipients themselves.
#' Distance ties are resolved by the underlying [nabor::knn()] search.
#'
#' @param x A numeric matrix of rowwise covariates used as matching variables
#'   for distance calculations.
#' @param recipients An integer vector of row indices in `x` identifying the
#'   reference (query) points.
#' @param donors An integer vector of row indices in `x` identifying the
#'   candidate points.
#' @param k An integer in `1:length(donors)` giving the number of ranked donors
#'   to return per recipient. Defaults to `length(donors)`, which ranks all of
#'   `donors`.
#' @param y A length `nrow(x)` numeric vector of responses indexed consistently
#'   with the rows of `x`. Only used by `order_donor_responses()`.
#'
#' @returns
#' Writing `n_recipients = length(recipients)` and `p = ncol(x)`:
#' * `order_donor_indices()`: an `n_recipients`-by-`k` index-valued matrix where
#'   each row ranks `donors` by their proximity to the corresponding recipient
#'   index and returns the `k` nearest donor indices.
#'
#' * `order_donor_covariates()`: a length `n_recipients` list of `k`-by-`p`
#'   numeric matrices. Each element in the list is a rowwise permuted subset of
#'   `x`, where the rowwise permutations correspond to the index permutation
#'   in `order_donor_indices()`.
#'
#' * `order_donor_responses()`: an `n_recipients`-by-`k` numeric matrix where
#'   each row is a permuted subset of `y` corresponding to the index permutation
#'   in `order_donor_indices()`.
#'
#' @seealso [nabor::knn()] for the underlying nearest-neighbour search.
#'
#' @name order_donor
NULL


#' @rdname order_donor
#' @export
order_donor_indices <- function(x, recipients, donors, k = length(donors)) {
  # drop = FALSE maintains matrix structure required for nabor::knn() when
  # length(recipients) = 1 or length(donors) = 1
  x_recipients <- x[recipients, , drop = FALSE]
  x_donors <- x[donors, , drop = FALSE]

  ordered_donor_idx <- nabor::knn(
    data = x_donors,
    query = x_recipients,
    k = k
  )$nn.idx

  # Map ranked donor sub-indices to actual donor indexes
  matrix(
    donors[ordered_donor_idx],
    nrow = nrow(ordered_donor_idx),
    ncol = ncol(ordered_donor_idx)
  )
}


#' @rdname order_donor
#' @export
order_donor_covariates <- function(x, recipients, donors, k = length(donors)) {
  ordered_donors <- order_donor_indices(
    x = x,
    recipients = recipients,
    donors = donors,
    k = k
  )

  lapply(
    asplit(ordered_donors, 1),
    function(idx) x[idx, , drop = FALSE]
  )
}



#'
#' @rdname order_donor
#' @export
order_donor_responses <- function(x,
                                  recipients,
                                  donors,
                                  y,
                                  k = length(donors)) {
  ordered_donors <- order_donor_indices(
    x = x,
    recipients = recipients,
    donors = donors,
    k = k
  )

  matrix(
    y[ordered_donors],
    nrow = nrow(ordered_donors),
    ncol = ncol(ordered_donors)
  )
}
