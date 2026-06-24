#' Order donors by their distance from recipients
#'
#' Rank a candidate set of points by their distance from a set of reference
#' points and return their associated indices, covariates, or responses.
#'
#' @details
#' `order_donor_*()` defines a family of functions that ranks a candidate set
#' of points by their distance to a set of reference points. The points are
#' the `ncol(covariates)`-dimensional vectors defined by the rows of
#' `covariates`. The disjoint index vectors `recipients` and `donors`
#' respectively define the indexes of the reference and candidate points. Donors
#' are independently ranked by proximity for each recipient, with respect to the
#' Euclidean metric via [nabor::knn()]. The associated ordered set of indices,
#' covariates, or responses are returned depending on which of `order_donor_*()`
#' is called.
#'
#' @param covariates A rowwise numeric matrix of points.
#' @param recipients An integer vector of row indices in `covariates`
#'   identifying the reference points. `recipients` must be disjoint to `donors`.
#' @param donors An integer vector of row indices in `covariates` identifying
#'   the candidate points.  `donors` must be disjoint to `recipients`.
#' @param k An integer scalar between 1 and `length(donors)` giving the number
#'   of donors to return per recipient. The default sets `k = length(donors)`
#'   and rank the full list of donors for each recipient.
#'
#' @returns
#' Writing `n_recipients = length(recipients)` and `p = ncol(covariates)`:
#'
#' * `order_donor_indices()`: an `n_recipients`-by-`k` integer matrix where
#'   each row is a permuted subset of `donors` where, for the corresponding
#'   recipient index, the permutation sorts the donors by their proximity from
#'   the recipient.
#'
#' * `order_donor_covariates()`: a length `n_recipients` list of `k`-by-`p`
#'   numeric matrices. Each element in the list is a rowwise permuted subset of
#'   `covariates` corresponding to the index permutation in
#'   `order_donor_indices()`.
#'
#' * `order_donor_responses()`: an `n_recipients`-by-`k` numeric matrix where
#'   each row is a permuted subset of `responses` corresponding to the index
#'   permutation in `order_donor_indices()`.
#'
#' @seealso [nabor::knn()] for the underlying nearest-neighbour search.
#'
#' @export
order_donor_indices <- function(covariates, recipients, donors, k = NULL) {
  # drop = FALSE maintains matrix structure required for nabor::knn() when
  # length(recipients) = 1 or length(donors) = 1
  recipient_covariates <- covariates[recipients, , drop = FALSE]
  donor_covariates <- covariates[donors, , drop = FALSE]

  if (is.null(k)) k <- length(donors)

  ordered_donor_idx <- nabor::knn(
    data = donor_covariates,
    query = recipient_covariates,
    k = k
  )$nn.idx

  matrix(
    donors[ordered_donor_idx],
    nrow = nrow(ordered_donor_idx),
    ncol = ncol(ordered_donor_idx)
  )
}


#' @rdname order_donor_indices
#' @export
order_donor_covariates <- function(covariates, recipients, donors, k = NULL) {
  ordered_donors <- order_donor_indices(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    k = k
  )

  lapply(
    asplit(ordered_donors, 1),
    function(idx) covariates[idx, , drop = FALSE]
  )
}


#' @param responses A numeric vector of responses indexed consistently with the
#'   rows of `covariates`.
#'
#' @rdname order_donor_indices
#' @export
order_donor_responses <- function(covariates,
                                  recipients,
                                  donors,
                                  responses,
                                  k = NULL) {
  ordered_donors <- order_donor_indices(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    k = k
  )

  matrix(
    responses[ordered_donors],
    nrow = nrow(ordered_donors),
    ncol = ncol(ordered_donors)
  )
}
