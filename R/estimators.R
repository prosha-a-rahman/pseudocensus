#' Compute the OLS estimator for a population subset
#'
#' @param x A numeric matrix of rowwise predictors.
#' @param y A numeric vector of responses of length `nrow(x)`.
#' @param donors An integer vector in the row indices of `x` that identifies the
#'   population subset.
#'
#'
#' @returns
#' A numeric vector of length `ncol(x)`
#'
#' @export
#'
fit_sample_ols <- function(x, y, donors) {
  x <- prepend_ones(x)

  x_donors <- x[donors, , drop = FALSE]

  y_donors <- y[donors]

  stats::.lm.fit(x_donors, y_donors)$coefficients
}




#' Compute the ordinary least squares estimator for a partially nearest
#' neighbours imputed population
#'
#' @details
#' Compute the ordinary least squares estimator in which a population subset of
#' `recipients` has had their responses imputed by weighted nearest neighbours.
#'
#'
#' @inheritParams merge
#'
#' @returns
#' A numeric vector of length `ncol(x)`.
#' @export
#'
fit_pseudocensus_ols <- function(x, y, donors, recipients, weights) {
  x <- prepend_ones(x)

  y_merged <- merge_responses(
    x = x,
    recipients = recipients,
    donors = donors,
    weights = weights,
    y = y
  )

  stats::.lm.fit(x, y_merged)$coefficients
}


#' Compute the bias corrected pseudocensus OLS estimator
#'
#' @details
#' Apply the bias correction matrix to the pseudocensus OLS estimator, removing
#' the bias induced by imputing the responses of `recipients` with weighted
#' nearest neighbours.
#'
#' @inheritParams merge
#'
#' @returns
#' A numeric vector of length `ncol(x)`.
#'
#' @seealso [fit_pseudocensus_ols()] for the uncorrected estimator and
#'   [debiasing_coefficient()] for the correction applied here.
#'
#' @export
fit_debiased_pseudocensus_ols <- function(x, y, donors, recipients, weights) {
  x <- prepend_ones(x)

  pseudocensus_ols <- fit_pseudocensus_ols(
    x = x,
    y = y,
    donors = donors,
    recipients = recipients,
    weights = weights
  )

  correction_matrix <- debiasing_coefficient(
    x = x,
    recipients = recipients,
    donors = donors,
    weights = weights
  )

  drop(correction_matrix %*% pseudocensus_ols)
}


