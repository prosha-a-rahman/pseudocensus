#' Compute the ordinary least squares estimator for a population subset
#'
#' @param covariates A numeric matrix of rowwise covariates.
#' @param responses A numeric vector of length `nrow(covariates)`.
#' @param donors An integer vector in the row indices of `covariates` that
#'   identifies the population subset.
#'
#'
#' @returns
#' A numeric vector of length `ncol(covariates)` that is the OLS estimate
#' between `covariates[donors, ]` and `responses[donors]`.
#'
#' @export
#'
fit_sample_ols <- function(covariates, responses, donors) {
  sample_covariates <- covariates[donors, ]

  sample_responses <- responses[donors]

  .lm.fit(sample_covariates, sample_responses)$coefficients
}




#' Compute the ordinary least squares estimator for a partially nearest
#' neighbours imputed population
#'
#' @details
#' Compute the ordinary least squares estimator in which a population subset of
#' `recipients` has had their responses imputed by weighted nearest neighbours.
#'
#'
#' @inheritParams merge_responses
#'
#' @returns
#' A numeric vector of length `ncol(covariates)`.
#' @export
#'
fit_pseudocensus_ols <- function(covariates,
                                 responses,
                                 donors,
                                 recipients,
                                 weights) {
  merged_responses <- merge_responses(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    weights = weights,
    responses = responses
  )

  .lm.fit(covariates, merged_responses)$coefficients
}


#' Compute the bias corrected pseudocensus OLS estimator
#'
#' @details
#' Apply the bias correction matrix to the pseudocensus OLS estimator, removing
#' the bias induced by imputing the responses of `recipients` with weighted
#' nearest neighbours.
#'
#' @inheritParams merge_responses
#'
#' @returns
#' A numeric vector of length `ncol(covariates)`.
#'
#' @seealso [fit_pseudocensus_ols()] for the uncorrected estimator and
#'   [bias_correction_matrix()] for the correction applied here.
#'
#' @export
fit_bias_corrected_pseudocensus_ols <- function(covariates,
                                                responses,
                                                donors,
                                                recipients,
                                                weights) {
  pseudocensus_ols <- fit_pseudocensus_ols(
    covariates = covariates,
    responses = responses,
    donors = donors,
    recipients = recipients,
    weights = weights
  )

  correction_matrix <- bias_correction_matrix(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    weights = weights
  )

  drop(correction_matrix %*% pseudocensus_ols)
}


