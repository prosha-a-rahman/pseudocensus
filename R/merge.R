#' Replace only recipient values using weighted nearest neighbours
#'
#' @inheritParams impute_covariates
#'
#' @returns
#' * `merge_covariates()`: A numeric matrix of the same dimension as
#'   `covariates`, where only rows indexed by `recipients` are replaced with
#'   their nearest neighbours estimate. Rows indexed by `donors` are unchanged.
#'
#' * `merge_responses()`: A numeric vector of the same length as `responses`
#'   where only elements indexed by `recipients` are replaced with their nearest
#'   neighbours estimate. Rows indexed by `donors` remain unchanged.
#'
#' @export
merge_covariates <- function(covariates, recipients, donors, weights) {
  imputed_covariates <- impute_covariates(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    weights = weights
  )

  covariates[recipients, ] <- imputed_covariates

  covariates
}


#' @rdname merge_covariates
#' @export
merge_responses <- function(covariates, recipients, donors, weights, responses) {
  imputed_responses <- impute_responses(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    weights = weights,
    responses = responses
  )

  responses[recipients] <- imputed_responses

  responses
}
