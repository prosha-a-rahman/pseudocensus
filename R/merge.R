#' Merge imputed and observed responses and covariates
#'
#' Combine observed and imputed responses (or covariates) into a single object.
#'
#' @inheritParams impute_response
#' @param recipients Integer vector of recipients whose responses, possibly
#'   missing, require imputation
#'
#' @return
#' * `merge_responses()`: returns `responses`, except where the responses for
#'   the recipient indices have been replaced with their weighted \eqn{k}
#'   nearest neighbours imputed equivalent.
#' * `merge_covariates()`: returns `covariates`, except where the rows for the
#'   recipient indices have been replaced with their weighted \eqn{k}
#'   nearest neighbours imputed equivalent.
#'
#' @export
merge_responses <- function(
  covariates,
  responses,
  weights,
  p = 2,
  recipients = NULL,
  donors = NULL
) {

  if (is.null(recipients)) { recipients <- identify_recipients(responses) }

  if (is.null(donors)) { donors <- identify_donors(responses) }

  imputation_wrapper <- function(recipient) {
    impute_response(
      recipient = recipient,
      covariates = covariates,
      responses = responses,
      weights = weights,
      donors = donors,
      p = p
    )
  }

  imputed_values <- vapply(recipients, imputation_wrapper, double(1))

  # responses[recipients] <- vapply(
  #   recipients,
  #   impute_response,
  #   double(1),
  #   donors,
  #   covariates,
  #   responses,
  #   weights,
  #   p
  # )

  responses[recipients] <- imputed_values

  responses

}

#' @rdname merge_responses
#' @export
merge_covariates <- function(
  covariates,
  weights,
  p = 2,
  recipients = NULL,
  donors = NULL,
  responses = NULL
) {

  if (is.null(recipients)) { recipients <- identify_recipients(responses) }

  if (is.null(donors)) { donors <- identify_donors(responses) }

  n_parameters <- ncol(covariates)

  # covariates[recipients, ] <- vapply(
  #   recipients,
  #   impute_covariate,
  #   double(number_parameters),
  #   donors,
  #   covariates,
  #   weights
  # ) |> t()

  imputation_wrapper <- function(recipient) {
    impute_covariate(
      recipient = recipient,
      donors = donors,
      covariates = covariates,
      weights = weights,
      p = p
    )
  }

  imputed_values <- vapply(recipients, imputation_wrapper, double(n_parameters))

  covariates[recipients, ] <- t(imputed_values)

  covariates

}
