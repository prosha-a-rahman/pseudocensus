#' Evaluate the bias correction matrix for the pseudocensus OLS estimator
#'
#' @inheritParams merge_covariates
#'
#' @returns A square numeric matrix with `ncol(covariates)` rows and columns.
#'
#'
#' @export
#'
bias_correction_matrix <- function(covariates, recipients, donors, weights) {
  merged_covariates <- merge_covariates(
    covariates = covariates,
    recipients = recipients,
    donors = donors,
    weights = weights
  )

  A <- crossprod(covariates)

  B <- crossprod(covariates, merged_covariates)

  solve(B, A)
}
