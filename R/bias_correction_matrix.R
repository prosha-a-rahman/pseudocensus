#' Evaluate the bias correction matrix for the pseudocensus OLS estimator
#'
#' @inheritParams merge
#'
#' @returns A square numeric matrix with `ncol(x)` rows and columns.
#'
#'
#' @export
#'
debiasing_coefficient <- function(x, recipients, donors, weights) {
  x_merged <- merge_covariates(
    x = x,
    recipients = recipients,
    donors = donors,
    weights = weights
  )

  A <- crossprod(x)

  B <- crossprod(x, x_merged)

  solve(B, A)
}
