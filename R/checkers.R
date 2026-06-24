# Internal argument checkers
#
# Each `check_*()` validates one argument and, on failure, aborts with an
# informative {cli} message attributed to the calling function via `call`.
# On success it returns its first argument invisibly so calls can be chained.
#
# The index checkers (`check_recipient()`, `check_donors()`, `check_rank()`,
# `check_recipients()`) assume `covariates` / `donors` have already been
# validated, e.g. by calling `check_covariates()` first.


#' Check that `covariates` is a numeric matrix
#'
#' @param covariates Object to validate.
#' @param call Environment used to attribute the error.
#' @returns `covariates`, invisibly.
#' @noRd
check_covariates <- function(covariates, call = parent.frame()) {
  if (!is.matrix(covariates) || !is.numeric(covariates)) {
    cli::cli_abort(
      "{.arg covariates} must be a numeric matrix, not {.obj_type_friendly {covariates}}.",
      call = call
    )
  }

  invisible(covariates)
}


#' Check that `recipient` is a single row index of `covariates`
#'
#' @param recipient Object to validate.
#' @param covariates A numeric matrix.
#' @param call Environment used to attribute the error.
#' @returns `recipient`, invisibly.
#' @noRd
check_recipient <- function(recipient, covariates, call = parent.frame()) {
  n_row <- nrow(covariates)

  valid <- is.numeric(recipient) &&
    length(recipient) == 1L &&
    !is.na(recipient) &&
    recipient == as.integer(recipient) &&
    recipient >= 1L &&
    recipient <= n_row

  if (!valid) {
    cli::cli_abort(
      c(
        "{.arg recipient} must be a single row index of {.arg covariates}.",
        "i" = "Valid recipients are the integers 1 to {n_row}.",
        "x" = "{.arg recipient} is {.val {recipient}}."
      ),
      call = call
    )
  }

  invisible(recipient)
}


#' Check that every element of `donors` is a row index of `covariates`
#'
#' @param donors Object to validate.
#' @param covariates A numeric matrix.
#' @param call Environment used to attribute the error.
#' @returns `donors`, invisibly.
#' @noRd
check_donors <- function(donors, covariates, call = parent.frame()) {
  n_row <- nrow(covariates)

  if (!is.numeric(donors) || length(donors) == 0L) {
    cli::cli_abort(
      "{.arg donors} must be a non-empty numeric vector of row indices of {.arg covariates}.",
      call = call
    )
  }

  invalid <- is.na(donors) |
    donors != as.integer(donors) |
    donors < 1L |
    donors > n_row

  if (any(invalid)) {
    cli::cli_abort(
      c(
        "Every element of {.arg donors} must be a row index of {.arg covariates}.",
        "i" = "Valid indices are the integers 1 to {n_row}.",
        "x" = "Invalid {.arg donors}: {.val {donors[invalid]}}."
      ),
      call = call
    )
  }

  invisible(donors)
}


#' Check that `proximity_rank` indexes into `donors`
#'
#' @param proximity_rank Object to validate.
#' @param donors An integer vector of donor indices.
#' @param call Environment used to attribute the error.
#' @returns `proximity_rank`, invisibly.
#' @noRd
check_proximity_rank <- function(proximity_rank, donors, call = parent.frame()) {
  n_donors <- length(donors)

  if (!is.numeric(proximity_rank) || length(proximity_rank) == 0L) {
    cli::cli_abort(
      "{.arg proximity_rank} must be a non-empty numeric vector.",
      call = call
    )
  }

  invalid <- is.na(proximity_rank) |
    proximity_rank != as.integer(proximity_rank) |
    proximity_rank < 1L |
    proximity_rank > n_donors

  if (any(invalid)) {
    cli::cli_abort(
      c(
        "Every element of {.arg proximity_rank} must index into {.arg donors}.",
        "i" = "Valid donor ranks are the integers 1 to {n_donors}.",
        "x" = "Invalid {.arg proximity_rank}: {.val {proximity_rank[invalid]}}."
      ),
      call = call
    )
  }

  invisible(proximity_rank)
}


#' Check that every element of `recipients` is a row index of `covariates`
#'
#' @param recipients Object to validate.
#' @param covariates A numeric matrix.
#' @param call Environment used to attribute the error.
#' @returns `recipients`, invisibly.
#' @noRd
check_recipients <- function(recipients, covariates, call = parent.frame()) {
  n_row <- nrow(covariates)

  if (!is.numeric(recipients) || length(recipients) == 0L) {
    cli::cli_abort(
      "{.arg recipients} must be a non-empty numeric vector of row indices of {.arg covariates}.",
      call = call
    )
  }

  invalid <- is.na(recipients) |
    recipients != as.integer(recipients) |
    recipients < 1L |
    recipients > n_row

  if (any(invalid)) {
    cli::cli_abort(
      c(
        "Every element of {.arg recipients} must be a row index of {.arg covariates}.",
        "i" = "Valid indices are the integers 1 to {n_row}.",
        "x" = "Invalid {.arg recipients}: {.val {recipients[invalid]}}."
      ),
      call = call
    )
  }

  invisible(recipients)
}
