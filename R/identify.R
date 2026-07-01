#' Identify donors and recipients in a population
#'
#' Identify donor indices for population elements with observed responses, and
#' recipient indices for those with missing responses requiring imputation.
#'
#' @param y A numeric vector, possibly containing `NA` values.
#'
#' @returns
#' * `identify_donors()`: Integer vector of indices where responses are not `NA`.
#'
#' * `identify_recipients()`: Integer vector of indices where responses are `NA`.
#'
#' @examples
#' # Donors have observed responses; recipients are missing (NA)
#' responses <- c(10, NA, 30, NA, 50)
#'
#' identify_donors(responses)
#' identify_recipients(responses)
#'
#' @name identify
NULL


#' @rdname identify
#' @export
identify_donors <- function(y) which(!is.na(y))


#' @rdname identify
#' @export
identify_recipients <- function(y) which(is.na(y))
