#' Identify donors and recipients in a population
#'
#' `identify_donors()` and `identify_recipients()` respectively identify the
#'   respondent (donors) and nonrespondent (recipient) subsets of a population
#'   based on missing responses.
#'
#' @param responses A numeric vector that may contain missing elements (`NA` values).
#'
#' @return
#' * `identify_donors()`: Integer vector of indices where responses are not `NA`.
#' * `identify_recipients()`: Integer vector of indices where responses are `NA`.
#'
#' @export
identify_donors <- function(responses) {

  which(!is.na(responses))

}

#' @rdname identify_donors
#' @export
identify_recipients <- function(responses) {

  which(is.na(responses))

}
