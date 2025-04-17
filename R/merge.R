#' Merge imputed and observed responses and locations
#'
#' Combine observed and imputed responses (or locations) into a single array.
#'
#' @inheritParams impute_response
#'
#' @return
#' * `merge_responses()`: returns `responses`, except where the responses for
#'   the recipient indices have been replaced with their imputed equivalent.

#' * `merge_locations()`: returns `locations`, except where the rows for the
#'   recipient indices have been replaced with their imputed equivalent.
#'
#' @export
merge_responses <- function(
    locations,
    responses,
    weights,
    p = 2,
    all_recipients = NULL,
    donors = NULL
) {

  if (is.null(all_recipients)) {

    all_recipients <- identify_recipients(responses)

  }

  if (is.null(donors)) { donors <- identify_donors(responses) }

  imputed_responses <- impute_all_recipient_responses(
    all_recipients = all_recipients,
    locations = locations,
    responses = responses,
    weights = weights,
    donors = donors,
    p = p
  )

  responses[all_recipients] <- imputed_responses

  responses

}




#' @rdname merge_responses
#' @export
merge_locations <- function(
    locations,
    weights,
    p = 2,
    all_recipients = NULL,
    donors = NULL,
    responses = NULL
) {

  if (is.null(all_recipients)) {

    all_recipients <- identify_recipients(responses)

  }

  if (is.null(donors)) { donors <- identify_donors(responses) }

  imputed_recipient_locations <- impute_all_recipient_locations(
    all_recipients = all_recipients,
    donors = donors,
    locations = locations,
    weights = weights,
    p = p
  )

  locations[all_recipients, ] <- imputed_recipient_locations

  locations

}
