# Data-generating processes for the pseudocensus bias simulation study.
#
# A population of size n is drawn from a superpopulation model
#   y_i = m(x_i) + e_i,   E[e_i | x_i] = 0,
# and a response indicator is drawn with propensity pi(x_i). Units with
# observed responses are donors; units with missing responses are recipients.
#
# The superpopulation estimand is the linear projection parameter
#   beta* = argmin_b E[(y - (1, x')b)^2],
# which coincides with the true coefficients when m is linear and remains a
# well-defined target under misspecification.


# ---- Covariates --------------------------------------------------------------

ar1_correlation <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

#' Draw n rowwise covariate vectors with AR(1) dependence.
#'
#' `gaussian` gives N(0, Sigma) margins; `skewed` pushes the same Gaussian
#' copula through standardised lognormal margins, producing a heavily
#' right-skewed design in which nearest-neighbour geometry is asymmetric.
draw_covariates <- function(n, p, covariate_design, rho = 0.5) {
  z <- matrix(rnorm(n * p), n, p) %*% chol(ar1_correlation(p, rho))

  switch(covariate_design,
    gaussian = z,
    skewed = (exp(z) - exp(0.5)) / sqrt(exp(2) - exp(1)),
    stop("unknown covariate_design: ", covariate_design)
  )
}


# ---- Conditional mean --------------------------------------------------------

#' Coefficients of the linear component, alternating in sign and decaying so
#' that no single covariate dominates: beta_j = 2 (-1)^(j+1) / j.
true_linear_coefficients <- function(p) 2 * (-1)^(seq_len(p) + 1) / seq_len(p)

#' Conditional mean m(x).
#'
#' * `linear`     -- the working linear model is correctly specified.
#' * `quadratic`  -- smooth misspecification: curvature and an interaction.
#' * `rough`      -- non-smooth misspecification: oscillation and a kink,
#'                   the hardest case for local averaging.
conditional_mean <- function(x, mean_spec) {
  linear <- 1 + drop(x %*% true_linear_coefficients(ncol(x)))

  switch(mean_spec,
    linear = linear,
    quadratic = linear + 1.5 * (x[, 1]^2 - 1) + x[, 1] * x[, 2],
    rough = linear + 2 * sin(2 * x[, 1]) + 1.5 * abs(x[, 2]),
    stop("unknown mean_spec: ", mean_spec)
  )
}


# ---- Errors ------------------------------------------------------------------

#' Mean-zero errors, scaled to unit (unconditional) variance in each case so
#' that signal-to-noise is comparable across specifications.
draw_errors <- function(x, error_spec) {
  n <- nrow(x)

  switch(error_spec,
    homoskedastic = rnorm(n),
    heteroskedastic = (0.5 + abs(x[, 1])) * rnorm(n),
    heavy_tailed = rt(n, df = 3) / sqrt(3),
    stop("unknown error_spec: ", error_spec)
  )
}


# ---- Missingness -------------------------------------------------------------

#' Response propensity pi(x) = P(response observed | x).
#'
#' * `mcar`  -- constant propensity; donors and recipients share one law.
#' * `mar`   -- logistic propensity in (x1, x2) with the intercept calibrated
#'              so the population response rate matches `response_rate`.
#' * `shift` -- near-complete response below the 0.6 quantile of x1 and sparse
#'              response above it: recipients concentrate exactly where donors
#'              are scarce, stressing the extrapolation behaviour of
#'              nearest-neighbour imputation.
response_propensity <- function(x, missing_spec, response_rate) {
  n <- nrow(x)

  switch(missing_spec,
    mcar = rep(response_rate, n),
    mar = {
      eta <- 1.5 * x[, 1] - x[, 2]
      intercept <- uniroot(
        function(a) mean(plogis(a + eta)) - response_rate,
        interval = c(-30, 30),
        tol = 1e-8
      )$root
      plogis(intercept + eta)
    },
    shift = {
      low <- x[, 1] <= quantile(x[, 1], 0.6)
      p_high <- (response_rate - 0.6 * 0.95) / 0.4
      ifelse(low, 0.95, p_high)
    },
    stop("unknown missing_spec: ", missing_spec)
  )
}


# ---- Population --------------------------------------------------------------

generate_population <- function(n,
                                p,
                                mean_spec,
                                error_spec,
                                missing_spec,
                                covariate_design,
                                rho = 0.5,
                                response_rate = 0.7) {
  x <- draw_covariates(n, p, covariate_design, rho)
  y <- conditional_mean(x, mean_spec) + draw_errors(x, error_spec)
  observed <- runif(n) < response_propensity(x, missing_spec, response_rate)

  list(x = x, y = y, observed = observed)
}


# ---- Estimand ----------------------------------------------------------------

#' Superpopulation projection parameter beta*, approximated by OLS on a large
#' noiseless draw. Errors are conditionally mean-zero, so projecting m(X)
#' identifies the same parameter as projecting Y, without error-induced Monte
#' Carlo noise. With n_mc = 2e6 the approximation error is O(1e-3) at most and
#' negligible relative to the Monte Carlo error of the study.
projection_target <- function(mean_spec,
                              covariate_design,
                              p,
                              rho = 0.5,
                              n_mc = 2e6,
                              seed = 20260702) {
  set.seed(seed)
  x <- draw_covariates(n_mc, p, covariate_design, rho)
  m <- conditional_mean(x, mean_spec)
  stats::.lm.fit(cbind(1, x), m)$coefficients
}
